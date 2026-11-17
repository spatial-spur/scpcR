//! Whole-operation exact SCPC inference pipeline.

use crate::basis::{scpc_basis, ScpcBasis};
use crate::calibration::{correlation_floor, solve_c0_from_avc};
use crate::covariance;
use crate::distance::{distance_matrix, lower_triangle, Metric};
use crate::linalg::{self, EigenWorkspace};
use crate::quadrature::RULE;
use crate::student_t::{ln_gamma, student_t_ppf};
use crate::threads::{is_serial, run_with_threads};
use crate::{checked_len, CoreError};
use rayon::prelude::*;
use std::f64::consts::{FRAC_2_PI, FRAC_PI_2};
use std::sync::OnceLock;

const MIN_AVC: f64 = 0.00001;
const CGRID_FAC: f64 = 1.2;
const DEFAULT_CV_TOL: f64 = 0.001;
const EXTRA_CV_LEVELS: [f64; 3] = [0.32, 0.10, 0.01];
const SCALE_REL_TOL: f64 = 1e-12;

/// Complete input bundle for one exact SCPC inference call.
///
/// All two-dimensional arrays are contiguous and row-major.
pub struct ExactScpcRequest<'a> {
    /// Coordinate matrix (`n x d`). For geodesic coordinates, select
    /// `Metric::GreatCircleFractionLonLat` or `Metric::GreatCircleFractionLatLon`.
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    pub metric: Metric,
    pub avc: f64,
    pub uncond: bool,
    pub cvs: bool,
    /// Estimated regression coefficients (`k_use`).
    pub coef: &'a [f64],
    /// Coefficient influence matrix (`n x k_use`), column `j` is `n * (B^{-1}[j,.] @ S') + coef[j]`.
    pub wj: &'a [f64],
    /// Conditional regressor direction (`n x k_use`), column `j` is `n * (B^{-1}[j,.] @ X_cond')`.
    /// Required when `uncond == false`.
    pub xj: Option<&'a [f64]>,
    /// Sign/scaling direction (`n x k_use`), for unclustered SCPC `sign(xj)`.
    /// Required when `uncond == false`.
    pub xjs: Option<&'a [f64]>,
    /// Left factor `U` (`n x proj_rank`) of conditional projection `A -> A - U (V' A)`.
    /// Required when `uncond == false` and `proj_rank > 0`.
    pub proj_u: Option<&'a [f64]>,
    /// Right factor `V` (`n x proj_rank`) of conditional projection `A -> A - U (V' A)`.
    /// Required when `uncond == false` and `proj_rank > 0`.
    pub proj_v: Option<&'a [f64]>,
    pub proj_rank: usize,
    /// Critical value solver tolerance in ordinary t-units (released default `0.001`).
    pub cv_tolerance: f64,
    /// Optional thread budget (`Some(1)` runs sequentially on caller thread).
    pub nthreads: Option<usize>,
}

/// Complete output bundle returned by [`run_exact_scpc`].
#[derive(Debug, Clone)]
pub struct ExactScpcOutput {
    /// Row-major `k_use x 6`: `[coef, se, tstat, pvalue, ci_lower, ci_upper]`.
    pub scpcstats: Vec<f64>,
    /// Optional row-major `k_use x 4`: two-sided critical values at `[32%, 10%, 5%, 1%]`.
    pub scpccvs: Option<Vec<f64>>,
    /// Selected spatial basis `n x (q + 1)` in row-major order.
    pub wfin: Vec<f64>,
    pub c0: f64,
    pub cmax: f64,
    pub cvfin: f64,
    pub q: usize,
    pub tied: bool,
}

fn map_indices<T: Send>(
    len: usize,
    f: impl Fn(usize) -> Result<T, CoreError> + Send + Sync,
) -> Result<Vec<T>, CoreError> {
    if is_serial() || len <= 1 {
        (0..len).map(f).collect()
    } else {
        (0..len).into_par_iter().map(f).collect()
    }
}

fn initial_qmax(avc: f64) -> usize {
    if avc >= 0.05 {
        10
    } else if avc >= 0.01 {
        20
    } else if avc >= 0.005 {
        60
    } else {
        120
    }
}

fn covariance_grids(
    dist: &[f64],
    n: usize,
    weights: &[&[f64]],
    widths: &[usize],
    scales: &[f64],
) -> Result<Vec<Vec<f64>>, CoreError> {
    if dist.len() != checked_len(&[n, n])? {
        return Err(CoreError::invalid_input(
            "distance matrix dimensions do not match",
        ));
    }
    covariance::project(n, weights, widths, scales, |i, j, ni, nj, tile| {
        for (row, values) in tile.chunks_exact_mut(nj).take(ni).enumerate() {
            let start = (i + row) * n + j;
            values.copy_from_slice(&dist[start..start + nj]);
        }
    })
}

fn log_one_plus_exp(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

fn angular_product(cos2: f64, ratios: &[f64], logs: &[f64]) -> f64 {
    if cos2 <= 0.0 {
        return 0.0;
    }
    let product: f64 = ratios.iter().map(|value| cos2 / (cos2 + value)).product();
    if product >= f64::MIN_POSITIVE {
        return product.sqrt();
    }
    let log_cos = cos2.ln();
    (-0.5
        * logs
            .iter()
            .map(|ratio| log_one_plus_exp(ratio - log_cos))
            .sum::<f64>())
    .exp()
}

fn complement_integrand(sin2: f64, denom: f64, ratios: &[f64], scale: f64) -> f64 {
    let mut product = 1.0;
    let mut deficit = 0.0;
    for &ratio in ratios {
        let extra = (ratio / scale) * denom;
        let inv = 1.0 / (sin2 + extra);
        deficit += product * extra * inv;
        product *= sin2 * inv;
    }
    deficit / ((1.0 + product.sqrt()) * denom)
}

fn panel(f: &impl Fn(f64, f64) -> f64, left: f64, right: f64) -> f64 {
    const HALF: f64 = FRAC_PI_2 * 0.5;
    type Rule = [(f64, f64, f64); RULE.len()];
    static COMMON: OnceLock<[Rule; 3]> = OnceLock::new();
    let index = match (left, right) {
        (0.0, FRAC_PI_2) => Some(0),
        (0.0, HALF) => Some(1),
        (HALF, FRAC_PI_2) => Some(2),
        _ => None,
    };
    if let Some(index) = index {
        let rules = COMMON.get_or_init(|| {
            [(0.0, FRAC_PI_2), (0.0, HALF), (HALF, FRAC_PI_2)].map(|(lo, hi)| {
                RULE.map(|(node, weight)| {
                    let (sin, cos) = (lo + (hi - lo) * node).sin_cos();
                    (sin * sin, cos * cos, (hi - lo) * weight)
                })
            })
        });
        return rules[index]
            .iter()
            .map(|&(sin2, cos2, weight)| weight * f(sin2, cos2))
            .sum();
    }
    let width = right - left;
    width
        * RULE
            .iter()
            .map(|&(node, weight)| {
                let (sin, cos) = (left + width * node).sin_cos();
                weight * f(sin * sin, cos * cos)
            })
            .sum::<f64>()
}

fn adaptive(
    f: &impl Fn(f64, f64) -> f64,
    lo: f64,
    hi: f64,
    whole: f64,
    tol: f64,
    depth: u8,
) -> f64 {
    let mid = 0.5 * (lo + hi);
    let left = panel(f, lo, mid);
    let right = panel(f, mid, hi);
    if depth == 20 || (left + right - whole).abs() <= tol.max(1e-11 * whole.abs()) {
        return left + right;
    }
    adaptive(f, lo, mid, left, tol * 0.5, depth + 1)
        + adaptive(f, mid, hi, right, tol * 0.5, depth + 1)
}

fn integrate(f: &impl Fn(f64, f64) -> f64) -> f64 {
    let whole = panel(f, 0.0, FRAC_PI_2);
    let tolerance = if whole.abs() < 1e-8 {
        1e-10 * whole.abs()
    } else {
        1e-11
    };
    adaptive(f, 0.0, FRAC_PI_2, whole, tolerance, 0) * FRAC_2_PI
}

struct ProbabilityWorkspace {
    eigensolver: EigenWorkspace,
    ratios: Vec<f64>,
    logs: Vec<f64>,
}

impl ProbabilityWorkspace {
    fn new(q: usize) -> Self {
        Self {
            eigensolver: EigenWorkspace::new(q + 1),
            ratios: Vec::with_capacity(q),
            logs: Vec::with_capacity(q),
        }
    }

    /// Evaluate tail-stable adaptive angular/complement rejection probability.
    fn evaluate(&mut self, omega: &[f64], width: usize, cv: f64) -> Result<f64, CoreError> {
        if !cv.is_finite() {
            return Err(CoreError::invalid_input("SCPC critical value must be finite"));
        }
        if cv <= 0.0 {
            return Ok(1.0);
        }
        let (positive, negative) = if cv > 1.0 {
            (cv.recip().powi(2), -1.0)
        } else {
            (1.0, -cv.powi(2))
        };
        let eigenvalues = self
            .eigensolver
            .weighted_eigenvalues(omega, width, positive, negative)?;
        let max_positive = eigenvalues.iter().copied().fold(0.0, f64::max);
        if !max_positive.is_finite() || max_positive <= 0.0 {
            return Ok(0.0);
        }
        self.ratios.clear();
        self.logs.clear();
        let max_pos_ln = max_positive.ln();
        for &value in eigenvalues.iter().filter(|&&value| value < 0.0) {
            self.ratios.push(-value / max_positive);
            self.logs.push((-value).ln() - max_pos_ln);
        }
        if self.ratios.is_empty() {
            return Ok(1.0);
        }
        if self.ratios.len() == 1 {
            return Ok(FRAC_2_PI * (-0.5 * self.logs[0]).exp().atan());
        }

        let scale: f64 = self.ratios.iter().sum();
        if scale == 0.0 {
            return Ok(1.0);
        }
        let probability = if scale < 0.0625 {
            let deficit = integrate(&|sin2, cos2| {
                complement_integrand(sin2, cos2 + scale * sin2, &self.ratios, scale)
            });
            1.0 - scale.sqrt() * deficit
        } else {
            integrate(&|_, cos2| angular_product(cos2, &self.ratios, &self.logs))
        };
        Ok(probability.clamp(0.0, 1.0))
    }

    fn maximum(
        &mut self,
        omegas: &[f64],
        width: usize,
        cv: f64,
        known: Option<(f64, usize)>,
    ) -> Result<(f64, usize), CoreError> {
        let mut best = (-1.0, 0);
        for (i, omega) in omegas.chunks_exact(width * width).enumerate() {
            let p = if let Some((p, known_i)) = known.filter(|&(_, index)| index == i) {
                debug_assert_eq!(i, known_i);
                p
            } else {
                self.evaluate(omega, width, cv)?
            };
            if p > best.0 {
                best = (p, i);
            }
        }
        Ok(best)
    }

    /// Released `get_cv(Oms, q, level)` solver: steps out by `1/sqrt(q)` and
    /// bisects to `tol / sqrt(q)`, returning the upper bracket endpoint in `t` units.
    fn critical(
        &mut self,
        omegas: &[f64],
        width: usize,
        q: usize,
        level: f64,
        tolerance: f64,
    ) -> Result<f64, CoreError> {
        let root_q = (q as f64).sqrt();
        let step = 1.0 / root_q;
        let tol = tolerance * step;
        let matrix_size = width * width;
        let grid_count = omegas.len() / matrix_size;

        let mut rp = 1.0;
        let mut active = 0usize;
        let mut cv0 = student_t_ppf(1.0 - 0.5 * level, q)? / root_q;
        let mut cv1 = cv0;

        for _ in 0..=grid_count + 4 {
            if rp <= level {
                break;
            }
            let omega = &omegas[active * matrix_size..(active + 1) * matrix_size];
            cv1 = cv0;
            for iter in 0..1024 {
                if self.evaluate(omega, width, cv1)? > level {
                    cv0 = cv1;
                    cv1 += if iter < 32 { step } else { cv1 };
                } else {
                    break;
                }
            }
            if self.evaluate(omega, width, cv1)? > level {
                return Err(CoreError::calibration(
                    "SCPC critical value could not be bracketed",
                ));
            }

            let mut p_upper = self.evaluate(omega, width, cv1)?;
            while cv1 - cv0 > tol {
                let mid = 0.5 * (cv0 + cv1);
                if mid == cv0 || mid == cv1 {
                    break;
                }
                let p_mid = self.evaluate(omega, width, mid)?;
                if p_mid > level {
                    cv0 = mid;
                } else {
                    cv1 = mid;
                    p_upper = p_mid;
                }
            }

            let (max_p, max_i) = self.maximum(omegas, width, cv1, Some((p_upper, active)))?;
            if max_i == active {
                break;
            }
            active = max_i;
            rp = max_p;
            cv0 = cv1;
        }
        Ok(cv1 * root_q)
    }
}

fn select_final_q(
    oms: &[f64],
    width: usize,
    qmax: usize,
    tolerance: f64,
) -> Result<(usize, f64), CoreError> {
    let results = map_indices(qmax, |idx| {
        let q = idx + 1;
        let mut ws = ProbabilityWorkspace::new(q);
        let cv = ws.critical(oms, width, q, 0.05, tolerance)?;
        let qf = q as f64;
        let gamma_ratio = (ln_gamma(0.5 * (qf + 1.0)) - ln_gamma(0.5 * qf)).exp();
        let length = cv * gamma_ratio / qf.sqrt();
        Ok((q, cv, length))
    })?;
    let best = results
        .into_iter()
        .min_by(|a, b| a.2.total_cmp(&b.2).then_with(|| a.0.cmp(&b.0)))
        .ok_or(CoreError::invalid_input("empty qmax search range"))?;
    Ok((best.0, best.1))
}

fn build_conditional_basis(
    wfin: &[f64],
    n: usize,
    q: usize,
    xj: &[f64],
    xjs: &[f64],
    k_use: usize,
    col_j: usize,
    proj_u: Option<&[f64]>,
    proj_v: Option<&[f64]>,
    proj_rank: usize,
) -> Result<Vec<f64>, CoreError> {
    let width = q + 1;
    let mut wx = vec![0.0; checked_len(&[n, width])?];
    for i in 0..n {
        let xj_i = xj[i * k_use + col_j];
        let xjs_i = xjs[i * k_use + col_j];
        wx[i * width] = wfin[i * width] * xj_i * xjs_i;
    }

    let mut pre_max = 0.0_f64;
    if q > 0 {
        // Denominator weights are S_x M_V D_x B_q.
        let mut rx = vec![0.0; checked_len(&[n, q])?];
        for c in 0..q {
            let mut col_sq = 0.0_f64;
            for i in 0..n {
                let val = wfin[i * width + 1 + c] * xj[i * k_use + col_j];
                rx[i * q + c] = val;
                col_sq += val * val;
            }
            pre_max = pre_max.max(col_sq.sqrt());
        }
        let rr = if proj_rank > 0 {
            let u = proj_u.ok_or(CoreError::invalid_input(
                "missing conditional projection U factor",
            ))?;
            let v = proj_v.ok_or(CoreError::invalid_input(
                "missing conditional projection V factor",
            ))?;
            let mut vt_rx = vec![0.0; checked_len(&[proj_rank, q])?];
            let mut u_vt_rx = vec![0.0; checked_len(&[n, q])?];
            linalg::multiply(v, &rx, &mut vt_rx, proj_rank, n, q, true);
            linalg::multiply(u, &vt_rx, &mut u_vt_rx, n, proj_rank, q, false);
            for (r, p) in rx.iter_mut().zip(u_vt_rx.iter()) {
                *r -= *p;
            }
            rx
        } else {
            rx
        };
        for i in 0..n {
            let xjs_i = xjs[i * k_use + col_j];
            for c in 0..q {
                wx[i * width + 1 + c] = rr[i * q + c] * xjs_i;
            }
        }
    }

    if wx.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::degenerate(
            "conditional projection produced non-finite W values",
        ));
    }
    if width > 1 {
        let tol = SCALE_REL_TOL * pre_max;
        let all_collapsed = !pre_max.is_finite()
            || pre_max == 0.0
            || (1..width).all(|c| {
                let col_norm = (0..n)
                    .map(|i| wx[i * width + c].powi(2))
                    .sum::<f64>()
                    .sqrt();
                !col_norm.is_finite() || col_norm <= tol
            });
        if all_collapsed {
            return Err(CoreError::degenerate(
                "conditional projection collapsed all non-constant W columns",
            ));
        }
    }
    Ok(wx)
}

fn validate_request(req: &ExactScpcRequest<'_>) -> Result<f64, CoreError> {
    if req.n < 2 || req.d == 0 || req.coef.is_empty() {
        return Err(CoreError::invalid_input(
            "invalid SCPC sample or coefficient dimensions",
        ));
    }
    let k_use = req.coef.len();
    if req.coords.len() != checked_len(&[req.n, req.d])?
        || req.wj.len() != checked_len(&[req.n, k_use])?
    {
        return Err(CoreError::invalid_input(
            "SCPC input array dimensions do not match",
        ));
    }
    if req.coef.iter().any(|v| !v.is_finite()) || req.wj.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input("SCPC inputs must be finite"));
    }
    if !req.avc.is_finite() || req.avc <= 0.001 || req.avc >= 0.99 {
        return Err(CoreError::invalid_input(
            "avc must lie strictly in (0.001, 0.99)",
        ));
    }
    let tol = if req.cv_tolerance > 0.0 && req.cv_tolerance.is_finite() {
        req.cv_tolerance
    } else {
        DEFAULT_CV_TOL
    };
    if !req.uncond {
        let xj = req.xj.ok_or(CoreError::invalid_input(
            "conditional SCPC requires xj direction matrix",
        ))?;
        let xjs = req.xjs.ok_or(CoreError::invalid_input(
            "conditional SCPC requires xjs scaling matrix",
        ))?;
        if xj.len() != checked_len(&[req.n, k_use])? || xjs.len() != checked_len(&[req.n, k_use])? {
            return Err(CoreError::invalid_input(
                "SCPC conditional xj/xjs dimensions do not match",
            ));
        }
        if xj.iter().any(|v| !v.is_finite()) || xjs.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input(
                "SCPC conditional directions must be finite",
            ));
        }
        if req.proj_rank > 0 {
            let expected = checked_len(&[req.n, req.proj_rank])?;
            let u = req.proj_u.ok_or(CoreError::invalid_input(
                "conditional SCPC requires projection U factor",
            ))?;
            let v = req.proj_v.ok_or(CoreError::invalid_input(
                "conditional SCPC requires projection V factor",
            ))?;
            if u.len() != expected || v.len() != expected {
                return Err(CoreError::invalid_input(
                    "conditional projection factor dimensions do not match",
                ));
            }
            if u.iter().any(|val| !val.is_finite()) || v.iter().any(|val| !val.is_finite()) {
                return Err(CoreError::invalid_input(
                    "conditional projection factors must be finite",
                ));
            }
        }
    }
    Ok(tol)
}

/// Run the full exact SCPC inference pipeline in one call.
pub fn run_exact_scpc(req: &ExactScpcRequest<'_>) -> Result<ExactScpcOutput, CoreError> {
    let cv_tolerance = validate_request(req)?;
    let n = req.n;
    let k_use = req.coef.len();

    run_with_threads(req.nthreads, || {
        let distmat = distance_matrix(req.coords, n, req.d, req.metric)?;
        let (c0, cmax) = {
            let dist_vec = lower_triangle(&distmat, n);
            let floor = correlation_floor(&dist_vec)?;
            let min_avc_eff = floor + (1.0 - floor) * MIN_AVC;
            let c0 = solve_c0_from_avc(&dist_vec, req.avc)?;
            let cmax = solve_c0_from_avc(&dist_vec, min_avc_eff)?;
            (c0, cmax)
        };

        let nc = ((cmax / c0).ln() / CGRID_FAC.ln()).ceil().max(2.0) as usize;
        let mut scales = Vec::with_capacity(nc);
        let mut c = c0;
        for _ in 0..nc {
            scales.push(c);
            c *= CGRID_FAC;
        }

        let mut qmax = initial_qmax(req.avc);
        let (basis, oms_full, q, cvfin): (ScpcBasis, Vec<f64>, usize, f64) = loop {
            qmax = qmax.min(n - 1);
            let basis = scpc_basis(&distmat, n, c0, qmax)?;
            let mut grids = covariance_grids(&distmat, n, &[&basis.weights], &[qmax + 1], &scales)?;
            let oms = grids.remove(0);
            let (q_opt, cv_opt) = select_final_q(&oms, qmax + 1, qmax, cv_tolerance)?;
            if q_opt < qmax || qmax == n - 1 {
                break (basis, oms, q_opt, cv_opt);
            }
            qmax = ((qmax as f64) + (qmax as f64) * 0.5).round() as usize;
        };

        let full_width = qmax + 1;
        let width = q + 1;
        let mut wfin = vec![0.0; checked_len(&[n, width])?];
        for (dst, src) in wfin
            .chunks_exact_mut(width)
            .zip(basis.weights.chunks_exact(full_width))
        {
            dst.copy_from_slice(&src[..width]);
        }

        // Unconditional score statistics: Wfin' * wj
        let mut projected = vec![0.0; checked_len(&[width, k_use])?];
        linalg::multiply(&wfin, req.wj, &mut projected, width, n, k_use, true);
        let q_root = (q as f64).sqrt();
        let se_denom = q_root * (n as f64).sqrt();

        let mut se_vec = Vec::with_capacity(k_use);
        let mut tau_vec = Vec::with_capacity(k_use);
        for j in 0..k_use {
            let wj_norm = (0..n)
                .map(|i| req.wj[i * k_use + j].powi(2))
                .sum::<f64>()
                .sqrt();
            let norm = (1..width)
                .map(|row| projected[row * k_use + j].powi(2))
                .sum::<f64>()
                .sqrt();
            if !norm.is_finite()
                || !wj_norm.is_finite()
                || wj_norm == 0.0
                || norm <= SCALE_REL_TOL * wj_norm
            {
                return Err(CoreError::degenerate(
                    "SCPC variance projection is degenerate",
                ));
            }
            tau_vec.push(q_root * projected[j] / norm);
            se_vec.push(norm / se_denom);
        }

        // Precompute conditional Oms before dropping the n x n distance matrix.
        let omsx_list = if req.uncond {
            drop(distmat);
            None
        } else {
            let xj = req.xj.expect("validated conditional xj");
            let xjs = req.xjs.expect("validated conditional xjs");
            let wx_list = map_indices(k_use, |j| {
                build_conditional_basis(
                    &wfin,
                    n,
                    q,
                    xj,
                    xjs,
                    k_use,
                    j,
                    req.proj_u,
                    req.proj_v,
                    req.proj_rank,
                )
            })?;
            let wx_refs: Vec<&[f64]> = wx_list.iter().map(Vec::as_slice).collect();
            let widths = vec![width; k_use];
            let grids = covariance_grids(&distmat, n, &wx_refs, &widths, &scales)?;
            drop(distmat);
            Some(grids)
        };

        // Unconditional rejection probabilities (and optional cvs_uncond reusing cvfin at 0.05)
        let p_u = map_indices(k_use, |j| {
            let mut ws = ProbabilityWorkspace::new(q);
            ws.maximum(&oms_full, full_width, tau_vec[j].abs() / q_root, None)
                .map(|(p, _)| p)
        })?;

        let cvs_uncond: Option<[f64; 4]> = if req.cvs {
            let vals = map_indices(EXTRA_CV_LEVELS.len(), |idx| {
                let mut ws = ProbabilityWorkspace::new(q);
                ws.critical(&oms_full, full_width, q, EXTRA_CV_LEVELS[idx], cv_tolerance)
            })?;
            Some([vals[0], vals[1], cvfin, vals[2]])
        } else {
            None
        };

        let mut scpcstats = vec![0.0; checked_len(&[k_use, 6])?];
        let mut scpccvs = if req.cvs {
            Some(vec![0.0; checked_len(&[k_use, 4])?])
        } else {
            None
        };

        if let Some(omsx_list) = omsx_list {
            let cond_results = map_indices(k_use, |j| {
                let omsx = &omsx_list[j];
                let mut ws = ProbabilityWorkspace::new(q);
                let (p_c, _) = ws.maximum(omsx, width, tau_vec[j].abs() / q_root, None)?;
                let cvx = ws.critical(omsx, width, q, 0.05, cv_tolerance)?;
                let cond_cvs = if req.cvs {
                    let extra = map_indices(EXTRA_CV_LEVELS.len(), |idx| {
                        let mut inner = ProbabilityWorkspace::new(q);
                        inner.critical(omsx, width, q, EXTRA_CV_LEVELS[idx], cv_tolerance)
                    })?;
                    Some([extra[0], extra[1], cvx, extra[2]])
                } else {
                    None
                };
                Ok((p_c, cvx, cond_cvs))
            })?;

            for j in 0..k_use {
                let (p_c, cvx, cond_cvs) = cond_results[j];
                let p_final = p_u[j].max(p_c);
                let cv = cvfin.max(cvx);
                let row = &mut scpcstats[j * 6..(j + 1) * 6];
                row[0] = req.coef[j];
                row[1] = se_vec[j];
                row[2] = tau_vec[j];
                row[3] = p_final;
                row[4] = req.coef[j] - cv * se_vec[j];
                row[5] = req.coef[j] + cv * se_vec[j];

                if let (Some(out_cvs), Some(uncond_cvs), Some(c_cvs)) = (
                    scpccvs.as_mut(),
                    cvs_uncond.as_ref(),
                    cond_cvs.as_ref(),
                ) {
                    let row = &mut out_cvs[j * 4..(j + 1) * 4];
                    for idx in 0..4 {
                        row[idx] = uncond_cvs[idx].max(c_cvs[idx]);
                    }
                }
            }
        } else {
            for j in 0..k_use {
                let row = &mut scpcstats[j * 6..(j + 1) * 6];
                row[0] = req.coef[j];
                row[1] = se_vec[j];
                row[2] = tau_vec[j];
                row[3] = p_u[j];
                row[4] = req.coef[j] - cvfin * se_vec[j];
                row[5] = req.coef[j] + cvfin * se_vec[j];
            }
            if let (Some(out_cvs), Some(uncond_cvs)) = (scpccvs.as_mut(), cvs_uncond.as_ref()) {
                for row in out_cvs.chunks_exact_mut(4) {
                    row.copy_from_slice(uncond_cvs);
                }
            }
        }

        Ok(ExactScpcOutput {
            scpcstats,
            scpccvs,
            wfin,
            c0,
            cmax,
            cvfin,
            q,
            tied: basis.tied,
        })
    })?
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::student_t::betainc;
    use crate::ErrorKind;

    fn solve_linear(a: &mut [f64], b: &mut [f64], p: usize) {
        for col in 0..p {
            let mut pivot = col;
            let mut max_val = a[col * p + col].abs();
            for row in (col + 1)..p {
                if a[row * p + col].abs() > max_val {
                    max_val = a[row * p + col].abs();
                    pivot = row;
                }
            }
            if pivot != col {
                for k in 0..p {
                    a.swap(col * p + k, pivot * p + k);
                }
                b.swap(col, pivot);
            }
            let diag = a[col * p + col];
            for k in col..p {
                a[col * p + k] /= diag;
            }
            b[col] /= diag;
            for row in 0..p {
                if row == col {
                    continue;
                }
                let factor = a[row * p + col];
                for k in col..p {
                    a[row * p + k] -= factor * a[col * p + k];
                }
                b[row] -= factor * b[col];
            }
        }
    }

    fn invert_symmetric(mat: &[f64], p: usize) -> Vec<f64> {
        let mut inv = vec![0.0; p * p];
        for col in 0..p {
            let mut a = mat.to_vec();
            let mut rhs = vec![0.0; p];
            rhs[col] = 1.0;
            solve_linear(&mut a, &mut rhs, p);
            for row in 0..p {
                inv[row * p + col] = rhs[row];
            }
        }
        inv
    }

    fn modified_gram_schmidt(x: &[f64], n: usize, p: usize) -> Vec<f64> {
        let mut q = vec![0.0; n * p];
        let mut rank = 0usize;
        for col in 0..p {
            let mut v: Vec<f64> = (0..n).map(|i| x[i * p + col]).collect();
            let orig_norm = v.iter().map(|val| val * val).sum::<f64>().sqrt();
            if orig_norm == 0.0 {
                continue;
            }
            for val in &mut v {
                *val /= orig_norm;
            }
            for _ in 0..2 {
                for r in 0..rank {
                    let dot: f64 = (0..n).map(|i| q[i * p + r] * v[i]).sum();
                    for i in 0..n {
                        v[i] -= dot * q[i * p + r];
                    }
                }
            }
            let norm = v.iter().map(|val| val * val).sum::<f64>().sqrt();
            if norm > 1e-12 {
                for i in 0..n {
                    q[i * p + rank] = v[i] / norm;
                }
                rank += 1;
            }
        }
        if rank < p {
            let mut compact = vec![0.0; n * rank];
            for i in 0..n {
                compact[i * rank..(i + 1) * rank].copy_from_slice(&q[i * p..i * p + rank]);
            }
            compact
        } else {
            q
        }
    }

    fn load_bench_csv(name: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let path = format!(
            "{}/../../../verification/data/{}",
            env!("CARGO_MANIFEST_DIR"),
            name
        );
        let text = std::fs::read_to_string(path).expect("fixture csv");
        let mut y = Vec::new();
        let mut x = Vec::new();
        let mut coords = Vec::new();
        for line in text.lines().skip(1) {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<f64> = line
                .split(',')
                .map(|s| s.trim().parse::<f64>().unwrap())
                .collect();
            y.push(cols[1]);
            x.extend_from_slice(&[1.0, cols[2], cols[3]]);
            coords.extend_from_slice(&[cols[6], cols[7]]);
        }
        (y, x, coords)
    }

    fn prepare_bench_ols(y: &[f64], x: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, usize) {
        let n = y.len();
        let p = 3usize;
        let mut xtx = vec![0.0; p * p];
        let mut xty = vec![0.0; p];
        for i in 0..n {
            for r in 0..p {
                xty[r] += x[i * p + r] * y[i];
                for c in 0..p {
                    xtx[r * p + c] += x[i * p + r] * x[i * p + c];
                }
            }
        }
        let inv_xtx = invert_symmetric(&xtx, p);
        let mut coef = xty.clone();
        let mut a = xtx.clone();
        solve_linear(&mut a, &mut coef, p);

        let mut wj = vec![0.0; n * p];
        let mut xj = vec![0.0; n * p];
        let mut xjs = vec![0.0; n * p];
        for i in 0..n {
            let yhat: f64 = (0..p).map(|k| x[i * p + k] * coef[k]).sum();
            let resid = y[i] - yhat;
            for j in 0..p {
                let x_dir: f64 =
                    (0..p).map(|k| (n as f64) * inv_xtx[j * p + k] * x[i * p + k]).sum();
                xj[i * p + j] = x_dir;
                xjs[i * p + j] = if x_dir > 0.0 {
                    1.0
                } else if x_dir < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                wj[i * p + j] = x_dir * resid + coef[j];
            }
        }
        let q_mat = modified_gram_schmidt(x, n, p);
        let proj_rank = q_mat.len() / n;
        (coef, wj, xj, xjs, q_mat, proj_rank)
    }

    #[test]
    fn conditional_basis_weights_apply_regressor_before_projection() {
        let x = [-2.0, -1.0, 1.0, 2.0];
        let sign = [-1.0, -1.0, 1.0, 1.0];
        let norm = 10.0_f64.sqrt();
        let w: Vec<_> = x.iter().flat_map(|value| [0.5, -value / norm]).collect();
        let projection: Vec<_> = x.iter().flat_map(|value| [0.5, value / norm]).collect();
        let weights = build_conditional_basis(
            &w,
            4,
            1,
            &x,
            &sign,
            1,
            0,
            Some(&projection),
            Some(&projection),
            2,
        )
        .unwrap();
        let a = 1.5 / norm;
        let expected = [1.0, a, 0.5, -a, 0.5, a, 1.0, -a];
        for (got, expected) in weights.iter().zip(expected) {
            assert!((got - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn identity_covariance_probability_matches_student_t_and_near_zero_tail() {
        for q in [1usize, 2, 5, 10] {
            let width = q + 1;
            let mut identity = vec![0.0; width * width];
            for i in 0..width {
                identity[i * width + i] = 1.0;
            }
            let mut ws = ProbabilityWorkspace::new(q);
            for t in [1e-6, 0.1, 1.0, 2.5, 10.0] {
                let cv_norm = t / (q as f64).sqrt();
                let got = ws.evaluate(&identity, width, cv_norm).unwrap();
                let expected = if t * t < q as f64 {
                    1.0 - betainc(0.5, 0.5 * q as f64, (t * t) / (q as f64 + t * t))
                } else {
                    betainc(0.5 * q as f64, 0.5, (q as f64) / (q as f64 + t * t))
                };
                assert!(
                    (got - expected).abs() < 1e-10,
                    "q={q} t={t}: got={got}, expected={expected}"
                );
            }
        }

        // Near-zero test statistic reproducer from rust_core_handoff.md:
        // coords=[[0.0],[1.0]], coef=[1e-6], wj=[[1+1e-6],[-1+1e-6]], avc=0.03, uncond=true
        let coords = [0.0, 1.0];
        let coef = [1e-6];
        let wj = [1.0 + 1e-6, -1.0 + 1e-6];
        let out = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n: 2,
            d: 1,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: true,
            cvs: false,
            coef: &coef,
            wj: &wj,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out.q, 1);
        let rho = 0.03_f64;
        let expected_p = 1.0 - FRAC_2_PI * (1e-6 * ((1.0 - rho) / (1.0 + rho)).sqrt()).atan();
        assert!(
            (out.scpcstats[3] - expected_p).abs() < 1e-12,
            "got p={}, expected={expected_p}",
            out.scpcstats[3]
        );
    }

    #[test]
    fn exact_scpc_matches_release_score_statistics_and_is_scale_invariant() {
        let (y, x, coords) = load_bench_csv("bench_200.csv");
        let n = y.len();
        let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

        let out_euclid = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();

        assert_eq!(out_euclid.q, 19);
        assert!(
            (out_euclid.c0 - 3.500409719601247).abs() < 1e-3,
            "c0={}",
            out_euclid.c0
        );
        let dist_euclid = distance_matrix(&coords, n, 2, Metric::Euclidean).unwrap();
        let dvec_euclid = lower_triangle(&dist_euclid, n);
        assert!(
            (crate::calibration::avg_correlation(out_euclid.c0, &dvec_euclid) - 0.03).abs() < 1e-12
        );
        assert!(
            (out_euclid.cvfin - 2.1008365544083087).abs() < 1e-3,
            "cvfin={}",
            out_euclid.cvfin
        );
        let expected_stats = [
            [1.4834167407175964, 0.11454779987636342, 12.95019845268713],
            [1.5473280763431876, 0.09650902048915827, 16.032989128897157],
            [-0.6589528089605742, 0.10738682321940547, -6.13625386435211],
        ];
        for (row, expected) in out_euclid.scpcstats.chunks_exact(6).zip(expected_stats) {
            for (got, exp) in row.iter().zip(expected) {
                assert!((got - exp).abs() <= 1e-3 * exp.abs().max(1.0), "got={got}, exp={exp}");
            }
        }
        // Check 5% CV reuse consistency between scpcstats CI and scpccvs[, 2]
        let cvs = out_euclid.scpccvs.as_ref().unwrap();
        for j in 0..3 {
            let se = out_euclid.scpcstats[j * 6 + 1];
            let half_w = out_euclid.scpcstats[j * 6 + 5] - out_euclid.scpcstats[j * 6];
            assert!((half_w - cvs[j * 4 + 2] * se).abs() < 1e-12);
        }

        // Single-threaded serial execution on caller thread matches multi-threaded exactly
        let out_one = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: false,
            coef: &coef,
            wj: &wj,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out_one.q, out_euclid.q);
        for (a, b) in out_one.scpcstats.iter().zip(out_euclid.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-12);
        }

        // Positive outcome rescaling (`y * 1e-12`) preserves q, t, p and scales SE, CI proportionally
        let y_scale = 1e-12;
        let y_scaled: Vec<f64> = y.iter().map(|&v| v * y_scale).collect();
        let (coef_ys, wj_ys, _, _, _, _) = prepare_bench_ols(&y_scaled, &x);
        let out_ys = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef_ys,
            wj: &wj_ys,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap();
        assert_eq!(out_ys.q, out_euclid.q);
        for j in 0..3 {
            assert!((out_ys.scpcstats[j * 6 + 1] - out_euclid.scpcstats[j * 6 + 1] * y_scale).abs() < 1e-22);
            assert!((out_ys.scpcstats[j * 6 + 2] - out_euclid.scpcstats[j * 6 + 2]).abs() < 1e-10);
            assert!((out_ys.scpcstats[j * 6 + 3] - out_euclid.scpcstats[j * 6 + 3]).abs() < 1e-10);
            assert!((out_ys.scpcstats[j * 6 + 4] - out_euclid.scpcstats[j * 6 + 4] * y_scale).abs() < 1e-22);
            assert!((out_ys.scpcstats[j * 6 + 5] - out_euclid.scpcstats[j * 6 + 5] * y_scale).abs() < 1e-22);
        }

        // Coordinate rescaling (`coords * 1e-6`) preserves q, se, t, p and scales c0 by 1e6
        let c_scale = 1e-6;
        let coords_scaled: Vec<f64> = coords.iter().map(|&v| v * c_scale).collect();
        let out_cs = run_exact_scpc(&ExactScpcRequest {
            coords: &coords_scaled,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: false,
            coef: &coef,
            wj: &wj,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_cs.q, out_euclid.q);
        assert!((out_cs.c0 * c_scale - out_euclid.c0).abs() / out_euclid.c0 < 1e-10);
        for (a, b) in out_cs.scpcstats.iter().zip(out_euclid.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-10, "a={a}, b={b}");
        }

        // Exact zero wj produces structured DegenerateVariance error
        let wj_zero = vec![0.0; n * 3];
        let err = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: true,
            cvs: false,
            coef: &coef,
            wj: &wj_zero,
            xj: None,
            xjs: None,
            proj_u: None,
            proj_v: None,
            proj_rank: 0,
            cv_tolerance: 0.001,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err.kind(), ErrorKind::DegenerateVariance);

        // Latlong conditional (`s_1 = lon, s_2 = lat`)
        let out_geo = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::GreatCircleFractionLonLat,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(2),
        })
        .unwrap();
        assert_eq!(out_geo.q, 18);
    }

    #[test]
    fn exact_scpc_matches_release_score_statistics_bench_1000() {
        let (y, x, coords) = load_bench_csv("bench_1000.csv");
        let n = y.len();
        let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

        let out = run_exact_scpc(&ExactScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj,
            xj: Some(&xj),
            xjs: Some(&xjs),
            proj_u: Some(&q_mat),
            proj_v: Some(&q_mat),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(4),
        })
        .unwrap();

        assert_eq!(out.q, 8);
        let dist_1000 = distance_matrix(&coords, n, 2, Metric::Euclidean).unwrap();
        let dvec_1000 = lower_triangle(&dist_1000, n);
        assert!((crate::calibration::avg_correlation(out.c0, &dvec_1000) - 0.03).abs() < 1e-12);
        assert!((out.c0 - 0.43896032977784594).abs() < 1e-3, "c0={}", out.c0);
        assert!(
            (out.cvfin - 2.365574447704166).abs() < 1e-3,
            "cvfin={}",
            out.cvfin
        );
        let expected_stats = [
            [1.5780823405600297, 0.13216765018777865, 11.940004519396016],
            [1.745699446088806, 0.0697326838936868, 25.03416401913148],
            [-0.9039531747843945, 0.037101498337623266, -24.3643306951765],
        ];
        for (row, expected) in out.scpcstats.chunks_exact(6).zip(expected_stats) {
            assert!((row[0] - expected[0]).abs() < 1e-12);
            for (got, exp) in row.iter().zip(expected) {
                assert!((got - exp).abs() <= 2e-3 * exp.abs().max(1.0), "got={got}, exp={exp}");
            }
        }
    }
}
