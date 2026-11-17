//! Shared `extendr` R FFI adapter for `scpcR` and `spuR`.

use extendr_api::prelude::*;
use spatial_core::{run_exact_scpc, CoreError, ErrorKind, ExactScpcRequest, Metric};

fn core_err(err: CoreError) -> Error {
    match err.kind() {
        ErrorKind::InvalidInput => Error::Other(err.message().to_string()),
        ErrorKind::UnsupportedCase => Error::Other(format!("unsupported case: {}", err.message())),
        ErrorKind::TiedEigenspace => Error::Other(format!("tied eigenspace: {}", err.message())),
        ErrorKind::CalibrationFailed => {
            Error::Other(format!("calibration failed: {}", err.message()))
        }
        ErrorKind::DegenerateVariance => {
            Error::Other(format!("degenerate variance: {}", err.message()))
        }
        ErrorKind::InternalFailure => Error::Other(format!("internal failure: {}", err.message())),
    }
}

/// Convert an R column-major matrix buffer (`nrows x ncols`) into a contiguous row-major `Vec<f64>`.
pub fn col_major_to_row_major(
    data: &[f64],
    nrows: usize,
    ncols: usize,
) -> extendr_api::Result<Vec<f64>> {
    let expected = nrows
        .checked_mul(ncols)
        .ok_or_else(|| Error::Other("matrix dimensions overflow usize".into()))?;
    if data.len() != expected {
        return Err(Error::Other(format!(
            "column-major buffer length {} does not match {} x {}",
            data.len(),
            nrows,
            ncols
        )));
    }
    let mut row_major = vec![0.0; expected];
    for col in 0..ncols {
        let src_col = &data[col * nrows..(col + 1) * nrows];
        for row in 0..nrows {
            row_major[row * ncols + col] = src_col[row];
        }
    }
    Ok(row_major)
}

/// Convert a contiguous row-major slice (`nrows x ncols`) into an R column-major `RMatrix<f64>`.
pub fn row_major_to_rmatrix(
    row_major: &[f64],
    nrows: usize,
    ncols: usize,
) -> extendr_api::Result<RMatrix<f64>> {
    let expected = nrows
        .checked_mul(ncols)
        .ok_or_else(|| Error::Other("matrix dimensions overflow usize".into()))?;
    if row_major.len() != expected {
        return Err(Error::Other(format!(
            "row-major buffer length {} does not match {} x {}",
            row_major.len(),
            nrows,
            ncols
        )));
    }
    if i32::try_from(nrows).is_err() || i32::try_from(ncols).is_err() {
        return Err(Error::Other(
            "matrix dimensions exceed R 32-bit dim limit".into(),
        ));
    }
    let mut mat = RMatrix::<f64>::new(nrows, ncols);
    let dst = mat.data_mut();
    for col in 0..ncols {
        let dst_col = &mut dst[col * nrows..(col + 1) * nrows];
        for row in 0..nrows {
            dst_col[row] = row_major[row * ncols + col];
        }
    }
    Ok(mat)
}

/// Bridge one whole exact SCPC inference call from R `SEXP` wrappers to `spatial_core::run_exact_scpc`.
///
/// All R column-major matrices are copied into row-major Rust buffers on the R main thread
/// prior to entering Rayon-parallel numerical routines, and R output matrices are allocated
/// back on the R main thread after `run_exact_scpc` completes.
#[allow(clippy::too_many_arguments)]
pub fn run_exact_scpc_r(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    let n = coords.nrows();
    let d = coords.ncols();
    let k_use = coef.len();

    if wj.nrows() != n || wj.ncols() != k_use {
        return Err(Error::Other(format!(
            "wj dimensions ({} x {}) do not match n = {} and k_use = {}",
            wj.nrows(),
            wj.ncols(),
            n,
            k_use
        )));
    }

    let nthreads = match nthreads {
        Nullable::Null => None,
        Nullable::NotNull(t) => {
            if t <= 0 {
                return Err(Error::Other(
                    "`nthreads` must be a positive integer.".into(),
                ));
            }
            Some(t as usize)
        }
    };

    let coef_slice = coef
        .as_robj()
        .as_real_slice()
        .ok_or_else(|| Error::ExpectedReal(coef.as_robj().clone()))?;

    let coords_rm = col_major_to_row_major(coords.data(), n, d)?;
    let wj_rm = col_major_to_row_major(wj.data(), n, k_use)?;

    let (xj_rm, xjs_rm, proj_u_rm, proj_v_rm, proj_rank) = if uncond {
        (None, None, None, None, 0usize)
    } else {
        let xj_mat = match xj {
            Nullable::NotNull(m) => m,
            Nullable::Null => return Err(Error::Other("conditional SCPC requires `xj` matrix".into())),
        };
        let xjs_mat = match xjs {
            Nullable::NotNull(m) => m,
            Nullable::Null => {
                return Err(Error::Other("conditional SCPC requires `xjs` matrix".into()))
            }
        };
        if xj_mat.nrows() != n || xj_mat.ncols() != k_use {
            return Err(Error::Other(format!(
                "xj dimensions ({} x {}) do not match n = {} and k_use = {}",
                xj_mat.nrows(),
                xj_mat.ncols(),
                n,
                k_use
            )));
        }
        if xjs_mat.nrows() != n || xjs_mat.ncols() != k_use {
            return Err(Error::Other(format!(
                "xjs dimensions ({} x {}) do not match n = {} and k_use = {}",
                xjs_mat.nrows(),
                xjs_mat.ncols(),
                n,
                k_use
            )));
        }
        let xj_rm = col_major_to_row_major(xj_mat.data(), n, k_use)?;
        let xjs_rm = col_major_to_row_major(xjs_mat.data(), n, k_use)?;

        match (proj_u, proj_v) {
            (Nullable::Null, Nullable::Null) => (Some(xj_rm), Some(xjs_rm), None, None, 0usize),
            (Nullable::NotNull(u_mat), Nullable::NotNull(v_mat)) => {
                let rank = u_mat.ncols();
                if u_mat.nrows() != n || v_mat.nrows() != n || v_mat.ncols() != rank {
                    return Err(Error::Other(
                        "conditional projection U/V dimensions do not match".into(),
                    ));
                }
                if rank == 0 {
                    (Some(xj_rm), Some(xjs_rm), None, None, 0usize)
                } else {
                    let u_slice = u_mat.data();
                    let v_slice = v_mat.data();
                    let u_rm = col_major_to_row_major(u_slice, n, rank)?;
                    let v_rm = if std::ptr::eq(u_slice.as_ptr(), v_slice.as_ptr()) {
                        None
                    } else {
                        Some(col_major_to_row_major(v_slice, n, rank)?)
                    };
                    (Some(xj_rm), Some(xjs_rm), Some(u_rm), v_rm, rank)
                }
            }
            _ => {
                return Err(Error::Other(
                    "proj_u and proj_v must both be provided or both be NULL".into(),
                ))
            }
        }
    };

    let metric = if latlong {
        Metric::GreatCircleFractionLonLat
    } else {
        Metric::Euclidean
    };

    let proj_u_ref = proj_u_rm.as_deref();
    let proj_v_ref = proj_v_rm.as_deref().or(proj_u_ref);

    let req = ExactScpcRequest {
        coords: &coords_rm,
        n,
        d,
        metric,
        avc,
        uncond,
        cvs,
        coef: coef_slice,
        wj: &wj_rm,
        xj: xj_rm.as_deref(),
        xjs: xjs_rm.as_deref(),
        proj_u: proj_u_ref,
        proj_v: proj_v_ref,
        proj_rank,
        cv_tolerance,
        nthreads,
    };

    let out = run_exact_scpc(&req).map_err(core_err)?;

    let scpcstats = row_major_to_rmatrix(&out.scpcstats, k_use, 6)?;
    let scpccvs: Robj = match out.scpccvs {
        Some(ref cvs_rm) => row_major_to_rmatrix(cvs_rm, k_use, 4)?.into(),
        None => r!(NULL),
    };
    let wfin = row_major_to_rmatrix(&out.wfin, n, out.q + 1)?;
    let q_i32 = i32::try_from(out.q)
        .map_err(|_| Error::Other("selected q exceeds 32-bit R integer".into()))?;

    List::from_names_and_values(
        [
            "scpcstats",
            "scpccvs",
            "wfin",
            "c0",
            "cmax",
            "cvfin",
            "q",
            "tied",
        ],
        [
            Robj::from(scpcstats),
            scpccvs,
            Robj::from(wfin),
            Robj::from(out.c0),
            Robj::from(out.cmax),
            Robj::from(out.cvfin),
            Robj::from(q_i32),
            Robj::from(out.tied),
        ],
    )
}

#[cfg(feature = "exports")]
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_exact_scpc(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    uncond: bool,
    cvs: bool,
    coef: Doubles,
    wj: RMatrix<f64>,
    xj: Nullable<RMatrix<f64>>,
    xjs: Nullable<RMatrix<f64>>,
    proj_u: Nullable<RMatrix<f64>>,
    proj_v: Nullable<RMatrix<f64>>,
    cv_tolerance: f64,
    nthreads: Nullable<i32>,
) -> extendr_api::Result<List> {
    run_exact_scpc_r(
        coords,
        latlong,
        avc,
        uncond,
        cvs,
        coef,
        wj,
        xj,
        xjs,
        proj_u,
        proj_v,
        cv_tolerance,
        nthreads,
    )
}

#[cfg(feature = "exports")]
extendr_module! {
    mod spatial_r;
    fn rust_run_exact_scpc;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn col_major_to_row_major_transposes_correctly() {
        // R matrix 2 x 3:
        // [1 3 5]
        // [2 4 6]
        let col_major = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let row_major = col_major_to_row_major(&col_major, 2, 3).unwrap();
        assert_eq!(row_major, vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0]);
    }

    #[test]
    fn col_major_to_row_major_rejects_mismatched_length() {
        assert!(col_major_to_row_major(&[1.0, 2.0], 2, 2).is_err());
    }

    #[test]
    fn core_err_maps_all_error_kinds() {
        let cases = [
            (CoreError::invalid_input("bad arg"), "bad arg"),
            (CoreError::unsupported("approx"), "unsupported case: approx"),
            (CoreError::tied("ambiguous"), "tied eigenspace: ambiguous"),
            (CoreError::calibration("no root"), "calibration failed: no root"),
            (CoreError::degenerate("zero var"), "degenerate variance: zero var"),
            (CoreError::internal("lapack"), "internal failure: lapack"),
        ];
        for (err, expected) in cases {
            match core_err(err) {
                Error::Other(msg) => assert_eq!(msg, expected),
                other => panic!("unexpected Error variant: {other:?}"),
            }
        }
    }
}
