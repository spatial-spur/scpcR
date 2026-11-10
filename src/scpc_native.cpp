#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <limits>
#include <unordered_map>
#include <vector>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

#ifndef FCONE
#define FCONE
#endif

#ifdef __APPLE__
extern "C" {
void vvexp(double* y, const double* x, const int* n);
void vvsin(double* y, const double* x, const int* n);
void vvasin(double* y, const double* x, const int* n);
void vvsqrt(double* y, const double* x, const int* n);
}
#endif

namespace {

void symmetrize_lower_to_upper_blocked(double* S, int n) {
  constexpr int BLK = 64;
  for (int j0 = 0; j0 < n; j0 += BLK) {
    const int j_end = std::min(n, j0 + BLK);
    for (int i0 = j0; i0 < n; i0 += BLK) {
      const int i_end = std::min(n, i0 + BLK);
      if (i0 == j0) {
        for (int j = j0; j < j_end; ++j) {
          const double* src_col = S + static_cast<size_t>(j) * n;
          for (int i = j + 1; i < i_end; ++i) {
            S[j + static_cast<size_t>(i) * n] = src_col[i];
          }
        }
      } else {
        for (int j = j0; j < j_end; ++j) {
          const double* src_col = S + static_cast<size_t>(j) * n;
          for (int i = i0; i < i_end; ++i) {
            S[j + static_cast<size_t>(i) * n] = src_col[i];
          }
        }
      }
    }
  }
}

void exp_lower_triangle_inplace(const double* distmat, double* S, int n, double neg_c) {
  for (int i = 0; i < n; ++i) {
    const double* d_col = distmat + static_cast<size_t>(i) * n;
    double* s_col = S + static_cast<size_t>(i) * n;
    const double d_ii = d_col[i];
    s_col[i] = (d_ii == 0.0) ? 1.0 : std::exp(neg_c * d_ii);
    const int len = n - 1 - i;
    if (len > 0) {
      const double* src = d_col + i + 1;
      double* dst = s_col + i + 1;
#ifdef __APPLE__
      for (int k = 0; k < len; ++k) {
        dst[k] = neg_c * src[k];
      }
      vvexp(dst, dst, &len);
#else
      for (int k = 0; k < len; ++k) {
        dst[k] = std::exp(neg_c * src[k]);
      }
#endif
    }
  }
}

// Double-demean an r x c matrix in-place (row means removed first, then column means)
void demean_rect_inplace(double* mat, int r, int c) {
  if (r <= 0 || c <= 0) return;
  const double inv_c = 1.0 / static_cast<double>(c);
  const double inv_r = 1.0 / static_cast<double>(r);
  std::vector<double> row_mean(r, 0.0);

  for (int j = 0; j < c; ++j) {
    const double* col = mat + static_cast<size_t>(j) * r;
    for (int i = 0; i < r; ++i) {
      row_mean[i] += col[i];
    }
  }
  for (int i = 0; i < r; ++i) {
    row_mean[i] *= inv_c;
  }

  for (int j = 0; j < c; ++j) {
    double* col = mat + static_cast<size_t>(j) * r;
    double csum = 0.0;
    for (int i = 0; i < r; ++i) {
      const double v = col[i] - row_mean[i];
      col[i] = v;
      csum += v;
    }
    const double cmean = csum * inv_r;
    for (int i = 0; i < r; ++i) {
      col[i] -= cmean;
    }
  }
}

// Unit-sphere positions let the ordinary great-circle cross-kernel use
// arithmetic on chords instead of repeating sine/cosine for every pair.
void spherical_positions(const double* coords, int n, double* unit) {
  const double radians = M_PI / 180.0;
  for (int i = 0; i < n; ++i) {
    const double longitude = coords[i] * radians;
    const double latitude = coords[i + n] * radians;
    const double cosine = std::cos(latitude);
    unit[i] = cosine * std::cos(longitude);
    unit[i + n] = cosine * std::sin(longitude);
    unit[i + 2 * n] = std::sin(latitude);
  }
}

// Project a bounded row block with the same metric as the subsample covariance.
void compute_cross_kernel(const double* S, int n, int d, const double* sub,
                          int m, bool latlong, double neg_c, const double* unit,
                          const double* sub_unit, double* out, int row_begin, int rows) {
  const double radians = M_PI / 180.0;
  for (int j = 0; j < m; ++j) {
    double* col = out + static_cast<size_t>(j) * rows;
    if (latlong) {
      const double sx = sub_unit[j];
      const double sy = sub_unit[j + m];
      const double sz = sub_unit[j + 2 * m];
      for (int row = 0; row < rows; ++row) {
        const int i = row_begin + row;
        const double x = unit[i], y = unit[i + n], z = unit[i + 2 * n];
        const double dx = x - sx, dy = y - sy, dz = z - sz;
        const double chord2 = dx * dx + dy * dy + dz * dz;
        double distance;
        if (chord2 < 1e-8) {
          // Subtract raw degrees before conversion to retain accuracy for
          // co-located/nearby points, including points across the dateline.
          const double dlat = std::sin(0.5 * (S[i + n] - sub[j + m]) * radians);
          const double dlon = std::sin(0.5 * std::remainder(S[i] - sub[j], 360.0) * radians);
          const double h = dlat * dlat + std::cos(S[i + n] * radians) *
              std::cos(sub[j + m] * radians) * dlon * dlon;
          distance = std::asin(std::sqrt(std::max(0.0, std::min(1.0, h)))) / M_PI;
        } else if (chord2 > 4.0 - 1e-8) {
          // asin(chord/2) loses angular precision near antipodes.
          const double cx = y * sz - z * sy;
          const double cy = z * sx - x * sz;
          const double cz = x * sy - y * sx;
          distance = std::atan2(std::sqrt(cx * cx + cy * cy + cz * cz),
                                x * sx + y * sy + z * sz) / (2.0 * M_PI);
        } else {
          distance = std::asin(0.5 * std::sqrt(chord2)) / M_PI;
        }
        col[row] = neg_c * distance;
      }
    } else if (d == 2) {
      const double* x = S + row_begin;
      const double* y = S + n + row_begin;
      const double sx = sub[j];
      const double sy = sub[j + m];
      for (int row = 0; row < rows; ++row) {
        const double dx = x[row] - sx;
        const double dy = y[row] - sy;
        col[row] = neg_c * std::sqrt(dx * dx + dy * dy);
      }
    } else {
      for (int row = 0; row < rows; ++row) {
        const int i = row_begin + row;
        double squared = 0.0;
        for (int k = 0; k < d; ++k) {
          const double diff = S[i + static_cast<size_t>(k) * n] - sub[j + static_cast<size_t>(k) * m];
          squared += diff * diff;
        }
        col[row] = neg_c * std::sqrt(squared);
      }
    }
#ifdef __APPLE__
    vvexp(col, col, &rows);
#else
    for (int row = 0; row < rows; ++row) col[row] = std::exp(col[row]);
#endif
  }
}

void compute_distmat_raw(const double* s, int n, int d, bool latlong, double* out) {
  if (latlong) {
    // In scpcR, column 0 is lon and column 1 is lat (in degrees)
    const double* lon = s;
    const double* lat = s + n;
    const double deg2rad = M_PI / 180.0;
    const double inv_pi = 1.0 / M_PI;

    std::vector<double> lon_rad(n);
    std::vector<double> lat_rad(n);
    std::vector<double> cos_lat(n);
    for (int i = 0; i < n; ++i) {
      lon_rad[i] = lon[i] * deg2rad;
      lat_rad[i] = lat[i] * deg2rad;
      cos_lat[i] = std::cos(lat_rad[i]);
    }

    for (int i = 0; i < n; ++i) {
      double* out_col = out + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      const double lon_i = lon_rad[i];
      const double lat_i = lat_rad[i];
      const double cos_i = cos_lat[i];
      for (int j = i + 1; j < n; ++j) {
        const double dlat = 0.5 * (lat_rad[j] - lat_i);
        const double dlon = 0.5 * (lon_rad[j] - lon_i);
        const double sdlat = std::sin(dlat);
        const double sdlon = std::sin(dlon);
        double inside = sdlat * sdlat + cos_i * cos_lat[j] * (sdlon * sdlon);
        if (inside > 1.0) inside = 1.0;
        if (inside < 0.0) inside = 0.0;
        out_col[j] = std::asin(std::sqrt(inside)) * inv_pi;
      }
    }
  } else if (d == 2) {
    const double* s0 = s;
    const double* s1 = s + n;
    for (int i = 0; i < n; ++i) {
      double* out_col = out + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      const double x0 = s0[i];
      const double y0 = s1[i];
      for (int j = i + 1; j < n; ++j) {
        const double dx = s0[j] - x0;
        const double dy = s1[j] - y0;
        out_col[j] = std::sqrt(dx * dx + dy * dy);
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double* out_col = out + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      for (int j = i + 1; j < n; ++j) {
        double sum_sq = 0.0;
        for (int k = 0; k < d; ++k) {
          const double diff = s[j + static_cast<size_t>(k) * n] - s[i + static_cast<size_t>(k) * n];
          sum_sq += diff * diff;
        }
        out_col[j] = std::sqrt(sum_sq);
      }
    }
  }
  symmetrize_lower_to_upper_blocked(out, n);
}

struct RootBracket {
  double lower;
  double upper;
  double value_upper;
};

// Illinois interpolation with a bisection safeguard. Retain the enclosing
// bracket and its conservative upper endpoint; interpolation changes the work,
// not the requested absolute/relative accuracy.
template <typename Evaluate>
RootBracket refine_monotone_root(double lower, double upper,
                                  double value_lower, double value_upper,
                                  double target, double absolute_tolerance,
                                  double relative_tolerance, Evaluate evaluate) {
  if (value_lower == target) return {lower, lower, value_lower};
  if (value_upper == target) return {upper, upper, value_upper};
  double weight_lower = value_lower - target;
  double weight_upper = value_upper - target;
  double checkpoint_width = upper - lower;
  int interpolation_steps = 0;
  int last_side = 0;
  while (upper - lower > absolute_tolerance + relative_tolerance * upper) {
    bool bisect = false;
    if (interpolation_steps == 3) {
      bisect = upper - lower > 0.5 * checkpoint_width;
      checkpoint_width = upper - lower;
      interpolation_steps = 0;
    }
    double candidate = bisect ? lower + 0.5 * (upper - lower) :
        lower + (upper - lower) * (weight_lower / (weight_lower - weight_upper));
    if (!(candidate > lower && candidate < upper)) candidate = lower + 0.5 * (upper - lower);
    // Adjacent floating-point endpoints cannot be refined further.
    if (!(candidate > lower && candidate < upper)) break;
    const double value = evaluate(candidate);
    if (!std::isfinite(value)) throw std::runtime_error("Non-finite value during monotone solve.");
    if (value == target) return {candidate, candidate, value};
    if (value > target) {
      lower = candidate;
      value_lower = value;
      weight_lower = value - target;
      weight_upper = last_side == 1 ? 0.5 * weight_upper : value_upper - target;
      last_side = 1;
    } else {
      upper = candidate;
      value_upper = value;
      weight_upper = value - target;
      weight_lower = last_side == -1 ? 0.5 * weight_lower : value_lower - target;
      last_side = -1;
    }
    ++interpolation_steps;
  }
  return {lower, upper, value_upper};
}

struct HistBin {
  double count;
  double mean;
  double lo;
  double hi;
  double w_lo;
  double w_hi;
};

struct EvalPoint {
  double c;
  double v;
};

void build_hist_bins(const double* vd, size_t m, double min_d, double max_d, std::vector<HistBin>& bins) {
  bins.clear();
  if (m < 4096) return;
  constexpr int NBINS = 1024;
  size_t b_count[NBINS] = {0};
  double b_sum[NBINS] = {0.0};
  double b_min[NBINS];
  double b_max[NBINS];
  for (int b = 0; b < NBINS; ++b) {
    b_min[b] = std::numeric_limits<double>::infinity();
    b_max[b] = -std::numeric_limits<double>::infinity();
  }
  const double span = max_d - min_d;
  const double inv_span = (span > 0.0) ? (static_cast<double>(NBINS) / span) : 0.0;
  for (size_t k = 0; k < m; ++k) {
    const double d = vd[k];
    int b = static_cast<int>((d - min_d) * inv_span);
    if (b < 0) b = 0;
    else if (b >= NBINS) b = NBINS - 1;
    b_count[b]++;
    b_sum[b] += d;
    if (d < b_min[b]) b_min[b] = d;
    if (d > b_max[b]) b_max[b] = d;
  }
  bins.reserve(NBINS);
  for (int b = 0; b < NBINS; ++b) {
    if (b_count[b] > 0) {
      const double cnt = static_cast<double>(b_count[b]);
      const double mu = b_sum[b] / cnt;
      const double lo = b_min[b];
      const double hi = b_max[b];
      double w_lo = 1.0;
      double w_hi = 0.0;
      if (hi > lo) {
        w_hi = (mu - lo) / (hi - lo);
        w_lo = (hi - mu) / (hi - lo);
      }
      bins.push_back({cnt, mu, lo, hi, w_lo, w_hi});
    }
  }
}

void solve_c0_from_avc_core(const double* vd,
                            size_t m,
                            const std::vector<HistBin>& bins,
                            std::vector<EvalPoint>& history,
                            const double* avc_targets,
                            int n_targets,
                            double* out) {
  if (m == 0) throw std::runtime_error("At least two distinct locations are required.");
  double distance_sum = 0.0;
  size_t zero_count = 0;
  for (size_t i = 0; i < m; ++i) {
    if (!std::isfinite(vd[i]) || vd[i] < 0.0) throw std::runtime_error("Distances must be finite and non-negative.");
    distance_sum += vd[i];
    if (vd[i] == 0.0) ++zero_count;
  }
  const double mean_distance = distance_sum / static_cast<double>(m);
  if (!(mean_distance > 0.0) || !std::isfinite(mean_distance)) throw std::runtime_error("At least two distinct locations are required.");
  const double correlation_floor = static_cast<double>(zero_count) / static_cast<double>(m);
  const double inv_m = 1.0 / static_cast<double>(m);
  constexpr int CHUNK = 2048;
  double arg_buf[CHUNK];
  double exp_buf[CHUNK];

  auto eval_convex_bounds = [&](double c, double& lb_out, double& ub_out) {
    const int nb = static_cast<int>(bins.size());
    const double neg_c = -c;
    double lb_sum = 0.0;
    double ub_sum = 0.0;
    for (int idx = 0; idx < nb; ++idx) {
      const auto& bin = bins[idx];
      const double e_mu = std::exp(neg_c * bin.mean);
      lb_sum += bin.count * e_mu;
      if (bin.hi > bin.lo) {
        const double e_lo = std::exp(neg_c * bin.lo);
        const double e_hi = std::exp(neg_c * bin.hi);
        ub_sum += bin.count * (bin.w_lo * e_lo + bin.w_hi * e_hi);
      } else {
        ub_sum += bin.count * e_mu;
      }
    }
    lb_out = lb_sum * inv_m;
    ub_out = ub_sum * inv_m;
  };

  auto eval_mean_exp = [&](double c) -> double {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end() && it->c == c) {
      return it->v;
    }
    double sum = 0.0;
    const double neg_c = -c;
#ifdef __APPLE__
    size_t k = 0;
    while (k < m) {
      const int b = static_cast<int>(std::min<size_t>(CHUNK, m - k));
      for (int i = 0; i < b; ++i) {
        arg_buf[i] = neg_c * vd[k + i];
      }
      vvexp(exp_buf, arg_buf, &b);
      for (int i = 0; i < b; ++i) {
        sum += exp_buf[i];
      }
      k += b;
    }
#else
    for (size_t k = 0; k < m; ++k) {
      sum += std::exp(neg_c * vd[k]);
    }
#endif
    const double val = sum * inv_m;
    history.insert(it, {c, val});
    return val;
  };

  // Returns true iff mean(exp(-c * vd)) < target
  auto is_below = [&](double c, double target) -> bool {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end()) {
      if (it->c == c) return it->v < target;
      if (it->v >= target) return false;
    }
    if (it != history.begin()) {
      auto prev = it - 1;
      if (prev->v < target) return true;
    }
    if (!bins.empty()) {
      double lb = 0.0, ub = 1.0;
      eval_convex_bounds(c, lb, ub);
      if (ub + 1e-12 < target) return true;
      if (lb - 1e-12 > target) return false;
    }
    return eval_mean_exp(c) < target;
  };

  // Returns true iff mean(exp(-c * vd)) > target
  auto is_above = [&](double c, double target) -> bool {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end()) {
      if (it->c == c) return it->v > target;
      if (it->v > target) return true;
    }
    if (it != history.begin()) {
      auto prev = it - 1;
      if (prev->v <= target) return false;
    }
    if (!bins.empty()) {
      double lb = 0.0, ub = 1.0;
      eval_convex_bounds(c, lb, ub);
      if (lb - 1e-12 > target) return true;
      if (ub + 1e-12 < target) return false;
    }
    return eval_mean_exp(c) > target;
  };

  for (int t = 0; t < n_targets; ++t) {
    const double avc0 = avc_targets[t];
    if (!std::isfinite(avc0) || avc0 <= correlation_floor || avc0 >= 1.0) {
      throw std::runtime_error("Target correlation must lie above the co-located-pair fraction and below one.");
    }
    double c0 = 1.0 / mean_distance;
    double c1 = c0;
    while (is_below(c0, avc0)) {
      c1 = c0;
      c0 = 0.5 * c0;
    }
    while (is_above(c1, avc0)) {
      c0 = c1;
      c1 = 2.0 * c1;
      if (!std::isfinite(c1)) throw std::runtime_error("Spatial scale is not representable.");
    }
    const RootBracket root = refine_monotone_root(
        c0, c1, eval_mean_exp(c0), eval_mean_exp(c1), avc0, 0.0, 1e-12, eval_mean_exp);
    out[t] = root.upper;
  }
}

// Cache the first three angular panels across all root evaluations in a
// workspace. The rule itself is checked so callers can safely replace it.
struct AngularRuleCache {
  std::vector<double> nodes;
  std::vector<double> weights;
  std::vector<double> cosine_squared;
  std::vector<double> sine_squared;
  void prepare(const double* x, const double* w, int n) {
    if (static_cast<int>(nodes.size()) == n &&
        std::equal(nodes.begin(), nodes.end(), x) &&
        std::equal(weights.begin(), weights.end(), w)) return;
    nodes.assign(x, x + n);
    weights.assign(w, w + n);
    cosine_squared.resize(static_cast<size_t>(3) * n);
    sine_squared.resize(static_cast<size_t>(3) * n);
    const double end = 0.5 * M_PI;
    for (int panel = 0; panel < 3; ++panel) {
      const double left = panel == 2 ? 0.5 * end : 0.0;
      const double width = panel == 0 ? end : 0.5 * end;
      for (int j = 0; j < n; ++j) {
        const double angle = left + width * x[j];
        const double cosine = std::cos(angle), sine = std::sin(angle);
        cosine_squared[static_cast<size_t>(panel) * n + j] = cosine * cosine;
        sine_squared[static_cast<size_t>(panel) * n + j] = sine * sine;
      }
    }
  }
};

// Every factor is in [0, 1], so the product cannot overflow. If the product
// underflows, its square root is below about 1e-161, far below the integration
// error budget. This avoids a logarithm for every eigenvalue at every node.
inline double angular_integrand(double cosine2, const double* eigenvalues, int n) {
  double product = 1.0;
  for (int i = 0; i < n; ++i) product *= cosine2 / (cosine2 + eigenvalues[i]);
  return std::sqrt(product);
}

// Stable complement of the square root of a product of factors in [0, 1].
inline double angular_deficit(double sine2, double denominator,
                               const double* normalized_eigenvalues, int n) {
  double product = 1.0, deficit = 0.0;
  for (int i = 0; i < n; ++i) {
    const double extra = normalized_eigenvalues[i] * denominator;
    const double inverse = 1.0 / (sine2 + extra);
    deficit += product * extra * inverse;
    product *= sine2 * inverse;
  }
  return deficit / ((1.0 + std::sqrt(product)) * denominator);
}

struct RpWorkspace {
  std::vector<double> omx;
  std::vector<double> wr;
  std::vector<double> wi;
  std::vector<double> work;
  std::vector<double> evals;
  AngularRuleCache angular;

  void ensure(int k) {
    if (static_cast<int>(wr.size()) < k) {
      omx.resize(static_cast<size_t>(k) * k);
      wr.resize(k);
      wi.resize(k);
      evals.resize(k);
      work.resize(std::max(1, 8 * k));
    }
  }
};

double eval_getrp_submatrix(const double* om_full,
                            int ld_om,
                            int k,
                            double cv,
                            const double* gqx,
                            const double* gqw,
                            int n_gq,
                            RpWorkspace& ws) {
  ws.ensure(k);
  const double neg_cv2 = -(cv * cv);
  double* omx = ws.omx.data();

  for (int j = 0; j < k; ++j) {
    const double* src_col = om_full + static_cast<size_t>(j) * ld_om;
    double* dst_col = omx + static_cast<size_t>(j) * k;
    dst_col[0] = src_col[0];
    for (int i = 1; i < k; ++i) {
      dst_col[i] = neg_cv2 * src_col[i];
    }
  }

  const char jobvl = 'N';
  const char jobvr = 'N';
  int n_lapack = k;
  int lwork = static_cast<int>(ws.work.size());
  int info = 0;
  double dummy = 0.0;
  int ldummy = 1;

  F77_CALL(dgeev)(&jobvl, &jobvr, &n_lapack, omx, &n_lapack,
                  ws.wr.data(), ws.wi.data(),
                  &dummy, &ldummy, &dummy, &ldummy,
                  ws.work.data(), &lwork, &info FCONE FCONE);

  int n_neg = 0;
  double denom = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < k; ++i) {
    const double v = ws.wr[i];
    if (v > denom) denom = v;
    if (v < 0.0) {
      ws.evals[n_neg++] = -v;
    }
  }

  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  if (n_neg == 0) return 1.0;

  for (int i = 0; i < n_neg; ++i) ws.evals[i] /= denom;
  if (n_neg == 1) {
    // Ratio of two independent chi-square(1) variables: exact Cauchy tail.
    return (2.0 / M_PI) * std::atan(1.0 / std::sqrt(ws.evals[0]));
  }

  // The angular integral is smooth at its endpoints for ordinary cutoffs.
  // Near probability one, rescale its complement using cot(theta) = k*tan(phi),
  // k^2 = sum(eigenvalues). This resolves the otherwise invisible endpoint
  // boundary layer without subtracting nearly equal factors in the integrand.
  double scale2 = 0.0;
  for (int i = 0; i < n_neg; ++i) scale2 += ws.evals[i];
  if (scale2 == 0.0) return 1.0;
  const bool complement = scale2 < 0.0625;
  if (complement) for (int i = 0; i < n_neg; ++i) ws.evals[i] /= scale2;
  ws.angular.prepare(gqx, gqw, n_gq);
  auto integrate = [&](double left, double right, int panel = -1) {
    double total = 0.0;
    for (int j = 0; j < n_gq; ++j) {
      double cosine2, sine2 = 0.0;
      if (panel >= 0) {
        const size_t index = static_cast<size_t>(panel) * n_gq + j;
        cosine2 = ws.angular.cosine_squared[index];
        if (complement) sine2 = ws.angular.sine_squared[index];
      } else {
        const double angle = left + (right - left) * gqx[j];
        const double cosine = std::cos(angle);
        cosine2 = cosine * cosine;
        if (complement) { const double sine = std::sin(angle); sine2 = sine * sine; }
      }
      const double value = complement ?
          angular_deficit(sine2, cosine2 + scale2 * sine2, ws.evals.data(), n_neg) :
          angular_integrand(cosine2, ws.evals.data(), n_neg);
      total += gqw[j] * value;
    }
    return (right - left) * total;
  };
  struct Interval { double left, right, whole, tolerance; int depth; };
  const double end = 0.5 * M_PI;
  std::vector<Interval> pending{{0.0, end, integrate(0.0, end, 0), 1e-11, 0}};
  double total = 0.0;
  while (!pending.empty()) {
    const Interval current = pending.back();
    pending.pop_back();
    const double mid = 0.5 * (current.left + current.right);
    const double left = integrate(current.left, mid, current.depth == 0 ? 1 : -1);
    const double right = integrate(mid, current.right, current.depth == 0 ? 2 : -1);
    if (std::abs(left + right - current.whole) <= current.tolerance || current.depth >= 20) {
      total += left + right;
    } else {
      pending.push_back({mid, current.right, right, current.tolerance * 0.5, current.depth + 1});
      pending.push_back({current.left, mid, left, current.tolerance * 0.5, current.depth + 1});
    }
  }
  const double probability = complement ?
      1.0 - std::sqrt(scale2) * total * (2.0 / M_PI) : total * (2.0 / M_PI);
  return std::max(0.0, std::min(1.0, probability));
}

void eval_maxrp_core(const std::vector<const double*>& oms_ptrs,
                     const std::vector<int>& oms_lds,
                     int q,
                     double cv,
                     const double* gqx,
                     const double* gqw,
                     int n_gq,
                     RpWorkspace& ws,
                     int known_idx,
                     double known_rp,
                     double& max_rp_out,
                     int& max_idx_out) {
  const int nc = static_cast<int>(oms_ptrs.size());
  const int k = q + 1;
  double best_rp = -1.0;
  int best_i = 0;

  for (int m = 0; m < nc; ++m) {
    const double rp = (m == known_idx)
                          ? known_rp
                          : eval_getrp_submatrix(oms_ptrs[m], oms_lds[m], k, cv, gqx, gqw, n_gq, ws);
    if (m == 0 || rp > best_rp) {
      best_rp = rp;
      best_i = m;
    }
  }
  max_rp_out = best_rp;
  max_idx_out = best_i;
}

double eval_getcv_core(const std::vector<const double*>& oms_ptrs,
                       const std::vector<int>& oms_lds,
                       int q,
                       double level,
                       const double* gqx,
                       const double* gqw,
                       int n_gq,
                       RpWorkspace& ws) {
  const int k = q + 1;
  const double sqrt_q = std::sqrt(static_cast<double>(q));
  const double inv_sqrt_q = 1.0 / sqrt_q;
  const double tol = 1e-9 * inv_sqrt_q;

  double rp = 1.0;
  int i = 0;
  double cv0 = Rf_qt(1.0 - level * 0.5, static_cast<double>(q), 1, 0) * inv_sqrt_q;
  double cv1 = cv0;
  bool first_iter = true;

  while (rp > level) {
    cv1 = cv0;
    double rp_cv1 = rp;
    double rp_cv0 = rp;
    if (!first_iter) {
      // At cv1 == cv0, we already know getrp(Oms[i], cv1) == rp > level
      cv0 = cv1;
      cv1 = cv1 + inv_sqrt_q;
    }
    first_iter = false;

    while (true) {
      rp_cv1 = eval_getrp_submatrix(oms_ptrs[i], oms_lds[i], k, cv1, gqx, gqw, n_gq, ws);
      if (rp_cv1 > level) {
        rp_cv0 = rp_cv1;
        cv0 = cv1;
        cv1 = cv1 + inv_sqrt_q;
      } else {
        break;
      }
    }

    const RootBracket root = refine_monotone_root(
        cv0, cv1, rp_cv0, rp_cv1, level, tol, 0.0,
        [&](double cv) { return eval_getrp_submatrix(oms_ptrs[i], oms_lds[i], k, cv, gqx, gqw, n_gq, ws); });
    cv1 = root.upper;
    rp_cv1 = root.value_upper;

    double max_rp = 0.0;
    int max_i = 0;
    eval_maxrp_core(oms_ptrs, oms_lds, q, cv1, gqx, gqw, n_gq, ws, i, rp_cv1, max_rp, max_i);
    if (max_i == i) break;
    i = max_i;
    rp = max_rp;
    cv0 = cv1;
  }

  return cv1 * sqrt_q;
}

void compute_top_k_evecs_dsyevr(double* A, int n, int qmax, double* V_out) {
  const char jobz = 'V';
  const char range = 'I';
  const char uplo = 'L';
  const double vl = 0.0;
  const double vu = 0.0;
  const int il = std::max(1, n - qmax + 1);
  const int iu = n;
  const int k = iu - il + 1;
  const double abstol = 0.0;
  int m_found = 0;
  int info = 0;

  std::vector<double> evals(n);
  std::vector<double> Z(static_cast<size_t>(n) * k);
  std::vector<int> isuppz(2 * std::max(1, k));

  double work_query = 0.0;
  int iwork_query = 0;
  int lwork = -1;
  int liwork = -1;

  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A, &n, &vl, &vu, &il, &iu,
                   &abstol, &m_found, evals.data(), Z.data(), &n, isuppz.data(),
                   &work_query, &lwork, &iwork_query, &liwork, &info FCONE FCONE FCONE);

  lwork = std::max(26 * n, static_cast<int>(work_query) + 1);
  liwork = std::max(10 * n, iwork_query + 1);
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);

  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A, &n, &vl, &vu, &il, &iu,
                   &abstol, &m_found, evals.data(), Z.data(), &n, isuppz.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE FCONE);

  for (int c = 0; c < k; ++c) {
    const double* src = Z.data() + static_cast<size_t>(k - 1 - c) * n;
    double* dst = V_out + static_cast<size_t>(c) * n;
    std::memcpy(dst, src, static_cast<size_t>(n) * sizeof(double));
  }
}

void compute_top_k_evecs_sym(double* A, int n, int qmax, double* V_out) {
  const int k = std::min(qmax, n);
  const int m_k = std::min(n, std::max(2 * k + 16, 44));
  // Use direct LAPACK dsyevr for small/medium matrices or when k is a large fraction of n
  if (n <= 450 || m_k >= n / 2) {
    compute_top_k_evecs_dsyevr(A, n, k, V_out);
    return;
  }

  // Thick-Restart Lanczos with full DGKS reorthogonalization using BLAS dgemv/dgemm
  std::vector<double> V(static_cast<size_t>(n) * (m_k + 1), 0.0);
  std::vector<double> H(static_cast<size_t>(m_k) * m_k, 0.0);
  std::vector<double> H_copy(static_cast<size_t>(m_k) * m_k, 0.0);
  std::vector<double> w(n, 0.0);
  std::vector<double> h_vec(m_k + 1, 0.0);
  std::vector<double> c_vec(m_k + 1, 0.0);
  std::vector<double> ritz_vals(m_k, 0.0);
  std::vector<double> ritz_vecs(static_cast<size_t>(m_k) * m_k, 0.0);
  std::vector<double> U_keep(static_cast<size_t>(m_k) * m_k, 0.0);
  std::vector<double> V_new(static_cast<size_t>(n) * m_k, 0.0);

  double norm0 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double v = std::sin(0.1 + 1.7 * static_cast<double>(i + 1));
    V[i] = v;
    norm0 += v * v;
  }
  norm0 = 1.0 / std::sqrt(norm0);
  for (int i = 0; i < n; ++i) {
    V[i] *= norm0;
  }

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  const double neg_one = -1.0;
  const int inc1 = 1;

  int k_start = 0;
  double beta_last = 0.0;
  bool converged = false;
  constexpr int MAX_RESTARTS = 25;

  for (int restart = 0; restart < MAX_RESTARTS; ++restart) {
    for (int j = k_start; j < m_k; ++j) {
      const double* v_j = V.data() + static_cast<size_t>(j) * n;
      F77_CALL(dgemv)(&trans_n, &n, &n, &one, A, &n, v_j, &inc1, &zero, w.data(), &inc1 FCONE);

      const int cols_j = j + 1;
      F77_CALL(dgemv)(&trans_t, &n, &cols_j, &one, V.data(), &n, w.data(), &inc1, &zero, h_vec.data(), &inc1 FCONE);
      F77_CALL(dgemv)(&trans_n, &n, &cols_j, &neg_one, V.data(), &n, h_vec.data(), &inc1, &one, w.data(), &inc1 FCONE);

      // Second DGKS pass
      F77_CALL(dgemv)(&trans_t, &n, &cols_j, &one, V.data(), &n, w.data(), &inc1, &zero, c_vec.data(), &inc1 FCONE);
      F77_CALL(dgemv)(&trans_n, &n, &cols_j, &neg_one, V.data(), &n, c_vec.data(), &inc1, &one, w.data(), &inc1 FCONE);

      for (int r = 0; r <= j; ++r) {
        const double val = h_vec[r] + c_vec[r];
        H[r + static_cast<size_t>(j) * m_k] = val;
        H[j + static_cast<size_t>(r) * m_k] = val;
      }

      double norm_w = 0.0;
      for (int i = 0; i < n; ++i) {
        norm_w += w[i] * w[i];
      }
      beta_last = std::sqrt(norm_w);

      double* v_next = V.data() + static_cast<size_t>(j + 1) * n;
      if (beta_last < 1e-13) {
        beta_last = 0.0;
        for (int i = 0; i < n; ++i) {
          v_next[i] = std::cos(0.3 + 2.3 * static_cast<double>((j + 1) * (i + 1)));
        }
        F77_CALL(dgemv)(&trans_t, &n, &cols_j, &one, V.data(), &n, v_next, &inc1, &zero, h_vec.data(), &inc1 FCONE);
        F77_CALL(dgemv)(&trans_n, &n, &cols_j, &neg_one, V.data(), &n, h_vec.data(), &inc1, &one, v_next, &inc1 FCONE);
        double nrm = 0.0;
        for (int i = 0; i < n; ++i) nrm += v_next[i] * v_next[i];
        nrm = 1.0 / std::sqrt(nrm);
        for (int i = 0; i < n; ++i) v_next[i] *= nrm;
      } else {
        const double inv_beta = 1.0 / beta_last;
        for (int i = 0; i < n; ++i) {
          v_next[i] = w[i] * inv_beta;
        }
      }
    }

    // Solve m_k x m_k symmetric eigenproblem for H
    std::memcpy(H_copy.data(), H.data(), static_cast<size_t>(m_k) * m_k * sizeof(double));
    const char jobz = 'V';
    const char uplo = 'L';
    int lwork = std::max(1, 3 * m_k);
    std::vector<double> work(lwork);
    int info = 0;
    F77_CALL(dsyev)(&jobz, &uplo, &m_k, H_copy.data(), &m_k, ritz_vals.data(), work.data(), &lwork, &info FCONE FCONE);
    if (info != 0) break;

    const double max_eval = std::max(1e-12, std::abs(ritz_vals[m_k - 1]));
    const double tol = 1e-10 * max_eval;
    bool all_ok = true;
    for (int c = 0; c < k; ++c) {
      const int ritz_idx = m_k - 1 - c;
      const double last_comp = H_copy[(m_k - 1) + static_cast<size_t>(ritz_idx) * m_k];
      if (beta_last * std::abs(last_comp) > tol) {
        all_ok = false;
        break;
      }
    }

    if (all_ok) {
      for (int c = 0; c < k; ++c) {
        const double* src = H_copy.data() + static_cast<size_t>(m_k - 1 - c) * m_k;
        std::memcpy(U_keep.data() + static_cast<size_t>(c) * m_k, src, static_cast<size_t>(m_k) * sizeof(double));
      }
      F77_CALL(dgemm)(&trans_n, &trans_n, &n, &k, &m_k, &one,
                      V.data(), &n, U_keep.data(), &m_k, &zero, V_out, &n FCONE FCONE);
      converged = true;
      break;
    }

    const int k_keep = std::min(m_k - 2, k + std::max(4, (m_k - k) / 2));
    for (int c = 0; c < k_keep; ++c) {
      const double* src = H_copy.data() + static_cast<size_t>(m_k - 1 - c) * m_k;
      std::memcpy(U_keep.data() + static_cast<size_t>(c) * m_k, src, static_cast<size_t>(m_k) * sizeof(double));
    }
    F77_CALL(dgemm)(&trans_n, &trans_n, &n, &k_keep, &m_k, &one,
                    V.data(), &n, U_keep.data(), &m_k, &zero, V_new.data(), &n FCONE FCONE);
    std::memcpy(V.data(), V_new.data(), static_cast<size_t>(n) * k_keep * sizeof(double));
    std::memcpy(V.data() + static_cast<size_t>(k_keep) * n,
                V.data() + static_cast<size_t>(m_k) * n,
                static_cast<size_t>(n) * sizeof(double));

    std::fill(H.begin(), H.end(), 0.0);
    for (int c = 0; c < k_keep; ++c) {
      H[c + static_cast<size_t>(c) * m_k] = ritz_vals[m_k - 1 - c];
    }
    k_start = k_keep;
  }

  if (!converged) {
    compute_top_k_evecs_dsyevr(A, n, k, V_out);
  }
}

}  // namespace

extern "C" {

SEXP cpp_lvech(SEXP mat_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(mat_sexp);
  const double* mat = REAL(mat_sexp);
  const R_xlen_t m = static_cast<R_xlen_t>(n) * (n - 1) / 2;
  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, m));
  double* out = REAL(out_sexp);

  R_xlen_t pos = 0;
  for (int j = 0; j < n - 1; ++j) {
    const double* col = mat + static_cast<size_t>(j) * n + (j + 1);
    const int len = n - 1 - j;
    std::memcpy(out + pos, col, static_cast<size_t>(len) * sizeof(double));
    pos += len;
  }

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_demeanmat(SEXP mat_sexp) {
  const char* error_message = nullptr;
  try {
  const int r = Rf_nrows(mat_sexp);
  const int c = Rf_ncols(mat_sexp);
  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, r, c));
  std::memcpy(REAL(out_sexp), REAL(mat_sexp), static_cast<size_t>(r) * c * sizeof(double));
  demean_rect_inplace(REAL(out_sexp), r, c);
  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_sigma_demean(SEXP distmat_sexp, SEXP c0_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double c0 = Rf_asReal(c0_sexp);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);

  exp_lower_triangle_inplace(distmat, out, n, -c0);
  symmetrize_lower_to_upper_blocked(out, n);
  demean_rect_inplace(out, n, n);

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_getW(SEXP distmat_sexp, SEXP c0_sexp, SEXP qmax_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double c0 = Rf_asReal(c0_sexp);
  const int qmax = std::min(Rf_asInteger(qmax_sexp), std::max(1, n - 1));

  std::vector<double> sig_d(static_cast<size_t>(n) * n);
  exp_lower_triangle_inplace(distmat, sig_d.data(), n, -c0);
  symmetrize_lower_to_upper_blocked(sig_d.data(), n);
  demean_rect_inplace(sig_d.data(), n, n);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, qmax + 1));
  double* out = REAL(out_sexp);
  const double inv_sqrt_n = 1.0 / std::sqrt(static_cast<double>(n));
  for (int i = 0; i < n; ++i) {
    out[i] = inv_sqrt_n;
  }

  compute_top_k_evecs_sym(sig_d.data(), n, qmax, out + n);

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_ln_subset_evecs(SEXP distmat_sexp, SEXP c0_sexp, SEXP qmax_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double c0 = Rf_asReal(c0_sexp);
  const int qmax = std::min(Rf_asInteger(qmax_sexp), n);

  std::vector<double> sig_d(static_cast<size_t>(n) * n);
  exp_lower_triangle_inplace(distmat, sig_d.data(), n, -c0);
  symmetrize_lower_to_upper_blocked(sig_d.data(), n);
  demean_rect_inplace(sig_d.data(), n, n);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, qmax));
  compute_top_k_evecs_sym(sig_d.data(), n, qmax, REAL(out_sexp));

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_distmat(SEXP s_sexp, SEXP latlong_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const double* s = REAL(s_sexp);
  const bool latlong = Rf_asLogical(latlong_sexp) == TRUE;

  if (latlong && d != 2) {
    throw std::runtime_error("Internal error: geodesic coordinates must have exactly two columns.");
  }

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  compute_distmat_raw(s, n, d, latlong, REAL(out_sexp));
  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_distvec(SEXP s1_sexp, SEXP s2_sexp, SEXP latlong_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s1_sexp);
  const int d = Rf_ncols(s1_sexp);
  const double* s1 = REAL(s1_sexp);
  const double* s2 = REAL(s2_sexp);
  const bool latlong = Rf_asLogical(latlong_sexp) == TRUE;

  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, n));
  double* out = REAL(out_sexp);

  if (latlong) {
    const double* lon1 = s1;
    const double* lat1 = s1 + n;
    const double* lon2 = s2;
    const double* lat2 = s2 + n;
    const double deg2rad = M_PI / 180.0;
    const double inv_pi = 1.0 / M_PI;
    for (int i = 0; i < n; ++i) {
      const double la1 = lat1[i] * deg2rad;
      const double la2 = lat2[i] * deg2rad;
      const double dlon = 0.5 * (lon1[i] - lon2[i]) * deg2rad;
      const double dlat = 0.5 * (la1 - la2);
      const double sdlat = std::sin(dlat);
      const double sdlon = std::sin(dlon);
      double inside = sdlat * sdlat + std::cos(la1) * std::cos(la2) * (sdlon * sdlon);
      if (inside > 1.0) inside = 1.0;
      if (inside < 0.0) inside = 0.0;
      out[i] = std::asin(std::sqrt(inside)) * inv_pi;
    }
  } else if (d == 2) {
    const double* x1 = s1;
    const double* y1 = s1 + n;
    const double* x2 = s2;
    const double* y2 = s2 + n;
    for (int i = 0; i < n; ++i) {
      const double dx = x1[i] - x2[i];
      const double dy = y1[i] - y2[i];
      out[i] = std::sqrt(dx * dx + dy * dy);
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double sum_sq = 0.0;
      for (int k = 0; k < d; ++k) {
        const double diff = s1[i + static_cast<size_t>(k) * n] - s2[i + static_cast<size_t>(k) * n];
        sum_sq += diff * diff;
      }
      out[i] = std::sqrt(sum_sq);
    }
  }

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_getavc(SEXP c_sexp, SEXP dist_sexp) {
  const char* error_message = nullptr;
  try {
  const double neg_c = -Rf_asReal(c_sexp);
  const R_xlen_t m = Rf_xlength(dist_sexp);
  const double* vd = REAL(dist_sexp);
  if (m == 0) {
    return Rf_ScalarReal(NA_REAL);
  }
  double sum = 0.0;
#ifdef __APPLE__
  constexpr int CHUNK = 2048;
  double arg_buf[CHUNK];
  double exp_buf[CHUNK];
  R_xlen_t k = 0;
  while (k < m) {
    const int b = static_cast<int>(std::min<R_xlen_t>(CHUNK, m - k));
    for (int i = 0; i < b; ++i) {
      arg_buf[i] = neg_c * vd[k + i];
    }
    vvexp(exp_buf, arg_buf, &b);
    for (int i = 0; i < b; ++i) {
      sum += exp_buf[i];
    }
    k += b;
  }
#else
  for (R_xlen_t k = 0; k < m; ++k) {
    sum += std::exp(neg_c * vd[k]);
  }
#endif
  return Rf_ScalarReal(sum / static_cast<double>(m));
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_getc0_from_avc(SEXP dist_sexp, SEXP avc_vec_sexp) {
  const char* error_message = nullptr;
  try {
  const size_t m = static_cast<size_t>(Rf_xlength(dist_sexp));
  const double* vd = REAL(dist_sexp);
  const int n_targets = Rf_length(avc_vec_sexp);
  const double* avc_targets = REAL(avc_vec_sexp);

  std::vector<HistBin> bins;
  if (m >= 4096) {
    double min_d = std::numeric_limits<double>::infinity();
    double max_d = -std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < m; ++k) {
      const double d = vd[k];
      if (d < min_d) min_d = d;
      if (d > max_d) max_d = d;
    }
    build_hist_bins(vd, m, min_d, max_d, bins);
  }

  std::vector<EvalPoint> history;
  history.reserve(64);

  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, n_targets));
  solve_c0_from_avc_core(vd, m, bins, history, avc_targets, n_targets, REAL(out_sexp));
  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_getc0_from_distmat(SEXP distmat_sexp, SEXP avc_vec_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);
  const int n_targets = Rf_length(avc_vec_sexp);
  const double* avc_targets = REAL(avc_vec_sexp);

  const size_t m = static_cast<size_t>(n) * (n - 1) / 2;
  std::vector<double> vd(m);
  double min_d = std::numeric_limits<double>::infinity();
  double max_d = -std::numeric_limits<double>::infinity();

  size_t pos = 0;
  for (int j = 0; j < n - 1; ++j) {
    const double* col = distmat + static_cast<size_t>(j) * n + (j + 1);
    const size_t len = static_cast<size_t>(n - 1 - j);
    std::memcpy(vd.data() + pos, col, len * sizeof(double));
    if (m >= 4096) {
      for (size_t k = 0; k < len; ++k) {
        const double d = col[k];
        if (d < min_d) min_d = d;
        if (d > max_d) max_d = d;
      }
    }
    pos += len;
  }

  std::vector<HistBin> bins;
  if (m >= 4096) {
    build_hist_bins(vd.data(), m, min_d, max_d, bins);
  }

  std::vector<EvalPoint> history;
  history.reserve(64);

  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, n_targets));
  solve_c0_from_avc_core(vd.data(), m, bins, history, avc_targets, n_targets, REAL(out_sexp));
  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_getrp(SEXP om_sexp, SEXP cv_sexp, SEXP gqx_sexp, SEXP gqw_sexp) {
  const char* error_message = nullptr;
  try {
  const int k = Rf_nrows(om_sexp);
  const double* om = REAL(om_sexp);
  const double cv = Rf_asReal(cv_sexp);
  const double* gqx = REAL(gqx_sexp);
  const double* gqw = REAL(gqw_sexp);
  const int n_gq = Rf_length(gqx_sexp);

  RpWorkspace ws;
  const double rp = eval_getrp_submatrix(om, k, k, cv, gqx, gqw, n_gq, ws);
  return Rf_ScalarReal(rp);
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_maxrp(SEXP oms_sexp, SEXP q_sexp, SEXP cv_sexp, SEXP gqx_sexp, SEXP gqw_sexp) {
  const char* error_message = nullptr;
  try {
  const int nc = Rf_length(oms_sexp);
  const int q = Rf_asInteger(q_sexp);
  const double cv = Rf_asReal(cv_sexp);
  const double* gqx = REAL(gqx_sexp);
  const double* gqw = REAL(gqw_sexp);
  const int n_gq = Rf_length(gqx_sexp);

  std::vector<const double*> oms_ptrs(nc);
  std::vector<int> oms_lds(nc);
  for (int m = 0; m < nc; ++m) {
    SEXP mat = VECTOR_ELT(oms_sexp, m);
    oms_ptrs[m] = REAL(mat);
    oms_lds[m] = Rf_nrows(mat);
  }

  RpWorkspace ws;
  double max_rp = 0.0;
  int max_idx = 0;
  eval_maxrp_core(oms_ptrs, oms_lds, q, cv, gqx, gqw, n_gq, ws, -1, 0.0, max_rp, max_idx);

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("max"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("i"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, Rf_ScalarReal(max_rp));
  SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarInteger(max_idx + 1));

  UNPROTECT(2);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_getcv(SEXP oms_sexp, SEXP q_sexp, SEXP level_sexp, SEXP gqx_sexp, SEXP gqw_sexp) {
  const char* error_message = nullptr;
  try {
  const int nc = Rf_length(oms_sexp);
  const int q = Rf_asInteger(q_sexp);
  const int n_lev = Rf_length(level_sexp);
  const double* levels = REAL(level_sexp);
  const double* gqx = REAL(gqx_sexp);
  const double* gqw = REAL(gqw_sexp);
  const int n_gq = Rf_length(gqx_sexp);

  std::vector<const double*> oms_ptrs(nc);
  std::vector<int> oms_lds(nc);
  for (int m = 0; m < nc; ++m) {
    SEXP mat = VECTOR_ELT(oms_sexp, m);
    oms_ptrs[m] = REAL(mat);
    oms_lds[m] = Rf_nrows(mat);
  }

  RpWorkspace ws;
  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, n_lev));
  double* out = REAL(out_sexp);
  for (int l = 0; l < n_lev; ++l) {
    out[l] = eval_getcv_core(oms_ptrs, oms_lds, q, levels[l], gqx, gqw, n_gq, ws);
  }

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_setfinalW(SEXP oms_sexp, SEXP w_sexp, SEXP qmax_sexp, SEXP gqx_sexp, SEXP gqw_sexp) {
  const char* error_message = nullptr;
  try {
  const int nc = Rf_length(oms_sexp);
  const int n = Rf_nrows(w_sexp);
  const double* W = REAL(w_sexp);
  const int qmax = Rf_asInteger(qmax_sexp);
  const double* gqx = REAL(gqx_sexp);
  const double* gqw = REAL(gqw_sexp);
  const int n_gq = Rf_length(gqx_sexp);

  std::vector<const double*> oms_ptrs(nc);
  std::vector<int> oms_lds(nc);
  for (int m = 0; m < nc; ++m) {
    SEXP mat = VECTOR_ELT(oms_sexp, m);
    oms_ptrs[m] = REAL(mat);
    oms_lds[m] = Rf_nrows(mat);
  }

  RpWorkspace ws;
  double best_len = std::numeric_limits<double>::infinity();
  double best_cv = 0.0;
  int q_opt = 1;

  for (int q = 1; q <= qmax; ++q) {
    const double cv_q = eval_getcv_core(oms_ptrs, oms_lds, q, 0.05, gqx, gqw, n_gq, ws);
    const double len_q = cv_q * std::tgamma(0.5 * (q + 1)) / (std::sqrt(static_cast<double>(q)) * std::tgamma(0.5 * q));
    if (q == 1 || len_q < best_len) {
      best_len = len_q;
      best_cv = cv_q;
      q_opt = q;
    }
  }

  const int ncol_out = q_opt + 1;
  SEXP w_out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, ncol_out));
  std::memcpy(REAL(w_out_sexp), W, static_cast<size_t>(n) * ncol_out * sizeof(double));

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("W"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("cv"));
  SET_STRING_ELT(names_sexp, 2, Rf_mkChar("q"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, w_out_sexp);
  SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarReal(best_cv));
  SET_VECTOR_ELT(out_sexp, 2, Rf_ScalarInteger(q_opt));

  UNPROTECT(3);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_get_oms(SEXP distmat_sexp, SEXP c0_sexp, SEXP cmax_sexp, SEXP w_sexp, SEXP cgridfac_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(distmat_sexp);
  const int k_cols = Rf_ncols(w_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double* W = REAL(w_sexp);
  const double c0 = Rf_asReal(c0_sexp);
  const double cmax = Rf_asReal(cmax_sexp);
  const double cgridfac = Rf_asReal(cgridfac_sexp);

  const int nc = std::max(2, static_cast<int>(std::ceil(std::log(cmax / c0) / std::log(cgridfac))));
  std::vector<double> c_vec(nc);
  double c_cur = c0;
  for (int m = 0; m < nc; ++m) {
    c_vec[m] = c_cur;
    c_cur *= cgridfac;
  }

  constexpr int BLK = 256;
  std::vector<double> SW_all(static_cast<size_t>(n) * k_cols * nc, 0.0);
  std::vector<double> D_blk(BLK * BLK);
  std::vector<double> S_blk(BLK * BLK);

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta_one = 1.0;
  const double beta_zero = 0.0;

  for (int j0 = 0; j0 < n; j0 += BLK) {
    const int bj = std::min(BLK, n - j0);
    for (int i0 = j0; i0 < n; i0 += BLK) {
      const int bi = std::min(BLK, n - i0);
      const int len = bi * bj;

      for (int j = 0; j < bj; ++j) {
        const double* d_col = distmat + static_cast<size_t>(j0 + j) * n + i0;
        std::memcpy(D_blk.data() + static_cast<size_t>(j) * bi, d_col, bi * sizeof(double));
      }

      const double* W_i = W + i0;
      const double* W_j = W + j0;
      const double* d_ptr = D_blk.data();
      double* s_ptr = S_blk.data();

      for (int m = 0; m < nc; ++m) {
        const double neg_c = -c_vec[m];
#ifdef __APPLE__
        for (int idx = 0; idx < len; ++idx) {
          s_ptr[idx] = neg_c * d_ptr[idx];
        }
        vvexp(s_ptr, s_ptr, &len);
#else
        for (int idx = 0; idx < len; ++idx) {
          s_ptr[idx] = std::exp(neg_c * d_ptr[idx]);
        }
#endif
        if (i0 == j0) {
          for (int d = 0; d < bi; ++d) {
            if (d_ptr[d + d * bi] == 0.0) s_ptr[d + d * bi] = 1.0;
          }
        }

        double* SW_m = SW_all.data() + static_cast<size_t>(m) * n * k_cols;
        if (i0 == j0) {
          F77_CALL(dgemm)(&trans_n, &trans_n, &bi, &k_cols, &bj, &alpha,
                          s_ptr, &bi, W_i, &n, &beta_one, SW_m + i0, &n FCONE FCONE);
        } else {
          F77_CALL(dgemm)(&trans_n, &trans_n, &bi, &k_cols, &bj, &alpha,
                          s_ptr, &bi, W_j, &n, &beta_one, SW_m + i0, &n FCONE FCONE);
          F77_CALL(dgemm)(&trans_t, &trans_n, &bj, &k_cols, &bi, &alpha,
                          s_ptr, &bi, W_i, &n, &beta_one, SW_m + j0, &n FCONE FCONE);
        }
      }
    }
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, nc));
  for (int m = 0; m < nc; ++m) {
    SEXP om_sexp = PROTECT(Rf_allocMatrix(REALSXP, k_cols, k_cols));
    const double* SW_m = SW_all.data() + static_cast<size_t>(m) * n * k_cols;
    F77_CALL(dgemm)(&trans_t, &trans_n, &k_cols, &k_cols, &n, &alpha,
                    W, &n, SW_m, &n, &beta_zero, REAL(om_sexp), &k_cols FCONE FCONE);
    SET_VECTOR_ELT(out_sexp, m, om_sexp);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

// Compute .getOms for a list of B matrices Wx[[1..B]] (each n x k_cols) in a single pass over distmat
SEXP cpp_get_oms_batch(SEXP distmat_sexp, SEXP c0_sexp, SEXP cmax_sexp, SEXP w_list_sexp, SEXP cgridfac_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(distmat_sexp);
  const int n_batch = Rf_length(w_list_sexp);
  if (n_batch == 0) {
    return Rf_allocVector(VECSXP, 0);
  }
  if (n_batch == 1) {
    SEXP res1 = PROTECT(cpp_get_oms(distmat_sexp, c0_sexp, cmax_sexp, VECTOR_ELT(w_list_sexp, 0), cgridfac_sexp));
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(out, 0, res1);
    UNPROTECT(2);
    return out;
  }

  const int k_cols = Rf_ncols(VECTOR_ELT(w_list_sexp, 0));
  const int total_cols = n_batch * k_cols;
  const double* distmat = REAL(distmat_sexp);
  const double c0 = Rf_asReal(c0_sexp);
  const double cmax = Rf_asReal(cmax_sexp);
  const double cgridfac = Rf_asReal(cgridfac_sexp);

  std::vector<double> W_stack(static_cast<size_t>(n) * total_cols);
  for (int b = 0; b < n_batch; ++b) {
    const double* W_b = REAL(VECTOR_ELT(w_list_sexp, b));
    std::memcpy(W_stack.data() + static_cast<size_t>(b) * n * k_cols,
                W_b,
                static_cast<size_t>(n) * k_cols * sizeof(double));
  }

  const int nc = std::max(2, static_cast<int>(std::ceil(std::log(cmax / c0) / std::log(cgridfac))));
  std::vector<double> c_vec(nc);
  double c_cur = c0;
  for (int m = 0; m < nc; ++m) {
    c_vec[m] = c_cur;
    c_cur *= cgridfac;
  }

  constexpr int BLK = 256;
  std::vector<double> SW_all(static_cast<size_t>(n) * total_cols * nc, 0.0);
  std::vector<double> D_blk(BLK * BLK);
  std::vector<double> S_blk(BLK * BLK);

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta_one = 1.0;
  const double beta_zero = 0.0;
  const double* W = W_stack.data();

  for (int j0 = 0; j0 < n; j0 += BLK) {
    const int bj = std::min(BLK, n - j0);
    for (int i0 = j0; i0 < n; i0 += BLK) {
      const int bi = std::min(BLK, n - i0);
      const int len = bi * bj;

      for (int j = 0; j < bj; ++j) {
        const double* d_col = distmat + static_cast<size_t>(j0 + j) * n + i0;
        std::memcpy(D_blk.data() + static_cast<size_t>(j) * bi, d_col, bi * sizeof(double));
      }

      const double* W_i = W + i0;
      const double* W_j = W + j0;
      const double* d_ptr = D_blk.data();
      double* s_ptr = S_blk.data();

      for (int m = 0; m < nc; ++m) {
        const double neg_c = -c_vec[m];
#ifdef __APPLE__
        for (int idx = 0; idx < len; ++idx) {
          s_ptr[idx] = neg_c * d_ptr[idx];
        }
        vvexp(s_ptr, s_ptr, &len);
#else
        for (int idx = 0; idx < len; ++idx) {
          s_ptr[idx] = std::exp(neg_c * d_ptr[idx]);
        }
#endif
        if (i0 == j0) {
          for (int d = 0; d < bi; ++d) {
            if (d_ptr[d + d * bi] == 0.0) s_ptr[d + d * bi] = 1.0;
          }
        }

        double* SW_m = SW_all.data() + static_cast<size_t>(m) * n * total_cols;
        if (i0 == j0) {
          F77_CALL(dgemm)(&trans_n, &trans_n, &bi, &total_cols, &bj, &alpha,
                          s_ptr, &bi, W_i, &n, &beta_one, SW_m + i0, &n FCONE FCONE);
        } else {
          F77_CALL(dgemm)(&trans_n, &trans_n, &bi, &total_cols, &bj, &alpha,
                          s_ptr, &bi, W_j, &n, &beta_one, SW_m + i0, &n FCONE FCONE);
          F77_CALL(dgemm)(&trans_t, &trans_n, &bj, &total_cols, &bi, &alpha,
                          s_ptr, &bi, W_i, &n, &beta_one, SW_m + j0, &n FCONE FCONE);
        }
      }
    }
  }

  SEXP batch_out = PROTECT(Rf_allocVector(VECSXP, n_batch));
  for (int b = 0; b < n_batch; ++b) {
    SEXP oms_b = PROTECT(Rf_allocVector(VECSXP, nc));
    const double* W_b = W + static_cast<size_t>(b) * n * k_cols;
    for (int m = 0; m < nc; ++m) {
      SEXP om_sexp = PROTECT(Rf_allocMatrix(REALSXP, k_cols, k_cols));
      const double* SW_bm = SW_all.data() + static_cast<size_t>(m) * n * total_cols + static_cast<size_t>(b) * n * k_cols;
      F77_CALL(dgemm)(&trans_t, &trans_n, &k_cols, &k_cols, &n, &alpha,
                      W_b, &n, SW_bm, &n, &beta_zero, REAL(om_sexp), &k_cols FCONE FCONE);
      SET_VECTOR_ELT(oms_b, m, om_sexp);
      UNPROTECT(1);
    }
    SET_VECTOR_ELT(batch_out, b, oms_b);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return batch_out;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_raninds(SEXP n_sexp, SEXP capM_sexp, SEXP random_t_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_asInteger(n_sexp);
  const int capM = Rf_asInteger(capM_sexp);
  uint64_t rt = static_cast<uint64_t>(Rf_asReal(random_t_sexp));

  SEXP v_sexp = PROTECT(Rf_allocVector(INTSXP, capM + 1));
  int* v = INTEGER(v_sexp);

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
  int j = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));

  for (int i = 0; i <= capM; ++i) {
    v[i] = j + 1;
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
    const int step = static_cast<int>(std::floor(static_cast<double>(n - 1) * (static_cast<double>(rt) * INV_2POW32)));
    j = (j + 1 + step) % n;
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("indices"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("state"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, v_sexp);
  SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarReal(static_cast<double>(rt)));

  UNPROTECT(3);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_jumble_s(SEXP s_sexp, SEXP m_sexp, SEXP random_t_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const int m = Rf_asInteger(m_sexp);
  uint64_t rt = static_cast<uint64_t>(Rf_asReal(random_t_sexp));

  SEXP r_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, d));
  double* r = REAL(r_sexp);
  std::memcpy(r, REAL(s_sexp), static_cast<size_t>(n) * d * sizeof(double));

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  for (int i = 0; i < m; ++i) {
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
    const int j = i + static_cast<int>(std::floor(static_cast<double>(n - i) * (static_cast<double>(rt) * INV_2POW32)));
    if (j != i) {
      for (int k = 0; k < d; ++k) {
        std::swap(r[i + static_cast<size_t>(k) * n], r[j + static_cast<size_t>(k) * n]);
      }
    }
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("coords"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("state"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, r_sexp);
  SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarReal(static_cast<double>(rt)));

  UNPROTECT(3);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_lnset_wc0_subsamples(SEXP s_sexp,
                              SEXP capN_sexp,
                              SEXP m_sexp,
                              SEXP random_t_sexp,
                              SEXP latlong_sexp,
                              SEXP avc0_sexp,
                              SEXP minavc_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const int capN = Rf_asInteger(capN_sexp);
  const int m = std::min(Rf_asInteger(m_sexp), n);
  uint64_t rt = static_cast<uint64_t>(Rf_asReal(random_t_sexp));
  const bool latlong = Rf_asLogical(latlong_sexp) == TRUE;
  const double avc0 = Rf_asReal(avc0_sexp);
  const double minavc = Rf_asReal(minavc_sexp);

  std::vector<double> r(static_cast<size_t>(n) * d);
  std::memcpy(r.data(), REAL(s_sexp), static_cast<size_t>(n) * d * sizeof(double));

  const size_t block_len = static_cast<size_t>(m) * (m - 1) / 2;
  const size_t total_len = static_cast<size_t>(capN) * block_len;
  std::vector<double> distvec(total_len);
  double min_d = std::numeric_limits<double>::infinity();
  double max_d = -std::numeric_limits<double>::infinity();

  SEXP coords_list = PROTECT(Rf_allocVector(VECSXP, capN));
  SEXP distmat_list = PROTECT(Rf_allocVector(VECSXP, capN));

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  for (int idx = 0; idx < capN; ++idx) {
    for (int i = 0; i < m; ++i) {
      rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
      const int j = i + static_cast<int>(std::floor(static_cast<double>(n - i) * (static_cast<double>(rt) * INV_2POW32)));
      if (j != i) {
        for (int k = 0; k < d; ++k) {
          std::swap(r[i + static_cast<size_t>(k) * n], r[j + static_cast<size_t>(k) * n]);
        }
      }
    }

    SEXP sub_coords = PROTECT(Rf_allocMatrix(REALSXP, m, d));
    double* sc_ptr = REAL(sub_coords);
    for (int k = 0; k < d; ++k) {
      std::memcpy(sc_ptr + static_cast<size_t>(k) * m,
                  r.data() + static_cast<size_t>(k) * n,
                  static_cast<size_t>(m) * sizeof(double));
    }
    SET_VECTOR_ELT(coords_list, idx, sub_coords);

    SEXP sub_dist = PROTECT(Rf_allocMatrix(REALSXP, m, m));
    double* sd_ptr = REAL(sub_dist);
    compute_distmat_raw(sc_ptr, m, d, latlong, sd_ptr);
    SET_VECTOR_ELT(distmat_list, idx, sub_dist);
    UNPROTECT(2);

    size_t pos = static_cast<size_t>(idx) * block_len;
    for (int j = 0; j < m - 1; ++j) {
      const double* col = sd_ptr + static_cast<size_t>(j) * m + (j + 1);
      const size_t len = static_cast<size_t>(m - 1 - j);
      std::memcpy(distvec.data() + pos, col, len * sizeof(double));
      for (size_t k = 0; k < len; ++k) {
        const double v = col[k];
        if (v < min_d) min_d = v;
        if (v > max_d) max_d = v;
      }
      pos += len;
    }
  }

  std::vector<HistBin> bins;
  build_hist_bins(distvec.data(), total_len, min_d, max_d, bins);
  std::vector<EvalPoint> history;
  history.reserve(64);

  const double floor_corr = static_cast<double>(std::count(distvec.begin(), distvec.end(), 0.0)) /
      static_cast<double>(total_len);
  const double targets[2] = {avc0, floor_corr + (1.0 - floor_corr) * minavc};
  double c_out[2] = {0.0, 0.0};
  solve_c0_from_avc_core(distvec.data(), total_len, bins, history, targets, 2, c_out);

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("coords_list"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("distmat_list"));
  SET_STRING_ELT(names_sexp, 2, Rf_mkChar("c0"));
  SET_STRING_ELT(names_sexp, 3, Rf_mkChar("cmax"));
  SET_STRING_ELT(names_sexp, 4, Rf_mkChar("random_t"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, coords_list);
  SET_VECTOR_ELT(out_sexp, 1, distmat_list);
  SET_VECTOR_ELT(out_sexp, 2, Rf_ScalarReal(c_out[0]));
  SET_VECTOR_ELT(out_sexp, 3, Rf_ScalarReal(c_out[1]));
  SET_VECTOR_ELT(out_sexp, 4, Rf_ScalarReal(static_cast<double>(rt)));

  UNPROTECT(4);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

// Compute Wall (n x (capN * qmax)) from S (n x d), coords_list (capN of m x d), W0_list (capN of m x qmax)
SEXP cpp_lnset_wc0_wall(SEXP s_sexp,
                        SEXP coords_list_sexp,
                        SEXP w0_list_sexp,
                        SEXP c0_sexp,
                        SEXP qmax_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const double* S = REAL(s_sexp);
  const int capN = Rf_length(coords_list_sexp);
  const int qmax = Rf_asInteger(qmax_sexp);
  const double neg_c0 = -Rf_asReal(c0_sexp);

  const int m = Rf_nrows(VECTOR_ELT(coords_list_sexp, 0));
  SEXP wall_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, capN * qmax));
  double* Wall = REAL(wall_sexp);

  std::vector<double> K_mat(static_cast<size_t>(n) * m);
  std::vector<double> Wx(static_cast<size_t>(n) * qmax);

  const char trans_n = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;
  const double inv_n = 1.0 / static_cast<double>(n);

  for (int idx = 0; idx < capN; ++idx) {
    const double* sub_coords = REAL(VECTOR_ELT(coords_list_sexp, idx));
    const double* W0 = REAL(VECTOR_ELT(w0_list_sexp, idx));

    if (d == 2) {
      const double* s0 = S;
      const double* s1 = S + n;
      const double* c0_col = sub_coords;
      const double* c1_col = sub_coords + m;
      for (int j = 0; j < m; ++j) {
        const double cx = c0_col[j];
        const double cy = c1_col[j];
        double* k_col = K_mat.data() + static_cast<size_t>(j) * n;
        for (int i = 0; i < n; ++i) {
          const double dx = s0[i] - cx;
          const double dy = s1[i] - cy;
          k_col[i] = neg_c0 * std::sqrt(dx * dx + dy * dy);
        }
#ifdef __APPLE__
        vvexp(k_col, k_col, &n);
#else
        for (int i = 0; i < n; ++i) {
          k_col[i] = std::exp(k_col[i]);
        }
#endif
      }
    } else {
      for (int j = 0; j < m; ++j) {
        double* k_col = K_mat.data() + static_cast<size_t>(j) * n;
        for (int i = 0; i < n; ++i) {
          double sum_sq = 0.0;
          for (int col_d = 0; col_d < d; ++col_d) {
            const double diff = S[i + static_cast<size_t>(col_d) * n] - sub_coords[j + static_cast<size_t>(col_d) * m];
            sum_sq += diff * diff;
          }
          k_col[i] = neg_c0 * std::sqrt(sum_sq);
        }
#ifdef __APPLE__
        vvexp(k_col, k_col, &n);
#else
        for (int i = 0; i < n; ++i) {
          k_col[i] = std::exp(k_col[i]);
        }
#endif
      }
    }

    F77_CALL(dgemm)(&trans_n, &trans_n, &n, &qmax, &m, &alpha,
                    K_mat.data(), &n, W0, &m, &beta, Wx.data(), &n FCONE FCONE);

    for (int q = 0; q < qmax; ++q) {
      double* wx_col = Wx.data() + static_cast<size_t>(q) * n;
      double csum = 0.0;
      for (int i = 0; i < n; ++i) {
        csum += wx_col[i];
      }
      const double cmean = csum * inv_n;
      double norm_sq = 0.0;
      for (int i = 0; i < n; ++i) {
        const double v = wx_col[i] - cmean;
        wx_col[i] = v;
        norm_sq += v * v;
      }
      double norm = std::sqrt(norm_sq);
      if (!std::isfinite(norm) || norm == 0.0) norm = 1.0;
      const double inv_norm = 1.0 / norm;

      double* dst_col = Wall + static_cast<size_t>(q * capN + idx) * n;
      for (int i = 0; i < n; ++i) {
        dst_col[i] = wx_col[i] * inv_norm;
      }
    }
  }

  UNPROTECT(1);
  return wall_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_lnset_wc0(SEXP s_sexp,
                   SEXP avc0_sexp,
                   SEXP qmax_sexp,
                   SEXP minavc_sexp,
                   SEXP latlong_sexp,
                   SEXP capN_sexp,
                   SEXP m_sexp,
                   SEXP random_t_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const double* S = REAL(s_sexp);
  const double avc0 = Rf_asReal(avc0_sexp);
  const int qmax = Rf_asInteger(qmax_sexp);
  const double minavc = Rf_asReal(minavc_sexp);
  const bool latlong = Rf_asLogical(latlong_sexp) == TRUE;
  const int capN = Rf_asInteger(capN_sexp);
  const int m = std::min(Rf_asInteger(m_sexp), n);
  uint64_t rt = static_cast<uint64_t>(Rf_asReal(random_t_sexp));

  std::vector<double> r(static_cast<size_t>(n) * d);
  std::memcpy(r.data(), S, static_cast<size_t>(n) * d * sizeof(double));

  const size_t block_len = static_cast<size_t>(m) * (m - 1) / 2;
  const size_t total_len = static_cast<size_t>(capN) * block_len;
  std::vector<double> distvec(total_len);
  double min_d = std::numeric_limits<double>::infinity();
  double max_d = -std::numeric_limits<double>::infinity();

  std::vector<double> sub_coords_all(static_cast<size_t>(capN) * m * d);
  std::vector<double> sub_dist_all(static_cast<size_t>(capN) * m * m);

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  for (int idx = 0; idx < capN; ++idx) {
    for (int i = 0; i < m; ++i) {
      rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
      const int j = i + static_cast<int>(std::floor(static_cast<double>(n - i) * (static_cast<double>(rt) * INV_2POW32)));
      if (j != i) {
        for (int k = 0; k < d; ++k) {
          std::swap(r[i + static_cast<size_t>(k) * n], r[j + static_cast<size_t>(k) * n]);
        }
      }
    }

    double* sc_ptr = sub_coords_all.data() + static_cast<size_t>(idx) * m * d;
    for (int k = 0; k < d; ++k) {
      std::memcpy(sc_ptr + static_cast<size_t>(k) * m,
                  r.data() + static_cast<size_t>(k) * n,
                  static_cast<size_t>(m) * sizeof(double));
    }

    double* sd_ptr = sub_dist_all.data() + static_cast<size_t>(idx) * m * m;
    compute_distmat_raw(sc_ptr, m, d, latlong, sd_ptr);

    size_t pos = static_cast<size_t>(idx) * block_len;
    for (int j = 0; j < m - 1; ++j) {
      const double* col = sd_ptr + static_cast<size_t>(j) * m + (j + 1);
      const size_t len = static_cast<size_t>(m - 1 - j);
      std::memcpy(distvec.data() + pos, col, len * sizeof(double));
      for (size_t k = 0; k < len; ++k) {
        const double v = col[k];
        if (v < min_d) min_d = v;
        if (v > max_d) max_d = v;
      }
      pos += len;
    }
  }

  std::vector<HistBin> bins;
  build_hist_bins(distvec.data(), total_len, min_d, max_d, bins);
  std::vector<EvalPoint> history;
  history.reserve(64);

  const double floor_corr = static_cast<double>(std::count(distvec.begin(), distvec.end(), 0.0)) /
      static_cast<double>(total_len);
  const double targets[2] = {avc0, floor_corr + (1.0 - floor_corr) * minavc};
  double c_out[2] = {0.0, 0.0};
  solve_c0_from_avc_core(distvec.data(), total_len, bins, history, targets, 2, c_out);
  const double c0 = c_out[0];
  const double cmax = c_out[1];
  const double neg_c0 = -c0;

  std::vector<double> unit(latlong ? static_cast<size_t>(n) * 3 : 0);
  std::vector<double> sub_unit(latlong ? static_cast<size_t>(m) * 3 : 0);
  if (latlong) spherical_positions(S, n, unit.data());

  const int total_wall_cols = capN * qmax;
  std::vector<double> Wall(static_cast<size_t>(n) * total_wall_cols);
  std::vector<double> sig_d(static_cast<size_t>(m) * m);
  std::vector<double> W0(static_cast<size_t>(m) * qmax);
  constexpr int KERNEL_ROWS = 256;
  std::vector<double> K_mat(static_cast<size_t>(std::min(n, KERNEL_ROWS)) * m);

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta = 0.0;
  const double inv_n = 1.0 / static_cast<double>(n);

  for (int idx = 0; idx < capN; ++idx) {
    const double* sub_coords = sub_coords_all.data() + static_cast<size_t>(idx) * m * d;
    const double* sub_dist = sub_dist_all.data() + static_cast<size_t>(idx) * m * m;

    exp_lower_triangle_inplace(sub_dist, sig_d.data(), m, neg_c0);
    symmetrize_lower_to_upper_blocked(sig_d.data(), m);
    demean_rect_inplace(sig_d.data(), m, m);
    compute_top_k_evecs_sym(sig_d.data(), m, qmax, W0.data());

    if (latlong) spherical_positions(sub_coords, m, sub_unit.data());

    // Write projections directly into their strided Wall columns.
    const int wall_stride = n * capN;
    for (int row0 = 0; row0 < n; row0 += KERNEL_ROWS) {
      const int rows = std::min(KERNEL_ROWS, n - row0);
      compute_cross_kernel(S, n, d, sub_coords, m, latlong, neg_c0,
                           unit.data(), sub_unit.data(), K_mat.data(), row0, rows);
      F77_CALL(dgemm)(&trans_n, &trans_n, &rows, &qmax, &m, &alpha,
                      K_mat.data(), &rows, W0.data(), &m, &beta,
                      Wall.data() + static_cast<size_t>(idx) * n + row0, &wall_stride FCONE FCONE);
    }

    for (int q = 0; q < qmax; ++q) {
      double* wx_col = Wall.data() + static_cast<size_t>(q * capN + idx) * n;
      double csum = 0.0;
      for (int i = 0; i < n; ++i) {
        csum += wx_col[i];
      }
      const double cmean = csum * inv_n;
      double norm_sq = 0.0;
      for (int i = 0; i < n; ++i) {
        const double v = wx_col[i] - cmean;
        wx_col[i] = v;
        norm_sq += v * v;
      }
      double norm = std::sqrt(norm_sq);
      if (!std::isfinite(norm) || norm == 0.0) norm = 1.0;
      const double inv_norm = 1.0 / norm;

      for (int i = 0; i < n; ++i) wx_col[i] *= inv_norm;
    }
  }

  SEXP w_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, qmax + 1));
  double* W_out = REAL(w_sexp);
  const double inv_sqrt_n = 1.0 / std::sqrt(static_cast<double>(n));
  for (int i = 0; i < n; ++i) {
    W_out[i] = inv_sqrt_n;
  }

  std::vector<double> G_mat(static_cast<size_t>(total_wall_cols) * total_wall_cols);
  std::vector<double> evec(total_wall_cols);
  std::vector<double> proj(total_wall_cols);
  const int inc1 = 1;
  const double neg_one = -1.0;

  for (int i = 1; i <= qmax; ++i) {
    const int ncol_wx = capN * i;
    F77_CALL(dgemm)(&trans_t, &trans_n, &ncol_wx, &ncol_wx, &n, &alpha,
                    Wall.data(), &n, Wall.data(), &n, &beta, G_mat.data(), &ncol_wx FCONE FCONE);
    compute_top_k_evecs_sym(G_mat.data(), ncol_wx, 1, evec.data());

    double* w_col = W_out + static_cast<size_t>(i) * n;
    F77_CALL(dgemv)(&trans_n, &n, &ncol_wx, &alpha,
                    Wall.data(), &n, evec.data(), &inc1, &beta, w_col, &inc1 FCONE);

    double norm_sq = 0.0;
    for (int r_idx = 0; r_idx < n; ++r_idx) {
      norm_sq += w_col[r_idx] * w_col[r_idx];
    }
    const double inv_norm = 1.0 / std::sqrt(norm_sq);
    for (int r_idx = 0; r_idx < n; ++r_idx) {
      w_col[r_idx] *= inv_norm;
    }

    if (i < qmax) {
      F77_CALL(dgemv)(&trans_t, &n, &total_wall_cols, &alpha,
                      Wall.data(), &n, w_col, &inc1, &beta, proj.data(), &inc1 FCONE);
      F77_CALL(dger)(&n, &total_wall_cols, &neg_one,
                     w_col, &inc1, proj.data(), &inc1, Wall.data(), &n);
    }
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("W"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("c0"));
  SET_STRING_ELT(names_sexp, 2, Rf_mkChar("cmax"));
  SET_STRING_ELT(names_sexp, 3, Rf_mkChar("random_t"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, w_sexp);
  SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarReal(c0));
  SET_VECTOR_ELT(out_sexp, 2, Rf_ScalarReal(cmax));
  SET_VECTOR_ELT(out_sexp, 3, Rf_ScalarReal(static_cast<double>(rt)));

  UNPROTECT(3);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

SEXP cpp_lnget_oms(SEXP s_sexp,
                   SEXP c0_sexp,
                   SEXP cmax_sexp,
                   SEXP w_sexp,
                   SEXP cgridfac_sexp,
                   SEXP capM_sexp,
                   SEXP random_t_sexp,
                   SEXP latlong_sexp) {
  const char* error_message = nullptr;
  try {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const int k_cols = Rf_ncols(w_sexp);
  const double* S = REAL(s_sexp);
  const double* W = REAL(w_sexp);
  const double c0 = Rf_asReal(c0_sexp);
  const double cmax = Rf_asReal(cmax_sexp);
  const double cgridfac = Rf_asReal(cgridfac_sexp);
  const int capM = Rf_asInteger(capM_sexp);
  uint64_t rt = static_cast<uint64_t>(Rf_asReal(random_t_sexp));
  const bool latlong = Rf_asLogical(latlong_sexp) == TRUE;

  const int nc = std::max(2, static_cast<int>(std::ceil(std::log(cmax / c0) / std::log(cgridfac))));
  std::vector<double> c_vec(nc);
  double c_cur = c0;
  for (int m = 0; m < nc; ++m) {
    c_vec[m] = c_cur;
    c_cur *= cgridfac;
  }

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  const double scale = 0.5 * (static_cast<double>(n) * (n - 1) / static_cast<double>(capM));
  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta_one = 1.0;
  const double beta_zero = 0.0;

  // The diagonal contribution is W'W; conditional weights are not orthonormal.
  std::vector<double> gram(static_cast<size_t>(k_cols) * k_cols);
  F77_CALL(dgemm)(&trans_t, &trans_n, &k_cols, &k_cols, &n, &alpha,
                  W, &n, W, &n, &beta_zero, gram.data(), &k_cols FCONE FCONE);

  // When n <= 1200, n*(n-1)/2 <= 719,400 < capM = 1,000,000, so accumulating
  // unordered pair counts into an n x n matrix reduces 1,000,000 pair evaluations to n x n!
  if (n <= 1200) {
    std::vector<int> pair_count(static_cast<size_t>(n) * n, 0);
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
    int j_prev = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));
    for (int m = 0; m < capM; ++m) {
      rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
      const int step = static_cast<int>(std::floor(static_cast<double>(n - 1) * (static_cast<double>(rt) * INV_2POW32)));
      const int j_next = (j_prev + 1 + step) % n;
      const int u = std::min(j_prev, j_next);
      const int v = std::max(j_prev, j_next);
      pair_count[v + static_cast<size_t>(u) * n]++;
      j_prev = j_next;
    }
    // One extra LCG step because .raninds generates capM + 1 indices (looping i = 1 .. capM + 1)
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;

    std::vector<double> distmat(static_cast<size_t>(n) * n);
    compute_distmat_raw(S, n, d, latlong, distmat.data());

    std::vector<double> A_mat(static_cast<size_t>(n) * n, 0.0);
    std::vector<double> AW(static_cast<size_t>(n) * k_cols);

    SEXP oms_sexp = PROTECT(Rf_allocVector(VECSXP, nc));
    for (int r = 0; r < nc; ++r) {
      const double neg_c = -c_vec[r];
      exp_lower_triangle_inplace(distmat.data(), A_mat.data(), n, neg_c);
      for (int u = 0; u < n; ++u) {
        A_mat[u + static_cast<size_t>(u) * n] = 0.0;
        double* a_col = A_mat.data() + static_cast<size_t>(u) * n;
        const int* cnt_col = pair_count.data() + static_cast<size_t>(u) * n;
        for (int v = u + 1; v < n; ++v) {
          a_col[v] *= static_cast<double>(cnt_col[v]);
        }
      }
      symmetrize_lower_to_upper_blocked(A_mat.data(), n);

      F77_CALL(dgemm)(&trans_n, &trans_n, &n, &k_cols, &n, &alpha,
                      A_mat.data(), &n, W, &n, &beta_zero, AW.data(), &n FCONE FCONE);

      SEXP om_sexp = PROTECT(Rf_allocMatrix(REALSXP, k_cols, k_cols));
      double* om = REAL(om_sexp);
      F77_CALL(dgemm)(&trans_t, &trans_n, &k_cols, &k_cols, &n, &scale,
                      W, &n, AW.data(), &n, &beta_zero, om, &k_cols FCONE FCONE);
      for (size_t a = 0; a < gram.size(); ++a) om[a] += gram[a];
      SET_VECTOR_ELT(oms_sexp, r, om_sexp);
      UNPROTECT(1);
    }

    SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names_sexp, 0, Rf_mkChar("Oms"));
    SET_STRING_ELT(names_sexp, 1, Rf_mkChar("state"));
    Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
    SET_VECTOR_ELT(out_sexp, 0, oms_sexp);
    SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarReal(static_cast<double>(rt)));

    UNPROTECT(3);
    return out_sexp;
  }

  // Large-n cache-tiled streaming path (for n > 1200)
  std::vector<double> lon_rad;
  std::vector<double> lat_rad;
  std::vector<double> cos_lat;
  if (latlong) {
    const double deg2rad = M_PI / 180.0;
    lon_rad.resize(n);
    lat_rad.resize(n);
    cos_lat.resize(n);
    const double* lon = S;
    const double* lat = S + n;
    for (int i = 0; i < n; ++i) {
      lon_rad[i] = lon[i] * deg2rad;
      lat_rad[i] = lat[i] * deg2rad;
      cos_lat[i] = std::cos(lat_rad[i]);
    }
  }

  constexpr int TILE = 4096;
  int idx1_tile[TILE];
  int idx2_tile[TILE];
  double dist_tile[TILE];
  double cd_tile[TILE];
#ifdef __APPLE__
  double dlat_tile[TILE];
  double dlon_tile[TILE];
#endif
  std::vector<double> W1_tile(static_cast<size_t>(TILE) * k_cols);
  std::vector<double> W2_tile(static_cast<size_t>(TILE) * k_cols);
  std::vector<double> W2_scaled(static_cast<size_t>(TILE) * k_cols);
  std::vector<double> acc_oms(static_cast<size_t>(k_cols) * k_cols * nc, 0.0);

  rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
  int j_prev = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));
  const double inv_pi = 1.0 / M_PI;

  for (int m0 = 0; m0 < capM; m0 += TILE) {
    const int b = std::min(TILE, capM - m0);
    for (int t = 0; t < b; ++t) {
      idx1_tile[t] = j_prev;
      rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
      const int step = static_cast<int>(std::floor(static_cast<double>(n - 1) * (static_cast<double>(rt) * INV_2POW32)));
      j_prev = (j_prev + 1 + step) % n;
      idx2_tile[t] = j_prev;
    }

    if (latlong) {
#ifdef __APPLE__
      for (int t = 0; t < b; ++t) {
        const int u = idx1_tile[t];
        const int v = idx2_tile[t];
        dlat_tile[t] = 0.5 * (lat_rad[u] - lat_rad[v]);
        dlon_tile[t] = 0.5 * (lon_rad[u] - lon_rad[v]);
      }
      vvsin(dlat_tile, dlat_tile, &b);
      vvsin(dlon_tile, dlon_tile, &b);
      for (int t = 0; t < b; ++t) {
        const int u = idx1_tile[t];
        const int v = idx2_tile[t];
        double inside = dlat_tile[t] * dlat_tile[t] + cos_lat[u] * cos_lat[v] * (dlon_tile[t] * dlon_tile[t]);
        if (inside > 1.0) inside = 1.0;
        if (inside < 0.0) inside = 0.0;
        dist_tile[t] = inside;
      }
      vvsqrt(dist_tile, dist_tile, &b);
      vvasin(dist_tile, dist_tile, &b);
      for (int t = 0; t < b; ++t) {
        dist_tile[t] *= inv_pi;
      }
#else
      for (int t = 0; t < b; ++t) {
        const int u = idx1_tile[t];
        const int v = idx2_tile[t];
        const double dlat = 0.5 * (lat_rad[u] - lat_rad[v]);
        const double dlon = 0.5 * (lon_rad[u] - lon_rad[v]);
        const double sdlat = std::sin(dlat);
        const double sdlon = std::sin(dlon);
        double inside = sdlat * sdlat + cos_lat[u] * cos_lat[v] * (sdlon * sdlon);
        if (inside > 1.0) inside = 1.0;
        if (inside < 0.0) inside = 0.0;
        dist_tile[t] = std::asin(std::sqrt(inside)) * inv_pi;
      }
#endif
    } else if (d == 2) {
      const double* s0 = S;
      const double* s1 = S + n;
      for (int t = 0; t < b; ++t) {
        const int u = idx1_tile[t];
        const int v = idx2_tile[t];
        const double dx = s0[u] - s0[v];
        const double dy = s1[u] - s1[v];
        dist_tile[t] = std::sqrt(dx * dx + dy * dy);
      }
    } else {
      for (int t = 0; t < b; ++t) {
        const int u = idx1_tile[t];
        const int v = idx2_tile[t];
        double sum_sq = 0.0;
        for (int col_d = 0; col_d < d; ++col_d) {
          const double diff = S[u + static_cast<size_t>(col_d) * n] - S[v + static_cast<size_t>(col_d) * n];
          sum_sq += diff * diff;
        }
        dist_tile[t] = std::sqrt(sum_sq);
      }
    }

    for (int col = 0; col < k_cols; ++col) {
      const double* w_col = W + static_cast<size_t>(col) * n;
      double* w1_col = W1_tile.data() + static_cast<size_t>(col) * b;
      double* w2_col = W2_tile.data() + static_cast<size_t>(col) * b;
      for (int t = 0; t < b; ++t) {
        w1_col[t] = w_col[idx1_tile[t]];
        w2_col[t] = w_col[idx2_tile[t]];
      }
    }

    for (int r = 0; r < nc; ++r) {
      const double neg_c = -c_vec[r];
#ifdef __APPLE__
      for (int t = 0; t < b; ++t) {
        cd_tile[t] = neg_c * dist_tile[t];
      }
      vvexp(cd_tile, cd_tile, &b);
#else
      for (int t = 0; t < b; ++t) {
        cd_tile[t] = std::exp(neg_c * dist_tile[t]);
      }
#endif
      for (int col = 0; col < k_cols; ++col) {
        const double* w2_col = W2_tile.data() + static_cast<size_t>(col) * b;
        double* w2s_col = W2_scaled.data() + static_cast<size_t>(col) * b;
        for (int t = 0; t < b; ++t) {
          w2s_col[t] = cd_tile[t] * w2_col[t];
        }
      }

      double* acc_r = acc_oms.data() + static_cast<size_t>(r) * k_cols * k_cols;
      F77_CALL(dgemm)(&trans_t, &trans_n, &k_cols, &k_cols, &b, &alpha,
                      W1_tile.data(), &b, W2_scaled.data(), &b, &beta_one, acc_r, &k_cols FCONE FCONE);
    }
  }
  // One final LCG step matching .raninds(n, capM, random_t) loop of capM + 1 iterations
  rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;

  SEXP oms_sexp = PROTECT(Rf_allocVector(VECSXP, nc));
  for (int r = 0; r < nc; ++r) {
    SEXP om_sexp = PROTECT(Rf_allocMatrix(REALSXP, k_cols, k_cols));
    double* om = REAL(om_sexp);
    const double* acc_r = acc_oms.data() + static_cast<size_t>(r) * k_cols * k_cols;
    for (int col = 0; col < k_cols; ++col) {
      for (int row = 0; row < k_cols; ++row) {
        const double sym_sum = acc_r[row + static_cast<size_t>(col) * k_cols] +
                               acc_r[col + static_cast<size_t>(row) * k_cols];
        om[row + static_cast<size_t>(col) * k_cols] = gram[row + static_cast<size_t>(col) * k_cols] + scale * sym_sum;
      }
    }
    SET_VECTOR_ELT(oms_sexp, r, om_sexp);
    UNPROTECT(1);
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names_sexp, 0, Rf_mkChar("Oms"));
  SET_STRING_ELT(names_sexp, 1, Rf_mkChar("state"));
  Rf_setAttrib(out_sexp, R_NamesSymbol, names_sexp);
  SET_VECTOR_ELT(out_sexp, 0, oms_sexp);
  SET_VECTOR_ELT(out_sexp, 1, Rf_ScalarReal(static_cast<double>(rt)));

  UNPROTECT(3);
  return out_sexp;
  } catch (const std::exception& error) {
    const size_t length = std::strlen(error.what()) + 1;
    char* message = static_cast<char*>(R_alloc(length, sizeof(char)));
    std::memcpy(message, error.what(), length);
    error_message = message;
  }
  // Leave the C++ exception handler before R performs its non-local jump.
  Rf_error("Native calculation failed: %s", error_message);
  return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"cpp_lvech", (DL_FUNC)&cpp_lvech, 1},
    {"cpp_demeanmat", (DL_FUNC)&cpp_demeanmat, 1},
    {"cpp_sigma_demean", (DL_FUNC)&cpp_sigma_demean, 2},
    {"cpp_getW", (DL_FUNC)&cpp_getW, 3},
    {"cpp_ln_subset_evecs", (DL_FUNC)&cpp_ln_subset_evecs, 3},
    {"cpp_distmat", (DL_FUNC)&cpp_distmat, 2},
    {"cpp_distvec", (DL_FUNC)&cpp_distvec, 3},
    {"cpp_getavc", (DL_FUNC)&cpp_getavc, 2},
    {"cpp_getc0_from_avc", (DL_FUNC)&cpp_getc0_from_avc, 2},
    {"cpp_getc0_from_distmat", (DL_FUNC)&cpp_getc0_from_distmat, 2},
    {"cpp_getrp", (DL_FUNC)&cpp_getrp, 4},
    {"cpp_maxrp", (DL_FUNC)&cpp_maxrp, 5},
    {"cpp_getcv", (DL_FUNC)&cpp_getcv, 5},
    {"cpp_setfinalW", (DL_FUNC)&cpp_setfinalW, 5},
    {"cpp_get_oms", (DL_FUNC)&cpp_get_oms, 5},
    {"cpp_get_oms_batch", (DL_FUNC)&cpp_get_oms_batch, 5},
    {"cpp_raninds", (DL_FUNC)&cpp_raninds, 3},
    {"cpp_jumble_s", (DL_FUNC)&cpp_jumble_s, 3},
    {"cpp_lnset_wc0_subsamples", (DL_FUNC)&cpp_lnset_wc0_subsamples, 7},
    {"cpp_lnset_wc0_wall", (DL_FUNC)&cpp_lnset_wc0_wall, 5},
    {"cpp_lnset_wc0", (DL_FUNC)&cpp_lnset_wc0, 8},
    {"cpp_lnget_oms", (DL_FUNC)&cpp_lnget_oms, 8},
    {NULL, NULL, 0}};

void R_init_scpcR(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

}  // extern "C"
