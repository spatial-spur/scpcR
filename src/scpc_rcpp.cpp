// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include <Rmath.h>
#include <R_ext/RS.h>

extern "C" {
void F77_NAME(dger)(const int* m, const int* n, const double* alpha,
                    const double* x, const int* incx,
                    const double* y, const int* incy,
                    double* a, const int* lda);
void F77_NAME(dsyevr)(const char* jobz, const char* range, const char* uplo,
                      const int* n, double* a, const int* lda,
                      const double* vl, const double* vu,
                      const int* il, const int* iu,
                      const double* abstol, int* m, double* w,
                      double* z, const int* ldz, int* isuppz,
                      double* work, const int* lwork,
                      int* iwork, const int* liwork, int* info);
}

#ifdef __APPLE__
extern "C" {
void vvexp(double* y, const double* x, const int* n);
void vvsin(double* y, const double* x, const int* n);
void vvasin(double* y, const double* x, const int* n);
void vvsqrt(double* y, const double* x, const int* n);
}
#endif

namespace {

inline int resolve_scpc_threads(int requested = -1) {
  if (requested > 0) return requested;
  const char* cran = std::getenv("_R_CHECK_LIMIT_CORES_");
  if (cran != nullptr && *cran != '\0') return 2;
  const char* env = std::getenv("SCPCR_NUM_THREADS");
  if (env != nullptr && *env != '\0') {
    char* end;
    long val = std::strtol(env, &end, 10);
    if (end != env && *end == '\0' && val > 0) return static_cast<int>(val);
  }
  env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  if (env != nullptr && *env != '\0') {
    char* end;
    long val = std::strtol(env, &end, 10);
    if (end != env && *end == '\0' && val > 0) return static_cast<int>(val);
  }
  return -1;
}

template <typename WorkerType>
inline void scpc_parallel_for(std::size_t begin, std::size_t end, WorkerType& worker, std::size_t grainSize = 1, int nthreads = -1) {
  const int nth = resolve_scpc_threads(nthreads);
  if (nth == 1) {
    worker(begin, end);
  } else {
    RcppParallel::parallelFor(begin, end, worker, grainSize, nth);
  }
}

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

struct ExpLowerTriangleWorker : public RcppParallel::Worker {
  const double* distmat;
  double* S;
  int n;
  double neg_c;

  ExpLowerTriangleWorker(const double* distmat_, double* S_, int n_, double neg_c_)
      : distmat(distmat_), S(S_), n(n_), neg_c(neg_c_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t idx = begin; idx < end; ++idx) {
      const int i = static_cast<int>(idx);
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
};

void exp_lower_triangle_inplace_serial(const double* distmat, double* S, int n, double neg_c) {
  ExpLowerTriangleWorker worker(distmat, S, n, neg_c);
  worker(0, static_cast<std::size_t>(n));
}

void exp_lower_triangle_inplace(const double* distmat, double* S, int n, double neg_c, int nthreads = -1) {
  if (n >= 512) {
    ExpLowerTriangleWorker worker(distmat, S, n, neg_c);
    scpc_parallel_for(0, static_cast<std::size_t>(n), worker, 32, nthreads);
  } else {
    exp_lower_triangle_inplace_serial(distmat, S, n, neg_c);
  }
}

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

struct DistMatLowerWorker : public RcppParallel::Worker {
  const double* s;
  int n;
  int d;
  bool latlong;
  const double* lon_rad;
  const double* lat_rad;
  const double* cos_lat;
  double* out;

  DistMatLowerWorker(const double* s_,
                     int n_,
                     int d_,
                     bool latlong_,
                     const double* lon_rad_,
                     const double* lat_rad_,
                     const double* cos_lat_,
                     double* out_)
      : s(s_),
        n(n_),
        d(d_),
        latlong(latlong_),
        lon_rad(lon_rad_),
        lat_rad(lat_rad_),
        cos_lat(cos_lat_),
        out(out_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    if (latlong) {
      const double inv_pi = 1.0 / M_PI;
#ifdef __APPLE__
      std::vector<double> dlat_buf(n);
      std::vector<double> dlon_buf(n);
#endif
      for (std::size_t idx = begin; idx < end; ++idx) {
        const int i = static_cast<int>(idx);
        double* out_col = out + static_cast<size_t>(i) * n;
        out_col[i] = 0.0;
        const int len = n - 1 - i;
        if (len <= 0) continue;
        const double lon_i = lon_rad[i];
        const double lat_i = lat_rad[i];
        const double cos_i = cos_lat[i];
#ifdef __APPLE__
        for (int k = 0; k < len; ++k) {
          const int j = i + 1 + k;
          dlat_buf[k] = 0.5 * (lat_rad[j] - lat_i);
          dlon_buf[k] = 0.5 * (lon_rad[j] - lon_i);
        }
        vvsin(dlat_buf.data(), dlat_buf.data(), &len);
        vvsin(dlon_buf.data(), dlon_buf.data(), &len);
        double* dst = out_col + i + 1;
        for (int k = 0; k < len; ++k) {
          const int j = i + 1 + k;
          double inside = dlat_buf[k] * dlat_buf[k] + cos_i * cos_lat[j] * (dlon_buf[k] * dlon_buf[k]);
          if (inside > 1.0) inside = 1.0;
          if (inside < 0.0) inside = 0.0;
          dst[k] = inside;
        }
        vvsqrt(dst, dst, &len);
        vvasin(dst, dst, &len);
        for (int k = 0; k < len; ++k) {
          dst[k] *= inv_pi;
        }
#else
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
#endif
      }
    } else if (d == 2) {
      const double* s0 = s;
      const double* s1 = s + n;
      for (std::size_t idx = begin; idx < end; ++idx) {
        const int i = static_cast<int>(idx);
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
      for (std::size_t idx = begin; idx < end; ++idx) {
        const int i = static_cast<int>(idx);
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
  }
};

void compute_distmat_raw_serial(const double* s, int n, int d, bool latlong, double* out) {
  std::vector<double> lon_rad;
  std::vector<double> lat_rad;
  std::vector<double> cos_lat;
  if (latlong) {
    const double* lon = s;
    const double* lat = s + n;
    const double deg2rad = M_PI / 180.0;
    lon_rad.resize(n);
    lat_rad.resize(n);
    cos_lat.resize(n);
    for (int i = 0; i < n; ++i) {
      lon_rad[i] = lon[i] * deg2rad;
      lat_rad[i] = lat[i] * deg2rad;
      cos_lat[i] = std::cos(lat_rad[i]);
    }
  }
  DistMatLowerWorker worker(
      s, n, d, latlong,
      latlong ? lon_rad.data() : nullptr,
      latlong ? lat_rad.data() : nullptr,
      latlong ? cos_lat.data() : nullptr,
      out);
  worker(0, static_cast<std::size_t>(n));
  symmetrize_lower_to_upper_blocked(out, n);
}

void compute_distmat_raw(const double* s, int n, int d, bool latlong, double* out, int nthreads = -1) {
  if (n < 384) {
    compute_distmat_raw_serial(s, n, d, latlong, out);
    return;
  }
  std::vector<double> lon_rad;
  std::vector<double> lat_rad;
  std::vector<double> cos_lat;
  if (latlong) {
    const double* lon = s;
    const double* lat = s + n;
    const double deg2rad = M_PI / 180.0;
    lon_rad.resize(n);
    lat_rad.resize(n);
    cos_lat.resize(n);
    for (int i = 0; i < n; ++i) {
      lon_rad[i] = lon[i] * deg2rad;
      lat_rad[i] = lat[i] * deg2rad;
      cos_lat[i] = std::cos(lat_rad[i]);
    }
  }
  DistMatLowerWorker worker(
      s, n, d, latlong,
      latlong ? lon_rad.data() : nullptr,
      latlong ? lat_rad.data() : nullptr,
      latlong ? cos_lat.data() : nullptr,
      out);
  scpc_parallel_for(0, static_cast<std::size_t>(n), worker, 16, nthreads);
  symmetrize_lower_to_upper_blocked(out, n);
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
    double c0 = 10.0;
    double c1 = 10.0;
    while (is_below(c0, avc0)) {
      c1 = c0;
      c0 = 0.5 * c0;
    }
    while (c1 < 5000.0 && is_above(c1, avc0)) {
      c0 = c1;
      c1 = 2.0 * c1;
    }
    double c = std::sqrt(c0 * c1);
    while (true) {
      c = std::sqrt(c0 * c1);
      if (is_above(c, avc0)) {
        c0 = c;
      } else {
        c1 = c;
      }
      if (c1 - c0 < 0.001) break;
    }
    out[t] = c;
  }
}

struct RpWorkspace {
  std::vector<double> omx;
  std::vector<double> wr;
  std::vector<double> wi;
  std::vector<double> work;
  std::vector<double> evals;

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

  char jobvl = 'N';
  char jobvr = 'N';
  arma::blas_int n_lapack = k;
  arma::blas_int lwork = static_cast<arma::blas_int>(ws.work.size());
  arma::blas_int info = 0;
  double dummy = 0.0;
  arma::blas_int ldummy = 1;

  arma::lapack::geev<double>(&jobvl, &jobvr, &n_lapack, omx, &n_lapack,
                             ws.wr.data(), ws.wi.data(),
                             &dummy, &ldummy, &dummy, &ldummy,
                             ws.work.data(), &lwork, &info);

  int n_neg = 0;
  double denom = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < k; ++i) {
    const double v = ws.wr[i];
    if (v > denom) denom = v;
    if (v < 0.0) {
      ws.evals[n_neg++] = -v;
    }
  }

  if (n_neg == 0 || !std::isfinite(denom) || denom <= 0.0) {
    return 0.0;
  }

  const double inv_denom = 1.0 / denom;
  for (int i = 0; i < n_neg; ++i) {
    ws.evals[i] *= inv_denom;
  }

  double tot = 0.0;
  const double* ev = ws.evals.data();
  for (int j = 0; j < n_gq; ++j) {
    const double u = gqx[j];
    const double one_minus_u2 = 1.0 - u * u;
    const double inv_one_minus_u2 = 1.0 / one_minus_u2;
    double prod = 1.0;
    for (int l = 0; l < n_neg; ++l) {
      prod *= (1.0 + ev[l] * inv_one_minus_u2);
    }
    double arg = one_minus_u2 * prod;
    if (!std::isfinite(arg)) {
      double sum_log = 0.0;
      for (int l = 0; l < n_neg; ++l) {
        sum_log += std::log1p(ev[l] * inv_one_minus_u2);
      }
      arg = one_minus_u2 * std::exp(sum_log);
    }
    if (!std::isfinite(arg) || arg <= 0.0) continue;
    tot += gqw[j] / std::sqrt(arg);
  }
  return tot * (2.0 / M_PI);
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
                       double cv0_init,
                       const double* gqx,
                       const double* gqw,
                       int n_gq,
                       RpWorkspace& ws) {
  const int k = q + 1;
  const double sqrt_q = std::sqrt(static_cast<double>(q));
  const double inv_sqrt_q = 1.0 / sqrt_q;
  const double tol = 0.001 * inv_sqrt_q;

  double rp = 1.0;
  int i = 0;
  double cv0 = cv0_init;
  double cv1 = cv0;
  bool first_iter = true;

  while (rp > level) {
    cv1 = cv0;
    double rp_cv1 = rp;
    if (!first_iter) {
      cv0 = cv1;
      cv1 = cv1 + inv_sqrt_q;
    }
    first_iter = false;

    while (true) {
      rp_cv1 = eval_getrp_submatrix(oms_ptrs[i], oms_lds[i], k, cv1, gqx, gqw, n_gq, ws);
      if (rp_cv1 > level) {
        cv0 = cv1;
        cv1 = cv1 + inv_sqrt_q;
      } else {
        break;
      }
    }

    while (cv1 - cv0 > tol) {
      const double cv = 0.5 * (cv0 + cv1);
      const double rp_mid = eval_getrp_submatrix(oms_ptrs[i], oms_lds[i], k, cv, gqx, gqw, n_gq, ws);
      if (rp_mid > level) {
        cv0 = cv;
      } else {
        cv1 = cv;
        rp_cv1 = rp_mid;
      }
    }

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
                   &work_query, &lwork, &iwork_query, &liwork, &info);

  lwork = std::max(26 * n, static_cast<int>(work_query) + 1);
  liwork = std::max(10 * n, iwork_query + 1);
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);

  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, A, &n, &vl, &vu, &il, &iu,
                   &abstol, &m_found, evals.data(), Z.data(), &n, isuppz.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info);

  for (int c = 0; c < k; ++c) {
    const double* src = Z.data() + static_cast<size_t>(k - 1 - c) * n;
    double* dst = V_out + static_cast<size_t>(c) * n;
    std::memcpy(dst, src, static_cast<size_t>(n) * sizeof(double));
  }
}

void compute_top_k_evecs_sym(double* A, int n, int qmax, double* V_out) {
  const int k = std::min(qmax, n);
  const int m_k = std::min(n, std::max(2 * k + 16, 44));
  if (n <= 450 || m_k >= n / 2) {
    compute_top_k_evecs_dsyevr(A, n, k, V_out);
    return;
  }

  std::vector<double> V(static_cast<size_t>(n) * (m_k + 1), 0.0);
  std::vector<double> H(static_cast<size_t>(m_k) * m_k, 0.0);
  std::vector<double> H_copy(static_cast<size_t>(m_k) * m_k, 0.0);
  std::vector<double> w(n, 0.0);
  std::vector<double> h_vec(m_k + 1, 0.0);
  std::vector<double> c_vec(m_k + 1, 0.0);
  std::vector<double> ritz_vals(m_k, 0.0);
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
  const arma::blas_int inc1 = 1;
  const arma::blas_int n_b = n;
  arma::blas_int mk_b = m_k;
  const arma::blas_int k_b = k;

  int k_start = 0;
  double beta_last = 0.0;
  bool converged = false;
  constexpr int MAX_RESTARTS = 25;

  for (int restart = 0; restart < MAX_RESTARTS; ++restart) {
    for (int j = k_start; j < m_k; ++j) {
      const double* v_j = V.data() + static_cast<size_t>(j) * n;
      arma::blas::gemv<double>(&trans_n, &n_b, &n_b, &one, A, &n_b, v_j, &inc1, &zero, w.data(), &inc1);

      const arma::blas_int cols_j = j + 1;
      arma::blas::gemv<double>(&trans_t, &n_b, &cols_j, &one, V.data(), &n_b, w.data(), &inc1, &zero, h_vec.data(), &inc1);
      arma::blas::gemv<double>(&trans_n, &n_b, &cols_j, &neg_one, V.data(), &n_b, h_vec.data(), &inc1, &one, w.data(), &inc1);

      arma::blas::gemv<double>(&trans_t, &n_b, &cols_j, &one, V.data(), &n_b, w.data(), &inc1, &zero, c_vec.data(), &inc1);
      arma::blas::gemv<double>(&trans_n, &n_b, &cols_j, &neg_one, V.data(), &n_b, c_vec.data(), &inc1, &one, w.data(), &inc1);

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
        arma::blas::gemv<double>(&trans_t, &n_b, &cols_j, &one, V.data(), &n_b, v_next, &inc1, &zero, h_vec.data(), &inc1);
        arma::blas::gemv<double>(&trans_n, &n_b, &cols_j, &neg_one, V.data(), &n_b, h_vec.data(), &inc1, &one, v_next, &inc1);
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

    std::memcpy(H_copy.data(), H.data(), static_cast<size_t>(m_k) * m_k * sizeof(double));
    char jobz = 'V';
    char uplo = 'L';
    arma::blas_int lwork = std::max(1, 3 * m_k);
    std::vector<double> work(lwork);
    arma::blas_int info = 0;
    arma::lapack::syev<double>(&jobz, &uplo, &mk_b, H_copy.data(), &mk_b, ritz_vals.data(), work.data(), &lwork, &info);
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
      arma::blas::gemm<double>(&trans_n, &trans_n, &n_b, &k_b, &mk_b, &one,
                               V.data(), &n_b, U_keep.data(), &mk_b, &zero, V_out, &n_b);
      converged = true;
      break;
    }

    const int k_keep = std::min(m_k - 2, k + std::max(4, (m_k - k) / 2));
    const arma::blas_int kk_b = k_keep;
    for (int c = 0; c < k_keep; ++c) {
      const double* src = H_copy.data() + static_cast<size_t>(m_k - 1 - c) * m_k;
      std::memcpy(U_keep.data() + static_cast<size_t>(c) * m_k, src, static_cast<size_t>(m_k) * sizeof(double));
    }
    arma::blas::gemm<double>(&trans_n, &trans_n, &n_b, &kk_b, &mk_b, &one,
                             V.data(), &n_b, U_keep.data(), &mk_b, &zero, V_new.data(), &n_b);
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

struct GetCvWorker : public RcppParallel::Worker {
  const std::vector<const double*>& oms_ptrs;
  const std::vector<int>& oms_lds;
  int q;
  const double* levels;
  const double* cv0_inits;
  const double* gqx;
  const double* gqw;
  int n_gq;
  double* out;

  GetCvWorker(const std::vector<const double*>& oms_ptrs_,
              const std::vector<int>& oms_lds_,
              int q_,
              const double* levels_,
              const double* cv0_inits_,
              const double* gqx_,
              const double* gqw_,
              int n_gq_,
              double* out_)
      : oms_ptrs(oms_ptrs_),
        oms_lds(oms_lds_),
        q(q_),
        levels(levels_),
        cv0_inits(cv0_inits_),
        gqx(gqx_),
        gqw(gqw_),
        n_gq(n_gq_),
        out(out_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    RpWorkspace ws;
    for (std::size_t l = begin; l < end; ++l) {
      out[l] = eval_getcv_core(oms_ptrs, oms_lds, q, levels[l], cv0_inits[l], gqx, gqw, n_gq, ws);
    }
  }
};

struct SetFinalWWorker : public RcppParallel::Worker {
  const std::vector<const double*>& oms_ptrs;
  const std::vector<int>& oms_lds;
  const double* cv0_inits;
  const double* gqx;
  const double* gqw;
  int n_gq;
  double* cvs_out;
  double* lens_out;

  SetFinalWWorker(const std::vector<const double*>& oms_ptrs_,
                  const std::vector<int>& oms_lds_,
                  const double* cv0_inits_,
                  const double* gqx_,
                  const double* gqw_,
                  int n_gq_,
                  double* cvs_out_,
                  double* lens_out_)
      : oms_ptrs(oms_ptrs_),
        oms_lds(oms_lds_),
        cv0_inits(cv0_inits_),
        gqx(gqx_),
        gqw(gqw_),
        n_gq(n_gq_),
        cvs_out(cvs_out_),
        lens_out(lens_out_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    RpWorkspace ws;
    for (std::size_t idx = begin; idx < end; ++idx) {
      const int q = static_cast<int>(idx) + 1;
      const double cv_q = eval_getcv_core(oms_ptrs, oms_lds, q, 0.05, cv0_inits[idx], gqx, gqw, n_gq, ws);
      const double len_q = cv_q * std::tgamma(0.5 * (q + 1)) /
                           (std::sqrt(static_cast<double>(q)) * std::tgamma(0.5 * q));
      cvs_out[idx] = cv_q;
      lens_out[idx] = len_q;
    }
  }
};

struct GetOmsBatchGridWorker : public RcppParallel::Worker {
  const double* distmat;
  const double* W;
  const double* c_vec;
  int n;
  int n_batch;
  int k_cols;
  int total_cols;
  double* SW_all;
  double* oms_all;

  GetOmsBatchGridWorker(const double* distmat_,
                        const double* W_,
                        const double* c_vec_,
                        int n_,
                        int n_batch_,
                        int k_cols_,
                        int total_cols_,
                        double* SW_all_,
                        double* oms_all_)
      : distmat(distmat_),
        W(W_),
        c_vec(c_vec_),
        n(n_),
        n_batch(n_batch_),
        k_cols(k_cols_),
        total_cols(total_cols_),
        SW_all(SW_all_),
        oms_all(oms_all_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    constexpr int BLK = 256;
    std::vector<double> D_blk(BLK * BLK);
    std::vector<double> S_blk(BLK * BLK);

    const char trans_n = 'N';
    const char trans_t = 'T';
    const double alpha = 1.0;
    const double beta_one = 1.0;
    const double beta_zero = 0.0;
    const arma::blas_int n_b = n;
    const arma::blas_int tc_b = total_cols;
    const arma::blas_int kc_b = k_cols;

    for (int j0 = 0; j0 < n; j0 += BLK) {
      const int bj = std::min(BLK, n - j0);
      const arma::blas_int bj_b = bj;
      for (int i0 = j0; i0 < n; i0 += BLK) {
        const int bi = std::min(BLK, n - i0);
        const arma::blas_int bi_b = bi;
        const int len = bi * bj;

        for (int j = 0; j < bj; ++j) {
          const double* d_col = distmat + static_cast<size_t>(j0 + j) * n + i0;
          std::memcpy(D_blk.data() + static_cast<size_t>(j) * bi, d_col, bi * sizeof(double));
        }

        const double* W_i = W + i0;
        const double* W_j = W + j0;
        const double* d_ptr = D_blk.data();
        double* s_ptr = S_blk.data();

        for (std::size_t m = begin; m < end; ++m) {
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

          double* SW_m = SW_all + m * static_cast<size_t>(n) * total_cols;
          if (i0 == j0) {
            arma::blas::gemm<double>(&trans_n, &trans_n, &bi_b, &tc_b, &bj_b, &alpha,
                                     s_ptr, &bi_b, W_i, &n_b, &beta_one, SW_m + i0, &n_b);
          } else {
            arma::blas::gemm<double>(&trans_n, &trans_n, &bi_b, &tc_b, &bj_b, &alpha,
                                     s_ptr, &bi_b, W_j, &n_b, &beta_one, SW_m + i0, &n_b);
            arma::blas::gemm<double>(&trans_t, &trans_n, &bj_b, &tc_b, &bi_b, &alpha,
                                     s_ptr, &bi_b, W_i, &n_b, &beta_one, SW_m + j0, &n_b);
          }
        }
      }
    }

    for (std::size_t m = begin; m < end; ++m) {
      for (int b = 0; b < n_batch; ++b) {
        const double* W_b = W + static_cast<size_t>(b) * n * k_cols;
        const double* SW_bm = SW_all + m * static_cast<size_t>(n) * total_cols + static_cast<size_t>(b) * n * k_cols;
        double* om_bm = oms_all + (m * static_cast<size_t>(n_batch) + static_cast<size_t>(b)) * k_cols * k_cols;
        arma::blas::gemm<double>(&trans_t, &trans_n, &kc_b, &kc_b, &n_b, &alpha,
                                 W_b, &n_b, SW_bm, &n_b, &beta_zero, om_bm, &kc_b);
      }
    }
  }
};

struct LnSubsampleDistWorker : public RcppParallel::Worker {
  const double* sub_coords_all;
  int m;
  int d;
  bool latlong;
  size_t block_len;
  double* sub_dist_all;
  double* distvec;
  double* min_d_vec;
  double* max_d_vec;

  LnSubsampleDistWorker(const double* sub_coords_all_,
                        int m_,
                        int d_,
                        bool latlong_,
                        size_t block_len_,
                        double* sub_dist_all_,
                        double* distvec_,
                        double* min_d_vec_,
                        double* max_d_vec_)
      : sub_coords_all(sub_coords_all_),
        m(m_),
        d(d_),
        latlong(latlong_),
        block_len(block_len_),
        sub_dist_all(sub_dist_all_),
        distvec(distvec_),
        min_d_vec(min_d_vec_),
        max_d_vec(max_d_vec_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t idx = begin; idx < end; ++idx) {
      const double* sc_ptr = sub_coords_all + idx * static_cast<size_t>(m) * d;
      double* sd_ptr = sub_dist_all + idx * static_cast<size_t>(m) * m;
      compute_distmat_raw_serial(sc_ptr, m, d, latlong, sd_ptr);

      double local_min = std::numeric_limits<double>::infinity();
      double local_max = -std::numeric_limits<double>::infinity();
      size_t pos = idx * block_len;
      for (int j = 0; j < m - 1; ++j) {
        const double* col = sd_ptr + static_cast<size_t>(j) * m + (j + 1);
        const size_t len = static_cast<size_t>(m - 1 - j);
        std::memcpy(distvec + pos, col, len * sizeof(double));
        for (size_t k = 0; k < len; ++k) {
          const double v = col[k];
          if (v < local_min) local_min = v;
          if (v > local_max) local_max = v;
        }
        pos += len;
      }
      min_d_vec[idx] = local_min;
      max_d_vec[idx] = local_max;
    }
  }
};

struct LnSubsampleWallWorker : public RcppParallel::Worker {
  const double* S;
  const double* sub_coords_all;
  const double* sub_dist_all;
  int n;
  int d;
  int m;
  int capN;
  int qmax;
  double neg_c0;
  double* Wall;

  LnSubsampleWallWorker(const double* S_,
                        const double* sub_coords_all_,
                        const double* sub_dist_all_,
                        int n_,
                        int d_,
                        int m_,
                        int capN_,
                        int qmax_,
                        double neg_c0_,
                        double* Wall_)
      : S(S_),
        sub_coords_all(sub_coords_all_),
        sub_dist_all(sub_dist_all_),
        n(n_),
        d(d_),
        m(m_),
        capN(capN_),
        qmax(qmax_),
        neg_c0(neg_c0_),
        Wall(Wall_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    std::vector<double> sig_d(static_cast<size_t>(m) * m);
    std::vector<double> W0(static_cast<size_t>(m) * qmax);
    std::vector<double> K_mat(static_cast<size_t>(n) * m);
    std::vector<double> Wx(static_cast<size_t>(n) * qmax);

    const char trans_n = 'N';
    const double alpha = 1.0;
    const double beta = 0.0;
    const double inv_n = 1.0 / static_cast<double>(n);

    for (std::size_t idx_st = begin; idx_st < end; ++idx_st) {
      const int idx = static_cast<int>(idx_st);
      const double* sub_coords = sub_coords_all + static_cast<size_t>(idx) * m * d;
      const double* sub_dist = sub_dist_all + static_cast<size_t>(idx) * m * m;

      exp_lower_triangle_inplace_serial(sub_dist, sig_d.data(), m, neg_c0);
      symmetrize_lower_to_upper_blocked(sig_d.data(), m);
      demean_rect_inplace(sig_d.data(), m, m);
      compute_top_k_evecs_sym(sig_d.data(), m, qmax, W0.data());

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
              const double diff =
                  S[i + static_cast<size_t>(col_d) * n] - sub_coords[j + static_cast<size_t>(col_d) * m];
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

      const arma::blas_int n_b = n;
      const arma::blas_int qmax_b = qmax;
      const arma::blas_int m_b = m;
      arma::blas::gemm<double>(&trans_n, &trans_n, &n_b, &qmax_b, &m_b, &alpha,
                               K_mat.data(), &n_b, W0.data(), &m_b, &beta, Wx.data(), &n_b);

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
  }
};

struct LnGetOmsPairTableWorker : public RcppParallel::Worker {
  const double* distmat;
  const int* pair_count;
  const double* W;
  const double* c_vec;
  int n;
  int k_cols;
  double scale;
  double* oms_all;

  LnGetOmsPairTableWorker(const double* distmat_,
                          const int* pair_count_,
                          const double* W_,
                          const double* c_vec_,
                          int n_,
                          int k_cols_,
                          double scale_,
                          double* oms_all_)
      : distmat(distmat_),
        pair_count(pair_count_),
        W(W_),
        c_vec(c_vec_),
        n(n_),
        k_cols(k_cols_),
        scale(scale_),
        oms_all(oms_all_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    std::vector<double> A_mat(static_cast<size_t>(n) * n, 0.0);
    std::vector<double> AW(static_cast<size_t>(n) * k_cols);
    const char trans_n = 'N';
    const char trans_t = 'T';
    const double alpha = 1.0;
    const double beta_zero = 0.0;

    for (std::size_t r = begin; r < end; ++r) {
      const double neg_c = -c_vec[r];
      exp_lower_triangle_inplace_serial(distmat, A_mat.data(), n, neg_c);
      for (int u = 0; u < n; ++u) {
        A_mat[u + static_cast<size_t>(u) * n] = 0.0;
        double* a_col = A_mat.data() + static_cast<size_t>(u) * n;
        const int* cnt_col = pair_count + static_cast<size_t>(u) * n;
        for (int v = u + 1; v < n; ++v) {
          a_col[v] *= static_cast<double>(cnt_col[v]);
        }
      }
      symmetrize_lower_to_upper_blocked(A_mat.data(), n);

      const arma::blas_int n_b = n;
      const arma::blas_int kc_b = k_cols;
      arma::blas::gemm<double>(&trans_n, &trans_n, &n_b, &kc_b, &n_b, &alpha,
                               A_mat.data(), &n_b, W, &n_b, &beta_zero, AW.data(), &n_b);

      double* om = oms_all + r * static_cast<size_t>(k_cols) * k_cols;
      arma::blas::gemm<double>(&trans_t, &trans_n, &kc_b, &kc_b, &n_b, &scale,
                               W, &n_b, AW.data(), &n_b, &beta_zero, om, &kc_b);
      for (int a = 0; a < k_cols; ++a) {
        om[a + static_cast<size_t>(a) * k_cols] += 1.0;
      }
    }
  }
};

struct LnGetOmsTileDistWorker : public RcppParallel::Worker {
  const double* S;
  const double* lon_rad;
  const double* lat_rad;
  const double* cos_lat;
  const int* idx1_all;
  const int* idx2_all;
  int n;
  int d;
  int capM;
  int tile_size;
  bool latlong;
  double* dist_all;

  LnGetOmsTileDistWorker(const double* S_,
                         const double* lon_rad_,
                         const double* lat_rad_,
                         const double* cos_lat_,
                         const int* idx1_all_,
                         const int* idx2_all_,
                         int n_,
                         int d_,
                         int capM_,
                         int tile_size_,
                         bool latlong_,
                         double* dist_all_)
      : S(S_),
        lon_rad(lon_rad_),
        lat_rad(lat_rad_),
        cos_lat(cos_lat_),
        idx1_all(idx1_all_),
        idx2_all(idx2_all_),
        n(n_),
        d(d_),
        capM(capM_),
        tile_size(tile_size_),
        latlong(latlong_),
        dist_all(dist_all_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    const double inv_pi = 1.0 / M_PI;
#ifdef __APPLE__
    std::vector<double> dlat_tile(tile_size);
    std::vector<double> dlon_tile(tile_size);
#endif
    for (std::size_t tile_idx = begin; tile_idx < end; ++tile_idx) {
      const int m0 = static_cast<int>(tile_idx) * tile_size;
      const int b = std::min(tile_size, capM - m0);
      const int* idx1_tile = idx1_all + m0;
      const int* idx2_tile = idx2_all + m0;
      double* dist_tile = dist_all + m0;

      if (latlong) {
#ifdef __APPLE__
        for (int t = 0; t < b; ++t) {
          const int u = idx1_tile[t];
          const int v = idx2_tile[t];
          dlat_tile[t] = 0.5 * (lat_rad[u] - lat_rad[v]);
          dlon_tile[t] = 0.5 * (lon_rad[u] - lon_rad[v]);
        }
        vvsin(dlat_tile.data(), dlat_tile.data(), &b);
        vvsin(dlon_tile.data(), dlon_tile.data(), &b);
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
    }
  }
};

struct LnGetOmsStreamGridWorker : public RcppParallel::Worker {
  const double* W;
  const int* idx1_all;
  const int* idx2_all;
  const double* dist_all;
  const double* c_vec;
  int n;
  int k_cols;
  int capM;
  int tile_size;
  double* acc_oms;

  LnGetOmsStreamGridWorker(const double* W_,
                           const int* idx1_all_,
                           const int* idx2_all_,
                           const double* dist_all_,
                           const double* c_vec_,
                           int n_,
                           int k_cols_,
                           int capM_,
                           int tile_size_,
                           double* acc_oms_)
      : W(W_),
        idx1_all(idx1_all_),
        idx2_all(idx2_all_),
        dist_all(dist_all_),
        c_vec(c_vec_),
        n(n_),
        k_cols(k_cols_),
        capM(capM_),
        tile_size(tile_size_),
        acc_oms(acc_oms_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    std::vector<double> cd_tile(tile_size);
    std::vector<double> W1_tile(static_cast<size_t>(tile_size) * k_cols);
    std::vector<double> W2_tile(static_cast<size_t>(tile_size) * k_cols);
    std::vector<double> W2_scaled(static_cast<size_t>(tile_size) * k_cols);

    const char trans_n = 'N';
    const char trans_t = 'T';
    const double alpha = 1.0;
    const double beta_one = 1.0;

    for (int m0 = 0; m0 < capM; m0 += tile_size) {
      const int b = std::min(tile_size, capM - m0);
      const int* idx1_tile = idx1_all + m0;
      const int* idx2_tile = idx2_all + m0;
      const double* dist_tile = dist_all + m0;

      for (int col = 0; col < k_cols; ++col) {
        const double* w_col = W + static_cast<size_t>(col) * n;
        double* w1_col = W1_tile.data() + static_cast<size_t>(col) * b;
        double* w2_col = W2_tile.data() + static_cast<size_t>(col) * b;
        for (int t = 0; t < b; ++t) {
          w1_col[t] = w_col[idx1_tile[t]];
          w2_col[t] = w_col[idx2_tile[t]];
        }
      }

      for (std::size_t r = begin; r < end; ++r) {
        const double neg_c = -c_vec[r];
#ifdef __APPLE__
        for (int t = 0; t < b; ++t) {
          cd_tile[t] = neg_c * dist_tile[t];
        }
        vvexp(cd_tile.data(), cd_tile.data(), &b);
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

        double* acc_r = acc_oms + r * static_cast<size_t>(k_cols) * k_cols;
        const arma::blas_int kc_b = k_cols;
        const arma::blas_int b_b = b;
        arma::blas::gemm<double>(&trans_t, &trans_n, &kc_b, &kc_b, &b_b, &alpha,
                                 W1_tile.data(), &b_b, W2_scaled.data(), &b_b, &beta_one, acc_r, &kc_b);
      }
    }
  }
};

}  // namespace

// [[Rcpp::export]]
arma::vec rcpp_lvech(const arma::mat& mat) {
  const int n = static_cast<int>(mat.n_rows);
  const size_t m = static_cast<size_t>(n) * (n - 1) / 2;
  arma::vec out(m);
  const double* src = mat.memptr();
  double* dst = out.memptr();

  size_t pos = 0;
  for (int j = 0; j < n - 1; ++j) {
    const double* col = src + static_cast<size_t>(j) * n + (j + 1);
    const size_t len = static_cast<size_t>(n - 1 - j);
    std::memcpy(dst + pos, col, len * sizeof(double));
    pos += len;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat rcpp_demeanmat(const arma::mat& mat) {
  arma::mat out(mat);
  demean_rect_inplace(out.memptr(), static_cast<int>(out.n_rows), static_cast<int>(out.n_cols));
  return out;
}

// [[Rcpp::export]]
arma::mat rcpp_sigma_demean(const arma::mat& distmat, double c0, int nthreads = -1) {
  const int n = static_cast<int>(distmat.n_rows);
  arma::mat out(n, n);
  exp_lower_triangle_inplace(distmat.memptr(), out.memptr(), n, -c0, nthreads);
  symmetrize_lower_to_upper_blocked(out.memptr(), n);
  demean_rect_inplace(out.memptr(), n, n);
  return out;
}

// [[Rcpp::export]]
arma::mat rcpp_getW(const arma::mat& distmat, double c0, int qmax, int nthreads = -1) {
  const int n = static_cast<int>(distmat.n_rows);
  const int q_use = std::min(qmax, std::max(1, n - 1));

  arma::mat sig_d(n, n);
  exp_lower_triangle_inplace(distmat.memptr(), sig_d.memptr(), n, -c0, nthreads);
  symmetrize_lower_to_upper_blocked(sig_d.memptr(), n);
  demean_rect_inplace(sig_d.memptr(), n, n);

  arma::mat out(n, q_use + 1);
  double* out_ptr = out.memptr();
  const double inv_sqrt_n = 1.0 / std::sqrt(static_cast<double>(n));
  for (int i = 0; i < n; ++i) {
    out_ptr[i] = inv_sqrt_n;
  }

  compute_top_k_evecs_sym(sig_d.memptr(), n, q_use, out_ptr + n);
  return out;
}

// [[Rcpp::export]]
arma::mat rcpp_ln_subset_evecs(const arma::mat& distmat, double c0, int qmax) {
  const int n = static_cast<int>(distmat.n_rows);
  const int q_use = std::min(qmax, n);

  arma::mat sig_d(n, n);
  exp_lower_triangle_inplace_serial(distmat.memptr(), sig_d.memptr(), n, -c0);
  symmetrize_lower_to_upper_blocked(sig_d.memptr(), n);
  demean_rect_inplace(sig_d.memptr(), n, n);

  arma::mat out(n, q_use);
  compute_top_k_evecs_sym(sig_d.memptr(), n, q_use, out.memptr());
  return out;
}

// [[Rcpp::export]]
arma::mat rcpp_distmat(const arma::mat& s, bool latlong, int nthreads = -1) {
  const int n = static_cast<int>(s.n_rows);
  const int d = static_cast<int>(s.n_cols);
  if (latlong && d != 2) {
    Rcpp::stop("Internal error: geodesic coordinates must have exactly two columns.");
  }
  arma::mat out(n, n);
  compute_distmat_raw(s.memptr(), n, d, latlong, out.memptr(), nthreads);
  return out;
}

// [[Rcpp::export]]
arma::vec rcpp_distvec(const arma::mat& s1, const arma::mat& s2, bool latlong) {
  const int n = static_cast<int>(s1.n_rows);
  const int d = static_cast<int>(s1.n_cols);
  const double* p1 = s1.memptr();
  const double* p2 = s2.memptr();

  arma::vec out(n);
  double* out_ptr = out.memptr();

  if (latlong) {
    const double* lon1 = p1;
    const double* lat1 = p1 + n;
    const double* lon2 = p2;
    const double* lat2 = p2 + n;
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
      out_ptr[i] = std::asin(std::sqrt(inside)) * inv_pi;
    }
  } else if (d == 2) {
    const double* x1 = p1;
    const double* y1 = p1 + n;
    const double* x2 = p2;
    const double* y2 = p2 + n;
    for (int i = 0; i < n; ++i) {
      const double dx = x1[i] - x2[i];
      const double dy = y1[i] - y2[i];
      out_ptr[i] = std::sqrt(dx * dx + dy * dy);
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double sum_sq = 0.0;
      for (int k = 0; k < d; ++k) {
        const double diff = p1[i + static_cast<size_t>(k) * n] - p2[i + static_cast<size_t>(k) * n];
        sum_sq += diff * diff;
      }
      out_ptr[i] = std::sqrt(sum_sq);
    }
  }

  return out;
}

// [[Rcpp::export]]
double rcpp_getavc(double c, const arma::vec& dist) {
  const double neg_c = -c;
  const size_t m = dist.n_elem;
  const double* vd = dist.memptr();
  if (m == 0) {
    return NA_REAL;
  }
  double sum = 0.0;
#ifdef __APPLE__
  constexpr int CHUNK = 2048;
  double arg_buf[CHUNK];
  double exp_buf[CHUNK];
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
  return sum / static_cast<double>(m);
}

// [[Rcpp::export]]
arma::vec rcpp_getc0_from_avc(const arma::vec& dist, const arma::vec& avc_vec) {
  const size_t m = dist.n_elem;
  const double* vd = dist.memptr();
  const int n_targets = static_cast<int>(avc_vec.n_elem);

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

  arma::vec out(n_targets);
  solve_c0_from_avc_core(vd, m, bins, history, avc_vec.memptr(), n_targets, out.memptr());
  return out;
}

// [[Rcpp::export]]
arma::vec rcpp_getc0_from_distmat(const arma::mat& distmat, const arma::vec& avc_vec) {
  const int n = static_cast<int>(distmat.n_rows);
  const double* dm = distmat.memptr();
  const int n_targets = static_cast<int>(avc_vec.n_elem);

  const size_t m = static_cast<size_t>(n) * (n - 1) / 2;
  std::vector<double> vd(m);
  double min_d = std::numeric_limits<double>::infinity();
  double max_d = -std::numeric_limits<double>::infinity();

  size_t pos = 0;
  for (int j = 0; j < n - 1; ++j) {
    const double* col = dm + static_cast<size_t>(j) * n + (j + 1);
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

  arma::vec out(n_targets);
  solve_c0_from_avc_core(vd.data(), m, bins, history, avc_vec.memptr(), n_targets, out.memptr());
  return out;
}

// [[Rcpp::export]]
double rcpp_getrp(const arma::mat& om, double cv, const arma::vec& gqx, const arma::vec& gqw) {
  const int k = static_cast<int>(om.n_rows);
  RpWorkspace ws;
  return eval_getrp_submatrix(om.memptr(), k, k, cv, gqx.memptr(), gqw.memptr(), static_cast<int>(gqx.n_elem), ws);
}

// [[Rcpp::export]]
Rcpp::List rcpp_maxrp(const Rcpp::List& oms, int q, double cv, const arma::vec& gqx, const arma::vec& gqw) {
  const int nc = oms.size();
  std::vector<const double*> oms_ptrs(nc);
  std::vector<int> oms_lds(nc);
  for (int m = 0; m < nc; ++m) {
    SEXP mat = oms[m];
    oms_ptrs[m] = REAL(mat);
    oms_lds[m] = Rf_nrows(mat);
  }

  RpWorkspace ws;
  double max_rp = 0.0;
  int max_idx = 0;
  eval_maxrp_core(oms_ptrs, oms_lds, q, cv, gqx.memptr(), gqw.memptr(), static_cast<int>(gqx.n_elem),
                  ws, -1, 0.0, max_rp, max_idx);

  return Rcpp::List::create(
      Rcpp::Named("max") = max_rp,
      Rcpp::Named("i") = max_idx + 1);
}

// [[Rcpp::export]]
arma::vec rcpp_getcv(const Rcpp::List& oms, int q, const arma::vec& levels, const arma::vec& gqx, const arma::vec& gqw, int nthreads = -1) {
  const int nc = oms.size();
  const int n_lev = static_cast<int>(levels.n_elem);

  std::vector<const double*> oms_ptrs(nc);
  std::vector<int> oms_lds(nc);
  for (int m = 0; m < nc; ++m) {
    SEXP mat = oms[m];
    oms_ptrs[m] = REAL(mat);
    oms_lds[m] = Rf_nrows(mat);
  }

  const double inv_sqrt_q = 1.0 / std::sqrt(static_cast<double>(q));
  std::vector<double> cv0_inits(n_lev);
  for (int l = 0; l < n_lev; ++l) {
    cv0_inits[l] = Rf_qt(1.0 - levels[l] * 0.5, static_cast<double>(q), 1, 0) * inv_sqrt_q;
  }

  arma::vec out(n_lev);
  GetCvWorker worker(oms_ptrs, oms_lds, q, levels.memptr(), cv0_inits.data(),
                     gqx.memptr(), gqw.memptr(), static_cast<int>(gqx.n_elem), out.memptr());
  if (n_lev > 1) {
    scpc_parallel_for(0, static_cast<std::size_t>(n_lev), worker, 1, nthreads);
  } else {
    worker(0, static_cast<std::size_t>(n_lev));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List rcpp_setfinalW(const Rcpp::List& oms,
                          const arma::mat& W,
                          int qmax,
                          const arma::vec& gqx,
                          const arma::vec& gqw,
                          int nthreads = -1) {
  const int nc = oms.size();
  std::vector<const double*> oms_ptrs(nc);
  std::vector<int> oms_lds(nc);
  for (int m = 0; m < nc; ++m) {
    SEXP mat = oms[m];
    oms_ptrs[m] = REAL(mat);
    oms_lds[m] = Rf_nrows(mat);
  }

  std::vector<double> cv0_inits(qmax);
  for (int q = 1; q <= qmax; ++q) {
    cv0_inits[q - 1] = Rf_qt(1.0 - 0.05 * 0.5, static_cast<double>(q), 1, 0) / std::sqrt(static_cast<double>(q));
  }

  std::vector<double> cvs_out(qmax, 0.0);
  std::vector<double> lens_out(qmax, 0.0);
  SetFinalWWorker worker(oms_ptrs, oms_lds, cv0_inits.data(),
                         gqx.memptr(), gqw.memptr(), static_cast<int>(gqx.n_elem),
                         cvs_out.data(), lens_out.data());
  if (qmax >= 4) {
    scpc_parallel_for(0, static_cast<std::size_t>(qmax), worker, 1, nthreads);
  } else {
    worker(0, static_cast<std::size_t>(qmax));
  }

  double best_len = std::numeric_limits<double>::infinity();
  double best_cv = 0.0;
  int q_opt = 1;
  for (int q = 1; q <= qmax; ++q) {
    const double len_q = lens_out[q - 1];
    if (q == 1 || len_q < best_len) {
      best_len = len_q;
      best_cv = cvs_out[q - 1];
      q_opt = q;
    }
  }

  arma::mat W_out = W.cols(0, q_opt);
  return Rcpp::List::create(
      Rcpp::Named("W") = W_out,
      Rcpp::Named("cv") = best_cv,
      Rcpp::Named("q") = q_opt);
}

// [[Rcpp::export]]
Rcpp::List rcpp_get_oms(const arma::mat& distmat,
                        double c0,
                        double cmax,
                        const arma::mat& W,
                        double cgridfac,
                        int nthreads = -1) {
  const int n = static_cast<int>(distmat.n_rows);
  const int k_cols = static_cast<int>(W.n_cols);
  const int nc = std::max(2, static_cast<int>(std::ceil(std::log(cmax / c0) / std::log(cgridfac))));

  std::vector<double> c_vec(nc);
  double c_cur = c0;
  for (int m = 0; m < nc; ++m) {
    c_vec[m] = c_cur;
    c_cur *= cgridfac;
  }

  std::vector<double> SW_all(static_cast<size_t>(n) * k_cols * nc, 0.0);
  std::vector<double> oms_all(static_cast<size_t>(k_cols) * k_cols * nc, 0.0);

  GetOmsBatchGridWorker worker(
      distmat.memptr(), W.memptr(), c_vec.data(),
      n, 1, k_cols, k_cols, SW_all.data(), oms_all.data());
  if (n >= 150) {
    scpc_parallel_for(0, static_cast<std::size_t>(nc), worker, std::max(1, nc / 8), nthreads);
  } else {
    worker(0, static_cast<std::size_t>(nc));
  }

  Rcpp::List out(nc);
  for (int m = 0; m < nc; ++m) {
    arma::mat om(k_cols, k_cols);
    std::memcpy(om.memptr(),
                oms_all.data() + static_cast<size_t>(m) * k_cols * k_cols,
                static_cast<size_t>(k_cols) * k_cols * sizeof(double));
    out[m] = om;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List rcpp_get_oms_batch(const arma::mat& distmat,
                              double c0,
                              double cmax,
                              const Rcpp::List& w_list,
                              double cgridfac,
                              int nthreads = -1) {
  const int n = static_cast<int>(distmat.n_rows);
  const int n_batch = w_list.size();
  if (n_batch == 0) {
    return Rcpp::List(0);
  }
  if (n_batch == 1) {
    arma::mat W0 = Rcpp::as<arma::mat>(w_list[0]);
    Rcpp::List out(1);
    out[0] = rcpp_get_oms(distmat, c0, cmax, W0, cgridfac, nthreads);
    return out;
  }

  SEXP w0_sexp = w_list[0];
  const int k_cols = Rf_ncols(w0_sexp);
  const int total_cols = n_batch * k_cols;

  std::vector<double> W_stack(static_cast<size_t>(n) * total_cols);
  for (int b = 0; b < n_batch; ++b) {
    SEXP wb_sexp = w_list[b];
    std::memcpy(W_stack.data() + static_cast<size_t>(b) * n * k_cols,
                REAL(wb_sexp),
                static_cast<size_t>(n) * k_cols * sizeof(double));
  }

  const int nc = std::max(2, static_cast<int>(std::ceil(std::log(cmax / c0) / std::log(cgridfac))));
  std::vector<double> c_vec(nc);
  double c_cur = c0;
  for (int m = 0; m < nc; ++m) {
    c_vec[m] = c_cur;
    c_cur *= cgridfac;
  }

  std::vector<double> SW_all(static_cast<size_t>(n) * total_cols * nc, 0.0);
  std::vector<double> oms_all(static_cast<size_t>(k_cols) * k_cols * n_batch * nc, 0.0);

  GetOmsBatchGridWorker worker(
      distmat.memptr(), W_stack.data(), c_vec.data(),
      n, n_batch, k_cols, total_cols, SW_all.data(), oms_all.data());
  if (n >= 150) {
    scpc_parallel_for(0, static_cast<std::size_t>(nc), worker, std::max(1, nc / 8), nthreads);
  } else {
    worker(0, static_cast<std::size_t>(nc));
  }

  Rcpp::List batch_out(n_batch);
  for (int b = 0; b < n_batch; ++b) {
    Rcpp::List oms_b(nc);
    for (int m = 0; m < nc; ++m) {
      arma::mat om(k_cols, k_cols);
      const double* src = oms_all.data() +
                          (static_cast<size_t>(m) * n_batch + static_cast<size_t>(b)) * k_cols * k_cols;
      std::memcpy(om.memptr(), src, static_cast<size_t>(k_cols) * k_cols * sizeof(double));
      oms_b[m] = om;
    }
    batch_out[b] = oms_b;
  }
  return batch_out;
}

// [[Rcpp::export]]
Rcpp::List rcpp_raninds(int n, int capM, double random_t) {
  uint64_t rt = static_cast<uint64_t>(random_t);
  Rcpp::IntegerVector v(capM + 1);
  int* v_ptr = v.begin();

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
  int j = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));

  for (int i = 0; i <= capM; ++i) {
    v_ptr[i] = j + 1;
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
    const int step = static_cast<int>(std::floor(static_cast<double>(n - 1) * (static_cast<double>(rt) * INV_2POW32)));
    j = (j + 1 + step) % n;
  }

  return Rcpp::List::create(
      Rcpp::Named("indices") = v,
      Rcpp::Named("state") = static_cast<double>(rt));
}

// [[Rcpp::export]]
Rcpp::List rcpp_jumble_s(const arma::mat& s, int m, double random_t) {
  const int n = static_cast<int>(s.n_rows);
  const int d = static_cast<int>(s.n_cols);
  uint64_t rt = static_cast<uint64_t>(random_t);

  arma::mat r(s);
  double* r_ptr = r.memptr();

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  for (int i = 0; i < m; ++i) {
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
    const int j = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));
    if (j != i) {
      for (int k = 0; k < d; ++k) {
        std::swap(r_ptr[i + static_cast<size_t>(k) * n], r_ptr[j + static_cast<size_t>(k) * n]);
      }
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("coords") = r,
      Rcpp::Named("state") = static_cast<double>(rt));
}

// [[Rcpp::export]]
Rcpp::List rcpp_lnset_wc0(const arma::mat& S,
                          double avc0,
                          int qmax,
                          double minavc,
                          bool latlong,
                          int capN,
                          int m_in,
                          double random_t,
                          int nthreads = -1) {
  const int n = static_cast<int>(S.n_rows);
  const int d = static_cast<int>(S.n_cols);
  const int m = std::min(m_in, n);
  uint64_t rt = static_cast<uint64_t>(random_t);

  std::vector<double> r(static_cast<size_t>(n) * d);
  std::memcpy(r.data(), S.memptr(), static_cast<size_t>(n) * d * sizeof(double));

  std::vector<double> sub_coords_all(static_cast<size_t>(capN) * m * d);
  constexpr double INV_2POW32 = 1.0 / 4294967296.0;

  // Step 1: Generate the capN coordinate subsamples sequentially to preserve exact LCG state
  for (int idx = 0; idx < capN; ++idx) {
    for (int i = 0; i < m; ++i) {
      rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
      const int j = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));
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
  }

  // Step 2: Parallelize the capN m x m distance matrices across RcppParallel workers
  const size_t block_len = static_cast<size_t>(m) * (m - 1) / 2;
  const size_t total_len = static_cast<size_t>(capN) * block_len;
  std::vector<double> sub_dist_all(static_cast<size_t>(capN) * m * m);
  std::vector<double> distvec(total_len);
  std::vector<double> min_d_vec(capN, std::numeric_limits<double>::infinity());
  std::vector<double> max_d_vec(capN, -std::numeric_limits<double>::infinity());

  LnSubsampleDistWorker dist_worker(
      sub_coords_all.data(), m, d, latlong, block_len,
      sub_dist_all.data(), distvec.data(), min_d_vec.data(), max_d_vec.data());
  if (m >= 100 && capN > 1) {
    scpc_parallel_for(0, static_cast<std::size_t>(capN), dist_worker, 1, nthreads);
  } else {
    dist_worker(0, static_cast<std::size_t>(capN));
  }

  double min_d = std::numeric_limits<double>::infinity();
  double max_d = -std::numeric_limits<double>::infinity();
  for (int idx = 0; idx < capN; ++idx) {
    if (min_d_vec[idx] < min_d) min_d = min_d_vec[idx];
    if (max_d_vec[idx] > max_d) max_d = max_d_vec[idx];
  }

  std::vector<HistBin> bins;
  build_hist_bins(distvec.data(), total_len, min_d, max_d, bins);
  std::vector<EvalPoint> history;
  history.reserve(64);

  const double targets[2] = {avc0, minavc};
  double c_out[2] = {0.0, 0.0};
  solve_c0_from_avc_core(distvec.data(), total_len, bins, history, targets, 2, c_out);
  const double c0 = c_out[0];
  const double cmax = c_out[1];
  const double neg_c0 = -c0;

  // Step 3: Parallelize the capN subsample eigendecompositions + n x m kernel projections
  const int total_wall_cols = capN * qmax;
  std::vector<double> Wall(static_cast<size_t>(n) * total_wall_cols);
  LnSubsampleWallWorker wall_worker(
      S.memptr(), sub_coords_all.data(), sub_dist_all.data(),
      n, d, m, capN, qmax, neg_c0, Wall.data());
  if (capN > 1) {
    scpc_parallel_for(0, static_cast<std::size_t>(capN), wall_worker, 1, nthreads);
  } else {
    wall_worker(0, static_cast<std::size_t>(capN));
  }

  // Step 4: Sequential Gram-matrix deflation across q = 1 .. qmax
  arma::mat W_out(n, qmax + 1);
  double* W_ptr = W_out.memptr();
  const double inv_sqrt_n = 1.0 / std::sqrt(static_cast<double>(n));
  for (int i = 0; i < n; ++i) {
    W_ptr[i] = inv_sqrt_n;
  }

  std::vector<double> G_mat(static_cast<size_t>(total_wall_cols) * total_wall_cols);
  std::vector<double> evec(total_wall_cols);
  std::vector<double> proj(total_wall_cols);
  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta = 0.0;
  const double neg_one = -1.0;
  const int inc1 = 1;

  const arma::blas_int n_b = n;
  const arma::blas_int twc_b = total_wall_cols;
  const arma::blas_int inc1_b = 1;

  for (int i = 1; i <= qmax; ++i) {
    const int ncol_wx = capN * i;
    const arma::blas_int nw_b = ncol_wx;
    arma::blas::gemm<double>(&trans_t, &trans_n, &nw_b, &nw_b, &n_b, &alpha,
                             Wall.data(), &n_b, Wall.data(), &n_b, &beta, G_mat.data(), &nw_b);
    compute_top_k_evecs_sym(G_mat.data(), ncol_wx, 1, evec.data());

    double* w_col = W_ptr + static_cast<size_t>(i) * n;
    arma::blas::gemv<double>(&trans_n, &n_b, &nw_b, &alpha,
                             Wall.data(), &n_b, evec.data(), &inc1_b, &beta, w_col, &inc1_b);

    double norm_sq = 0.0;
    for (int r_idx = 0; r_idx < n; ++r_idx) {
      norm_sq += w_col[r_idx] * w_col[r_idx];
    }
    const double inv_norm = 1.0 / std::sqrt(norm_sq);
    for (int r_idx = 0; r_idx < n; ++r_idx) {
      w_col[r_idx] *= inv_norm;
    }

    if (i < qmax) {
      arma::blas::gemv<double>(&trans_t, &n_b, &twc_b, &alpha,
                               Wall.data(), &n_b, w_col, &inc1_b, &beta, proj.data(), &inc1_b);
      F77_CALL(dger)(&n, &total_wall_cols, &neg_one,
                     w_col, &inc1, proj.data(), &inc1, Wall.data(), &n);
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("W") = W_out,
      Rcpp::Named("c0") = c0,
      Rcpp::Named("cmax") = cmax,
      Rcpp::Named("random_t") = static_cast<double>(rt));
}

// [[Rcpp::export]]
Rcpp::List rcpp_lnget_oms(const arma::mat& S,
                          double c0,
                          double cmax,
                          const arma::mat& W,
                          double cgridfac,
                          int capM,
                          double random_t,
                          bool latlong,
                          int nthreads = -1) {
  const int n = static_cast<int>(S.n_rows);
  const int d = static_cast<int>(S.n_cols);
  const int k_cols = static_cast<int>(W.n_cols);
  uint64_t rt = static_cast<uint64_t>(random_t);

  const int nc = std::max(2, static_cast<int>(std::ceil(std::log(cmax / c0) / std::log(cgridfac))));
  std::vector<double> c_vec(nc);
  double c_cur = c0;
  for (int m = 0; m < nc; ++m) {
    c_vec[m] = c_cur;
    c_cur *= cgridfac;
  }

  constexpr double INV_2POW32 = 1.0 / 4294967296.0;
  const double scale = 0.5 * (static_cast<double>(n) * (n - 1) / static_cast<double>(capM));

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
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;

    std::vector<double> distmat(static_cast<size_t>(n) * n);
    compute_distmat_raw(S.memptr(), n, d, latlong, distmat.data());

    std::vector<double> oms_all(static_cast<size_t>(k_cols) * k_cols * nc, 0.0);
    LnGetOmsPairTableWorker worker(
        distmat.data(), pair_count.data(), W.memptr(), c_vec.data(),
        n, k_cols, scale, oms_all.data());
    if (n >= 150) {
      scpc_parallel_for(0, static_cast<std::size_t>(nc), worker, std::max(1, nc / 8), nthreads);
    } else {
      worker(0, static_cast<std::size_t>(nc));
    }

    Rcpp::List oms_list(nc);
    for (int r = 0; r < nc; ++r) {
      arma::mat om(k_cols, k_cols);
      std::memcpy(om.memptr(),
                  oms_all.data() + static_cast<size_t>(r) * k_cols * k_cols,
                  static_cast<size_t>(k_cols) * k_cols * sizeof(double));
      oms_list[r] = om;
    }

    return Rcpp::List::create(
        Rcpp::Named("Oms") = oms_list,
        Rcpp::Named("state") = static_cast<double>(rt));
  }

  // Streaming tile path for n > 1200 with RcppParallel across tiles & grid points
  std::vector<double> lon_rad;
  std::vector<double> lat_rad;
  std::vector<double> cos_lat;
  if (latlong) {
    const double deg2rad = M_PI / 180.0;
    lon_rad.resize(n);
    lat_rad.resize(n);
    cos_lat.resize(n);
    const double* lon = S.memptr();
    const double* lat = S.memptr() + n;
    for (int i = 0; i < n; ++i) {
      lon_rad[i] = lon[i] * deg2rad;
      lat_rad[i] = lat[i] * deg2rad;
      cos_lat[i] = std::cos(lat_rad[i]);
    }
  }

  // Precompute capM index pairs sequentially from LCG (exact state preservation)
  std::vector<int> idx1_all(capM);
  std::vector<int> idx2_all(capM);
  rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
  int j_prev = static_cast<int>(std::floor(static_cast<double>(n) * (static_cast<double>(rt) * INV_2POW32)));
  for (int m = 0; m < capM; ++m) {
    idx1_all[m] = j_prev;
    rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;
    const int step = static_cast<int>(std::floor(static_cast<double>(n - 1) * (static_cast<double>(rt) * INV_2POW32)));
    j_prev = (j_prev + 1 + step) % n;
    idx2_all[m] = j_prev;
  }
  rt = (64389ULL * rt + 1ULL) & 0xFFFFFFFFULL;

  constexpr int TILE = 4096;
  const int n_tiles = (capM + TILE - 1) / TILE;
  std::vector<double> dist_all(capM);

  LnGetOmsTileDistWorker dist_worker(
      S.memptr(),
      latlong ? lon_rad.data() : nullptr,
      latlong ? lat_rad.data() : nullptr,
      latlong ? cos_lat.data() : nullptr,
      idx1_all.data(), idx2_all.data(),
      n, d, capM, TILE, latlong, dist_all.data());
  scpc_parallel_for(0, static_cast<std::size_t>(n_tiles), dist_worker, 8, nthreads);

  std::vector<double> acc_oms(static_cast<size_t>(k_cols) * k_cols * nc, 0.0);
  LnGetOmsStreamGridWorker grid_worker(
      W.memptr(), idx1_all.data(), idx2_all.data(), dist_all.data(), c_vec.data(),
      n, k_cols, capM, TILE, acc_oms.data());
  scpc_parallel_for(0, static_cast<std::size_t>(nc), grid_worker, std::max(1, nc / 8), nthreads);

  Rcpp::List oms_list(nc);
  for (int r = 0; r < nc; ++r) {
    arma::mat om(k_cols, k_cols);
    double* om_ptr = om.memptr();
    const double* acc_r = acc_oms.data() + static_cast<size_t>(r) * k_cols * k_cols;
    for (int col = 0; col < k_cols; ++col) {
      for (int row = 0; row < k_cols; ++row) {
        const double sym_sum = acc_r[row + static_cast<size_t>(col) * k_cols] +
                               acc_r[col + static_cast<size_t>(row) * k_cols];
        om_ptr[row + static_cast<size_t>(col) * k_cols] = (row == col ? 1.0 : 0.0) + scale * sym_sum;
      }
    }
    oms_list[r] = om;
  }

  return Rcpp::List::create(
      Rcpp::Named("Oms") = oms_list,
      Rcpp::Named("state") = static_cast<double>(rt));
}
