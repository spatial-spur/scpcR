#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifdef SPATIAL_USE_RUST
extern SEXP wrap__rust_build_info(void);
extern SEXP wrap__rust_run_model_scpc(
    SEXP coords, SEXP latlong, SEXP avc, SEXP cvs, SEXP coef, SEXP scores,
    SEXP bread, SEXP design, SEXP instruments, SEXP direction_design,
    SEXP cluster_id, SEXP fe_id, SEXP nthreads, SEXP seed
);

extern SEXP wrap__rust_run_exact_scpc(
    SEXP coords,
    SEXP latlong,
    SEXP avc,
    SEXP uncond,
    SEXP cvs,
    SEXP coef,
    SEXP wj,
    SEXP xj,
    SEXP xjs,
    SEXP proj_u,
    SEXP proj_v,
    SEXP cv_tolerance,
    SEXP nthreads
);

extern SEXP wrap__rust_run_exact_scpc_clustered(
    SEXP coords,
    SEXP latlong,
    SEXP avc,
    SEXP uncond,
    SEXP cvs,
    SEXP coef,
    SEXP wj,
    SEXP xj,
    SEXP xjs,
    SEXP proj_u,
    SEXP proj_v,
    SEXP cv_tolerance,
    SEXP nthreads,
    SEXP cluster_id
);

extern SEXP wrap__rust_run_exact_scpc_full(
    SEXP coords,
    SEXP latlong,
    SEXP avc,
    SEXP uncond,
    SEXP cvs,
    SEXP coef,
    SEXP wj,
    SEXP xj,
    SEXP xjs,
    SEXP proj_u,
    SEXP proj_v,
    SEXP cv_tolerance,
    SEXP nthreads,
    SEXP cluster_id,
    SEXP fe_id
);

extern SEXP wrap__rust_run_approx_scpc(
    SEXP coords,
    SEXP latlong,
    SEXP avc,
    SEXP uncond,
    SEXP cvs,
    SEXP coef,
    SEXP wj,
    SEXP xj,
    SEXP xjs,
    SEXP proj_u,
    SEXP proj_v,
    SEXP cv_tolerance,
    SEXP nthreads,
    SEXP cluster_id,
    SEXP fe_id,
    SEXP seed,
    SEXP landmarks,
    SEXP samples,
    SEXP pair_count
);

static const R_CallMethodDef CallEntries[] = {
    {"wrap__rust_run_model_scpc", (DL_FUNC) &wrap__rust_run_model_scpc, 14},
    {"wrap__rust_build_info", (DL_FUNC) &wrap__rust_build_info, 0},
    {"wrap__rust_run_exact_scpc", (DL_FUNC) &wrap__rust_run_exact_scpc, 13},
    {"wrap__rust_run_exact_scpc_clustered", (DL_FUNC) &wrap__rust_run_exact_scpc_clustered, 14},
    {"wrap__rust_run_exact_scpc_full", (DL_FUNC) &wrap__rust_run_exact_scpc_full, 15},
    {"wrap__rust_run_approx_scpc", (DL_FUNC) &wrap__rust_run_approx_scpc, 19},
    {NULL, NULL, 0}
};
#else
static const R_CallMethodDef CallEntries[] = {
    {NULL, NULL, 0}
};
#endif

/* Some macOS R builds leave attribute_visible empty even under Clang. */
#if defined(__GNUC__) && !defined(_WIN32) && !defined(HAVE_VISIBILITY_ATTRIBUTE)
__attribute__((visibility("default")))
#endif
void attribute_visible R_init_scpcR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
