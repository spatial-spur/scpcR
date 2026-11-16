#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifdef SPATIAL_USE_RUST
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

static const R_CallMethodDef CallEntries[] = {
    {"wrap__rust_run_exact_scpc", (DL_FUNC) &wrap__rust_run_exact_scpc, 13},
    {NULL, NULL, 0}
};
#else
static const R_CallMethodDef CallEntries[] = {
    {NULL, NULL, 0}
};
#endif

void R_init_scpcR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
