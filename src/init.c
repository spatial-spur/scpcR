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

static const R_CallMethodDef CallEntries[] = {
    {"wrap__rust_run_model_scpc", (DL_FUNC) &wrap__rust_run_model_scpc, 14},
    {"wrap__rust_build_info", (DL_FUNC) &wrap__rust_build_info, 0},
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
