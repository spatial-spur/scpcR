use extendr_api::prelude::*;

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
    spatial_r::run_exact_scpc_r(
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

extendr_module! {
    mod scpcR;
    fn rust_run_exact_scpc;
}
