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

#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_exact_scpc_clustered(
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
    cluster_id: Nullable<Integers>,
) -> extendr_api::Result<List> {
    spatial_r::run_exact_scpc_r_impl(
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
        cluster_id,
        Nullable::Null,
    )
}

#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_exact_scpc_full(
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
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
) -> extendr_api::Result<List> {
    spatial_r::run_exact_scpc_r_impl(
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
        cluster_id,
        fe_id,
    )
}

#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_approx_scpc(
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
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
    seed: f64,
    landmarks: i32,
    samples: i32,
    pair_count: i32,
) -> std::result::Result<List, spatial_r::ApproxScpcError> {
    spatial_r::run_approx_scpc_r_impl(
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
        cluster_id,
        fe_id,
        seed,
        landmarks,
        samples,
        pair_count,
    )
}

extendr_module! {
    mod scpcR;
    fn rust_run_exact_scpc;
    fn rust_run_exact_scpc_clustered;
    fn rust_run_exact_scpc_full;
    fn rust_run_approx_scpc;
}
