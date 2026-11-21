use extendr_api::prelude::*;

#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn rust_run_model_scpc(
    coords: RMatrix<f64>,
    latlong: bool,
    avc: f64,
    cvs: bool,
    coef: Doubles,
    scores: RMatrix<f64>,
    bread: RMatrix<f64>,
    design: Nullable<RMatrix<f64>>,
    instruments: Nullable<RMatrix<f64>>,
    direction_design: Nullable<RMatrix<f64>>,
    cluster_id: Nullable<Integers>,
    fe_id: Nullable<Integers>,
    nthreads: Nullable<i32>,
    seed: Nullable<f64>,
) -> std::result::Result<List, spatial_r::ScpcError> {
    spatial_r::rust_run_model_scpc(
        coords,
        latlong,
        avc,
        cvs,
        coef,
        scores,
        bread,
        design,
        instruments,
        direction_design,
        cluster_id,
        fe_id,
        nthreads,
        seed,
    )
}

#[extendr]
pub fn rust_build_info() -> List {
    spatial_r::build_info_r()
}

extendr_module! {
    mod scpcR;
    fn rust_run_model_scpc;
    fn rust_build_info;
}
