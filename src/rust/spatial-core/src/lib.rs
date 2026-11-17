//! Shared numerical core for spatial SCPC and SPUR.
//!
//! Matrices at this boundary are contiguous and row-major.

pub mod basis;
pub mod calibration;
mod covariance;
pub mod distance;
pub mod linalg;
pub mod nystrom;
mod quadrature;
pub mod scpc_op;
mod student_t;
pub mod threads;

pub use distance::Metric;
pub use scpc_op::{
    run_approx_scpc, run_exact_scpc, ApproxScpcRequest, ExactScpcOutput, ExactScpcRequest,
};

use std::fmt;

/// Categorized error kind so host adapters can distinguish recoverable states
/// (e.g. tied eigenspaces, degenerate variance, calibration failure) from bugs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorKind {
    InvalidInput,
    UnsupportedCase,
    TiedEigenspace,
    CalibrationFailed,
    DegenerateVariance,
    InternalFailure,
}

/// A numerical input or calculation error, translated by each host adapter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoreError {
    pub kind: ErrorKind,
    pub message: &'static str,
}

impl CoreError {
    pub const fn new(kind: ErrorKind, message: &'static str) -> Self {
        Self { kind, message }
    }

    pub const fn invalid_input(message: &'static str) -> Self {
        Self::new(ErrorKind::InvalidInput, message)
    }

    pub const fn unsupported(message: &'static str) -> Self {
        Self::new(ErrorKind::UnsupportedCase, message)
    }

    pub const fn tied(message: &'static str) -> Self {
        Self::new(ErrorKind::TiedEigenspace, message)
    }

    pub const fn calibration(message: &'static str) -> Self {
        Self::new(ErrorKind::CalibrationFailed, message)
    }

    pub const fn degenerate(message: &'static str) -> Self {
        Self::new(ErrorKind::DegenerateVariance, message)
    }

    pub const fn internal(message: &'static str) -> Self {
        Self::new(ErrorKind::InternalFailure, message)
    }

    pub const fn kind(&self) -> ErrorKind {
        self.kind
    }

    pub const fn message(&self) -> &'static str {
        self.message
    }
}

impl fmt::Display for CoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.message.fmt(f)
    }
}

impl std::error::Error for CoreError {}

pub(crate) fn checked_len(dimensions: &[usize]) -> Result<usize, CoreError> {
    dimensions
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
        .filter(|&size| size <= isize::MAX as usize / size_of::<f64>())
        .ok_or(CoreError::invalid_input("matrix dimensions are too large"))
}
