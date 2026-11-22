% Please edit documentation in R/core.R
\name{.scpc_accelerate_state}
\alias{.scpc_accelerate_state}
\title{Internal Accelerate fork-state helper}
\usage{
.scpc_accelerate_state(mark = FALSE)
}
\arguments{
\item{mark}{Logical; if \code{TRUE}, record the current process ID as having
initialized Apple Accelerate.}
}
\value{
Logical scalar; \code{TRUE} if Apple Accelerate was initialized in a
parent process before \code{fork()}.
}
\description{
Tracks whether Apple Accelerate was initialized in a parent process across companion packages.
}
\keyword{internal}
