\name{HiPLARb-package}
\alias{HiPLARb-package}
\alias{HiPLARb}
\alias{hiplarb}
\docType{package}
\title{
HiPLARb - High Performance Linear Algebra in R base
}
\description{
This package provides multi-core CPU and GPU support for linear algebra in R, by overriding functions in the R base library and the R evaluator with modified versions that use the PLASMA and MAGMA libraries.
}
\details{
\tabular{ll}{
Package: \tab HiPLARb\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2012-12-10\cr
License: \tab GPL (>=2)\cr
}

When loading the HiPLARb library it overrides functions in the R evaluator and in the R base library with versions using PLASMA and/or MAGMA. The actual function interfaces remain unmodified and R scripts can make use of multi-core CPUs and GPUs without modifications.

}
\author{
Vendel Szeremi, HiPLAR Team

Maintainer: Giovanni Montana  <support@hiplar.org>
}
\note{
Using HiPLARb on a machine with multiple concurrent users with an oversubscription of processes, may not give any benefits.
}
\references{
\url{http://icl.cs.utk.edu/plasma}

\url{http://icl.cs.utk.edu/magma}
}


\keyword{ package }
\seealso{
\code{\link[base:base-package]{base}} for help regarding the actual functions.

\code{\link{hiplarb_mode_plasma}},
\code{\link{hiplarb_mode_magma}}, and
\code{\link{hiplarb_mode_auto}} for changing the mode of operation when installed with both PLASMA and MAGMA.

\code{\link{hiplarb_do_xover_timings}},
\code{\link{hiplarb_get_param}}, and
\code{\link{hiplarb_set_param}} for changing HiPLARb parameters.
}
\examples{
library(HiPLARb)
}
