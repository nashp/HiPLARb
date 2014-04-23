\name{HiPLARb Operation Mode}
\alias{hiplarb_mode_plasma}
\alias{hiplarb_mode_magma}
\alias{hiplarb_mode_auto}
\title{
Change the HiPLARb Mode of Operation
}
\description{
Change the HiPLARb operation mode: force the use of PLASMA or MAGMA, or enable auto mode where HiPLARb will decide which library to use depending on the problem size.
}
\usage{
hiplarb_mode_plasma()
hiplarb_mode_magma()
hiplarb_mode_auto()
}
\details{
In auto mode HiPLARb will require timing data. This needs to be collected once with hiplarb_do_xover_timings.
}
\author{
Vendel Szeremi
}
\note{
These functions will only have an effect if HiPLARb is installed with both PLASMA and MAGMA.
}

\seealso{
\code{\link{hiplarb_do_xover_timings}}
\code{\link{hiplarb_get_param}}
\code{\link{hiplarb_set_param}}
}
\examples{
hiplarb_mode_plasma()
hiplarb_mode_magma()
hiplarb_mode_auto()
}
