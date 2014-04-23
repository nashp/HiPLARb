\name{HiPLARb Parameters}
\alias{hiplarb_get_param}
\alias{hiplarb_set_param}
\title{
Get and set HiPLARb Optimisation Parameters
}
\description{
Functions to get and set the crossover and maximum size optimisation parameters.
}
\usage{
hiplarb_get_param()
hiplarb_set_param(param)
}
\arguments{\item{param}{Parameter data frame as returnend by hiplarb_get_param.}}
\details{
hiplarb_get_param returns a data frame with the current parameters. This data frame can be modified and the new values set with hiplarb_set_param.
}
\value{
    hiplarb_get_param returns a data frame with the optimisation parameters.
}
\author{
Vendel Szeremi
}
\note{
The crossover parameters are determined by hiplarb_do_xover_timings and need not be modified.
}


\seealso{
\code{\link{hiplarb_do_xover_timings}}
}
\examples{
hiplarb_get_param()

p <- hiplarb_get_param()
p["dgesv", "maxmagma"] <- 10000
hiplarb_set_param(p)
}
