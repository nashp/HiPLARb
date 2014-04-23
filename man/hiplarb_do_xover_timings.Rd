\name{HiPLARb Parameter Timings}
\alias{hiplarb_do_xover_timings}
\title{
Perform Timings
}
\description{
For the auto operation mode HiPLARb needs timing data for each function, in order to decide which library to use. This is collected with the function hiplarb_do_xover_timings().
}
\usage{
hiplarb_do_xover_timings(min = 10, max = 3000, rep = 20)
}
\arguments{
  \item{min}{The minimum problem size to be used in the search.}
  \item{max}{The maximum problem size to be used in the search.}
  \item{rep}{The number of repetitions to perform at each step.}
}
\details{
To avoid noisy measurements each timing step is repeated a few times and an average is used. The number of repetitions can be adjusted with the rep argument. The results will be saved to the file ~/.hiplarb.csv and loaded automatically by the HiPLARb package.
}
\value{
Returns a data frame with the new HiPLARb parameters.
}
\author{
Vendel Szeremi
}
\note{
To reduce variability of the measurements, the machine should not otherwise be used while running this function.
}

\seealso{
\code{\link{hiplarb_get_param}}, and
\code{\link{hiplarb_set_param}}.
}
\examples{
hiplarb_do_xover_timings()
}
