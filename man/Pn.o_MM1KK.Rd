% File man/Pn.o_MM1KK.Rd
\name{Pn.o_MM1KK}
\alias{Pn.o_MM1KK}
\alias{Qn.o_MM1KK}
\title{Returns the probabilities of a M/M/1/K/K queueing model}
\description{
  Pn eeturns the probabilities of a M/M/1/K/K queueing model
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MM1KK}(x, \dots)
  \method{Qn}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/1/K/K queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1KK}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2, method=3)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Returns the probabilities
Pn(o_mm1kk)
Qn(o_mm1kk)
}
\keyword{M/M/1/K/K}
