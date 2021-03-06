% File man/SP.Rd
\name{SP}
\alias{SP}
\title{Returns the saturation point of a queueing model}
\description{
  Returns the saturation point of a queueing model
}
\usage{
  SP(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Returns the saturation point of a queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{SP.o_MM1KK}}
}


\examples{
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=4, method=3)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Returns the saturation point
SP(o_mm1kk)

}
\keyword{M/M/1/K/K}
