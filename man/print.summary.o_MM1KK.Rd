% File man/print.summary.o_MM1KK.Rd
\name{print.summary.o_MM1KK}
\alias{print.summary.o_MM1KK}
\title{Summary of the results of a M/M/1/K/K queueing model}
\description{
  Summary of the results of a M/M/1/K/K queueing model.
}
\usage{
  \method{print}{summary.o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class summary.o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Summaries a M/M/1/K/K queueing model}

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

## Report the results
print(summary(o_mm1kk))
}
\keyword{M/M/1/K/K}
