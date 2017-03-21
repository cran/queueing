% File man/print.summary.o_MMInf.Rd
\name{print.summary.o_MMInf}
\alias{print.summary.o_MMInf}
\title{Summary of the results of a M/M/Infinite queueing model}
\description{
  Summary of the results of a M/M/Infinite queueing model.
}
\usage{
  \method{print}{summary.o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class summary.o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Summaries a M/M/Infinite queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMInf}}.
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
o_mminf <- QueueingModel(i_mminf)

## Report the results
print(summary(o_mminf))
}
\keyword{M/M/Infinite}
