% File man/print.summary.o_MMCC.Rd
\name{print.summary.o_MMCC}
\alias{print.summary.o_MMCC}
\title{Summary of the results of a M/M/c/c queueing model}
\description{
  Summary of the results of a M/M/c/c queueing model.
}
\usage{
  \method{print}{summary.o_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class summary.o_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Summaries a M/M/c/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCC}}.
}


\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
o_mmcc <- QueueingModel(i_mmcc)

## Report the results
print(summary(o_mmcc))
}
\keyword{M/M/c/c}
