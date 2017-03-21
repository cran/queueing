% File man/print.summary.o_MMCKK.Rd
\name{print.summary.o_MMCKK}
\alias{print.summary.o_MMCKK}
\title{Summary of the results of a M/M/c/K/K queueing model}
\description{
  Summary of the results of a M/M/c/K/K queueing model.
}
\usage{
  \method{print}{summary.o_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class summary.o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Summaries a M/M/c/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKK}}.
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12, method=0)

## Build the model
o_mmckk <- QueueingModel(i_mmckk)

## Report the results
print(summary(o_mmckk))
}
\keyword{M/M/c/K/K}

