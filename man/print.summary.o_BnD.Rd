% File man/print.summary.o_BnD.Rd
\name{print.summary.o_BnD}
\alias{print.summary.o_BnD}
\title{Summary of the results of a generic Birth and Death process model}
\description{
  Summary of the results of a generic Birth and Death process model.
}
\usage{
  \method{print}{summary.o_BnD}(x, \dots)
}
\arguments{
  \item{x}{a object of class summary.o_BnD}
  \item{\dots}{aditional arguments}
}

\details{Summaries a generic Birth and Death process model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_BnD}}.
}


\examples{
## Generating a generic Birth and Death model with the same lambda and mu vectors as M/M/1 model
## create input parameters
lambda <- rep(1/4, 200)
mu <- rep(1/3, 200)

i_BnD <- NewInput.BnD(lambda=lambda, mu=mu)

## Build the model
o_BnD <- QueueingModel(i_BnD)

## Report the results
print(summary(o_BnD))
}
\keyword{BnD}
