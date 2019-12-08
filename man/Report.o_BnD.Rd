% File man/Report.o_BnD.Rd
\name{Report.o_BnD}
\alias{Report.o_BnD}
\title{Reports the results of a M/M/1 queueing model}
\description{
  Reports the results of a M/M/1 queueing model.
}
\usage{
  \method{Report}{o_BnD}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_BnD}
  \item{\dots}{aditional arguments}
}

\details{Generates a report of the queueing model received as parameter}

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
Report(o_BnD)
}
\keyword{BnD}
