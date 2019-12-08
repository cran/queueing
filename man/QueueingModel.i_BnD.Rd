% File man/QueueingModel.i_BnD.Rd
\name{QueueingModel.i_BnD}
\alias{QueueingModel.i_BnD}
\title{Builds a a generic Birth and Death process model}
\description{
  Builds a a generic Birth and Death process model
}
\usage{
  \method{QueueingModel}{i_BnD}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_BnD}
  \item{\dots}{aditional arguments}
}

\details{Build a generic Birth and Death process model. It also checks the input params calling the \link{CheckInput.i_BnD}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_BnD}}
}


\examples{
## Generating a generic Birth and Death model with the same lambda and mu vectors as M/M/1 model
## create input parameters
lambda <- rep(1/4, 200)
mu <- rep(1/3, 200)

i_BnD <- NewInput.BnD(lambda=lambda, mu=mu)

## Build the model
o_BnD <- QueueingModel(i_BnD)

}
\keyword{BnD}
