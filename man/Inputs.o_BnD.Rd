% File man/Inputs.o_BnD.Rd
\name{Inputs.o_BnD}
\alias{Inputs.o_BnD}
\title{Returns the input parameters of a generic Birth and Death process model}
\description{
  Returns the inputs parameters of a already built Birth and Death process model
}
\usage{
  \method{Inputs}{o_BnD}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_BnD}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a generic Birth and Death process model. The inputs parameters are created calling previously the \link{NewInput.BnD}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.BnD}}.
}


\examples{
## Generating a generic Birth and Death model with the same lambda and mu vectors as M/M/1 model
## create input parameters
lambda <- rep(1/4, 200)
mu <- rep(1/3, 200)

i_BnD <- NewInput.BnD(lambda=lambda, mu=mu)

## Build the model
o_BnD <- QueueingModel(i_BnD)

## It returns the Inputs
Inputs(o_BnD)
}

\keyword{M/M/1}
