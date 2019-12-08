% File man/CheckInput.i_BnD.Rd
\name{CheckInput.i_BnD}
\alias{CheckInput.i_BnD}
\title{Checks the input params of a generic Birth and Death process model}
\description{
  Checks the input params of a generic Birth and Death process model
}
\usage{
  \method{CheckInput}{i_BnD}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_BnD}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a generic Birth and Death process model. The inputs params are created calling previously the \link{NewInput.BnD}}

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

## Check the inputs
CheckInput(i_BnD)

}
\keyword{BnD}
