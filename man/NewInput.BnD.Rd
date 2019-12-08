% File man/NewInput.BnD.Rd
\name{NewInput.BnD}
\alias{NewInput.BnD}
\title{Define the inputs of a new generic Birth and Death process model}
\description{
  Define the inputs of a new generic Birth and Death process model
}

\usage{
NewInput.BnD(lambda=NULL, mu=NULL)
}

\arguments{
  \item{lambda}{vectors of arrival rate depending of the number of users in the system. Observe that in R, the vectors starts counting at 1, so lambda[1] is the arrival rate when the system has 0 users on it}
  \item{mu}{vectors of service rate depending of the number of users in the system. Observe that in R, the vectors starts counting at 1, so mu[1] is the arrival rate when the system has 1 users on it}
}


\details{Define the inputs of a new generic Birth and Death process model}

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
}

\keyword{BnD}
