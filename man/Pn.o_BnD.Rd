% File man/Pn.o_BnD.Rd
\name{Pn.o_BnD}
\alias{Pn.o_BnD}
\alias{Qn.o_BnD}
\title{Returns the probabilities of a generic Birth and Death process model}
\description{
  Pn returns the probabilities that a generic Birth and Death process model has n customers.
}
\usage{
  \method{Pn}{o_BnD}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_BnD}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a generic Birth and Death process model has n customers.}

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

## Returns the probabilities
Pn(o_BnD)

## Simulating M/M/1
lambda <- rep(1/4, 200)
mu <- rep(1/3, 200)

pn_bnd_mm1 <- Pn(QueueingModel(NewInput.BnD(lambda=lambda, mu=mu)))
pn_mm1 <- Pn(QueueingModel(NewInput.MM1(lambda=1/4, mu=1/3, n=200)))

## Simulating M/M/2
lambda <- rep(5, 200)
mu <- c(1*10, rep(2*10, 199))

pn_mmc <- Pn(QueueingModel(NewInput.MMC(lambda=5, mu=10, c=2, n=200, method=0)))
pn_bnd_mmc <- Pn(QueueingModel(NewInput.BnD(lambda=lambda, mu=mu)))


## Simulating M/M/1/K/K
lambda <- c(2*0.25, 0.25)
mu <- rep(4, 2)
pn_mm1kk <- Pn(QueueingModel(NewInput.MM1KK(lambda=0.25, mu=4, k=2, method=3)))
pn_bnd <- Pn(QueueingModel(NewInput.BnD(lambda=lambda, mu=mu)))

}
\keyword{BnD}
