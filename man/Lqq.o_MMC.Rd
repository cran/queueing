% File man/Lqq.o_MMC.Rd
\name{Lqq.o_MMC}
\alias{Lqq.o_MMC}
\title{Returns the mean number of customers in queue when there is queue in the M/M/c queueing model}
\description{
  Returns the mean number of customers in queue when there is queue in the M/M/c queueing model
}
\usage{
  \method{Lqq}{o_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMC}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of customers in queue when there is queue in the M/M/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMC}}.
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)

## Build the model
o_mmc <- QueueingModel(i_mmc)

## Returns the Lqq
Lqq(o_mmc)
}

\keyword{M/M/c}
