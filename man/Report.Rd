% File man/Report.Rd
\name{Report}
\alias{Report}
\title{Reports the results of a queueing model}
\description{
  Reports the results of a queueing model.
}

\usage{
  Report(x, \dots)
}

\arguments{
  \item{x}{i_MM1, i_MMC, i_MM1K, i_MMCK, i_MM1KK, i_MMCKK, i_MMCC, i_MMCKM, i_MMInfKK, i_MMInf, i_OJN, i_MCON}
  \item{\dots}{aditional arguments}
}

\details{Generic S3 method to report a queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel}}.
}


\examples{
## See example 10.3 in reference [Sixto2004] for more details.
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Report the results
Report(o_mm1)
}
\keyword{M/M/1}
