% File man/B_erlang.Rd
\name{B_erlang}
\alias{B_erlang}
\title{Returns the probability that all servers are busy}
\description{
  Returns the probability that all servers are busy
}
\usage{
  B_erlang(c=1, u=0)
}
\arguments{
  \item{c}{numbers of servers}
  \item{u}{lambda/mu, that is, ratio of rate of arrivals and rate of service}
}

\details{Returns the probability that all servers are busy}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.\cr
\cr
[Jagerman1974] Jagerman, D. L. (1974).\cr
\emph{Some properties of  the Erlang loss function}.\cr
Bell System Tech. J. (53), 525-551

}

\seealso{
  \code{\link{C_erlang}}
}


\examples{
## two servers
B_erlang(2, 0.5/0.7)
}
\keyword{B_erlang}
