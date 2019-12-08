############################################################
############################################################
## MODEL M/M/Infinite/K/K
############################################################
############################################################
NewInput.MMInfKK <- function(lambda=0, mu=0, k=1)
{
  res <- list(lambda = lambda, mu = mu, k = k)
  class(res) <- "i_MMInfKK"
  res
}


CheckInput.i_MMInfKK <- function(x, ...)
{
  MMInfKK_class <- "The class of the object x has to be M/M/Inf/K/K (i_MMInfKK)"
  MMInfKK_anomalous <- "Some value of lambda, mu, or n is anomalous. Check the values."

  if (!inherits(x, "i_MMInfKK"))
   	stop(MMInfKK_class)	

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$k))
    stop(MMInfKK_anomalous)

  if (x$mu <= 0)
 		stop(ALL_mu_positive)

 	if (x$lambda < 0)
		stop(ALL_lambda_zpositive)

  if (x$k < 0)
		stop(ALL_k_warning)

  if (!is.wholenumber(x$k))
   stop(ALL_k_integer)
}


MMInfKK_InitPn_Aprox_Aux <- function(n, lambda, mu, c, k, m)
{
  (n * (log(lambda) - log(mu))) + (lfactorial(k) - lfactorial(k-n) - lfactorial(n))
}


MMInfKK_InitPn <- function(x)
{
  ProbFactCalculus(
    x$lambda, x$mu, 1, x$k, x$k, x$k, MMInfKK_InitPn_Aprox_Aux, MMInfKK_InitPn_Aprox_Aux, MMInfKK_InitPn_Aprox_Aux
  )
}



QueueingModel.i_MMInfKK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMInfKK(x, ...)

  # we're going to calculate the probability distribution  
  Pn <- MMInfKK_InitPn(x)

  u <- x$lambda/x$mu

  # Calculate the output parameters of the model
  L <- (x$k * u)/(1 + u)
  Throughput <- x$lambda * (x$k - L) 
  W <- L / Throughput

  Lq  <- 0
  VNq <- 0
  Wq  <- 0
  VTq <- 0
  Wqq <- NA
  Lqq <- NA

  QnAux <- function(n){ Pn[n] * (x$k - (n-1)) / (x$k - L) }
  Qn <- sapply(1:x$k, QnAux)

  FW <- function(t){ exp(x$mu) }
  FWq <- function(t){ 0 }

  # if the sum(Pn) == 0, then too big K or lambda/mu is
  if (sum(Pn) == 0)
  {
    VN <- NA
  }
  else
  {
    VT <- ( ((0:x$k)^2) * Pn) - (L^2)
  }

  VN <- 1/(x$mu^2)

  # The result
  res <- list(
    Inputs=x, RO = L, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, W = W, Lqq = Lqq, Wqq = Wqq, Pn = Pn, Qn = Qn, FW = FW, FWq = FWq 
  )

  class(res) <- "o_MMInfKK"
  res

} 

Inputs.o_MMInfKK     <- function(x, ...) { x$Inputs }
L.o_MMInfKK          <- function(x, ...) { x$L }
VN.o_MMInfKK         <- function(x, ...) { x$VN }
W.o_MMInfKK          <- function(x, ...) { x$W }
VT.o_MMInfKK         <- function(x, ...) { x$VT }
RO.o_MMInfKK         <- function(x, ...) { x$RO }
Lq.o_MMInfKK         <- function(x, ...) { x$Lq }
VNq.o_MMInfKK        <- function(x, ...) { x$VNq }
Wq.o_MMInfKK         <- function(x, ...) { x$Wq }
VTq.o_MMInfKK        <- function(x, ...) { x$VTq }
Wqq.o_MMInfKK        <- function(x, ...) { x$Wqq }
Lqq.o_MMInfKK        <- function(x, ...) { x$Lqq }
Pn.o_MMInfKK         <- function(x, ...) { x$Pn }
Qn.o_MMInfKK         <- function(x, ...) { x$Qn }
Throughput.o_MMInfKK <- function(x, ...) { x$Throughput }


Report.o_MMInfKK <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MMInfKK <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMInfKK  <- function(x, ...)
{
  print_summary(x, ...)
}
