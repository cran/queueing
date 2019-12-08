############################################################
############################################################
## MODEL M/M/Infinite
############################################################
############################################################
NewInput.MMInf <- function(lambda=0, mu=0, n=0)
{
  res <- list(lambda = lambda, mu = mu, n = n)
  class(res) <- "i_MMInf"
  res
}


CheckInput.i_MMInf <- function(x, ...)
{
  MMInf_class <- "The class of the object x has to be M/M/Inf (i_MMInf)"
  MMInf_anomalous <- "Some value of lambda, mu, or n is anomalous. Check the values."

  if (!inherits(x, "i_MMInf"))
   	stop(MMInf_class)	

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$n))
    stop(MMInf_anomalous)

  if (x$mu <= 0)
 		stop(ALL_mu_positive)

 	if (x$lambda < 0)
		stop(ALL_lambda_zpositive)

  if (!is.wholenumber(x$n))
    stop(ALL_n_integer)
}


QueueingModel.i_MMInf <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMInf(x, ...)

  # Calculate the output parameters of the model 
  W <- 1 / x$mu
  L <- x$lambda * W
  
  Throughput <- x$lambda

  # we're going to calculate the probability distribution
  if (x$n < 0)
    Pn <- numeric()
  else
    Pn <- sapply(0:x$n, dpois, L)

  FW <- function(t){ exp(x$mu) }
  FWq <- function(t){ 0 }

  Lq  <- 0
  Wq  <- 0
  Lqq <- NA
  Wqq <- NA

  VN  <- L
  VNq <- 0
  VT  <- 1/(x$mu^2)
  VTq <- 0  

  # The result
  res <- list(
    Inputs=x, RO = L, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, VN = VN, W = W, VT = VT, Lqq = Lqq, Wqq = Wqq, Pn = Pn, Qn = Qn, FW = FW, FWq = FWq
  )

  class(res) <- "o_MMInf"
  res

} 

Inputs.o_MMInf <- function(x, ...) { x$Inputs }
L.o_MMInf          <- function(x, ...) { x$L }
VN.o_MMInf         <- function(x, ...) { x$VN }
W.o_MMInf          <- function(x, ...) { x$W }
VT.o_MMInf         <- function(x, ...) { x$VT }
RO.o_MMInf         <- function(x, ...) { x$RO }
Lq.o_MMInf         <- function(x, ...) { x$Lq }
VNq.o_MMInf        <- function(x, ...) { x$VNq }
Wq.o_MMInf         <- function(x, ...) { x$Wq }
VTq.o_MMInf        <- function(x, ...) { x$VTq }
Wqq.o_MMInf        <- function(x, ...) { x$Wqq }
Lqq.o_MMInf        <- function(x, ...) { x$Lqq }
Pn.o_MMInf         <- function(x, ...) { x$Pn }
Qn.o_MMInf         <- function(x, ...) { x$Qn }
Throughput.o_MMInf <- function(x, ...) { x$Throughput }

Report.o_MMInf <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MMInf <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMInf  <- function(x, ...)
{
  print_summary(x, ...)
}
