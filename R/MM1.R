############################################################
############################################################
## MODEL M/M/1
############################################################
############################################################
NewInput.MM1 <- function(lambda=0, mu=0, n=0)
{
  res <- list(lambda = lambda, mu = mu, n = n)
  class(res) <- "i_MM1"
  res
}

CheckInput.i_MM1 <- function(x, ...)
{
  MM1_ro_warning <- "ro is greater or equal to one!!"
  MM1_class <- "the class of the object x has to be M/M/1 (i_MM1)"
  MM1_anomalous <- "Some value of lambda, mu or n is anomalous. Check the values." 

 if (!inherits(x, "i_MM1"))
  stop(MM1_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$n))
    stop(MM1_anomalous)

 if (x$mu <= 0)
 	stop(ALL_mu_positive)

 if (x$lambda < 0)
	stop(ALL_lambda_zpositive)

 if (!is.wholenumber(x$n))
  stop(ALL_n_integer)

 ro <- x$lambda / x$mu
 if (ro >= 1)
 {
	 cat(paste("Throughput is ", x$mu, "\n", sep=""))
	 cat(paste("Utilization is ", ro * 100, "%\n", sep=""))
	 stop(MM1_ro_warning)
 }
}


QueueingModel.i_MM1 <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MM1(x, ...)

  # variables to improve the eficiency of the computing
  aux <- (x$mu - x$lambda)
  aux1 <- (x$mu * aux)
  
  RO <- x$lambda / x$mu
  Lq <- (x$lambda^2) / aux1
  Wq <- x$lambda / aux1
  L <- x$lambda / aux
  W <- 1 / aux
  Lqq <- x$mu / aux
  Throughput <- x$lambda

  # Variance
  VNq <- (RO^2 * (1 + RO - RO^2)) / ((1 - RO)^2)
  VTq <- (RO * (2 - RO)) / (x$mu^2 * ((1 - RO)^2))
  VN  <- RO / ((1 - RO)^2)
  VT  <- W^2

  if (x$n < 0)
    Pn <- numeric()
  else
    Pn <- sapply(seq(0, x$n, 1), function(i){dgeom(i, 1-RO)})

  # The distribution functions
  FWq <- function(t) { 1 - (RO * exp(-t/W)) }
  FW <- function(t) { 1 - exp(-t/W) }

  res <- list(
    Inputs = x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput, L = L, VN = VN,
    W = W, VT = VT, Wqq = W, Lqq = Lqq, Pn = Pn, Qn = Pn, FW = FW, FWq = FWq 
  )

  class(res) <- "o_MM1"
  res
} 

RO.o_MM1         <- function(x, ...) { x$RO }
Pn.o_MM1         <- function(x, ...) { x$Pn }
Qn.o_MM1         <- function(x, ...) { x$Qn }
Lq.o_MM1         <- function(x, ...) { x$Lq }
VNq.o_MM1        <- function(x, ...) { x$VNq }
Wq.o_MM1         <- function(x, ...) { x$Wq }
VTq.o_MM1        <- function(x, ...) { x$VTq }
L.o_MM1          <- function(x, ...) { x$L }
VN.o_MM1         <- function(x, ...) { x$VN }
W.o_MM1          <- function(x, ...) { x$W }
VT.o_MM1         <- function(x, ...) { x$VT }
Wqq.o_MM1        <- function(x, ...) { x$Wqq }
Lqq.o_MM1        <- function(x, ...) { x$Lqq }
Inputs.o_MM1     <- function(x, ...) { x$Inputs }
Throughput.o_MM1 <- function(x, ...) { x$Throughput }

Report.o_MM1 <- function(x, ...)
{ 
  reportAux(x)
}

summary.o_MM1 <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MM1  <- function(x, ...)
{
  print_summary(x, ...)
}
