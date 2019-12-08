###############################################################
###############################################################
## MODEL M/M/c/c - Capacity limited of the system, c servers.##
## truncated model, Erlang-B function #########################
###############################################################
###############################################################
NewInput.MMCC <- function(lambda=0, mu=0, c=1, method=1)
{
  res <- list(lambda = lambda, mu = mu, c = c, method = method)
  class(res) <- "i_MMCC"
  res
}

CheckInput.i_MMCC <- function(x, ...)
{
  MMCC_class     <- "the class of the object x has to be M/M/C/C (i_MMCC)"
  MMCC_anomalous <- "Some value of lambda, mu or c is anomalous. Check the values."
  MMCC_method    <- "method variable has to be 0 to be definiton calculus, 1 to be exact calculus"

  if (!inherits(x, "i_MMCC"))
    stop(MMCC_class)

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$c))
     stop(MMCC_anomalous)

  if (x$lambda < 0)
	 stop(ALL_lambda_zpositive)

  if (x$mu <= 0)
 	 stop(ALL_mu_positive)

  if (x$c < 1)
	 stop(ALL_c_warning)

  if (!is.wholenumber(x$c))
    stop(ALL_c_integer)

  if (x$method != 0 && x$method != 1)
   stop(MMCC_method)
}


MMCC_InitPn <- function(x)
{
  if (x$method == 0)
    MMCC_InitPn_def(x)
  else
    MMCC_InitPn_exact(x) 
}


MMCC_InitPn_exact <- function(x)
{
  u <- x$lambda / x$mu  
  tpoisson(0:x$c, x$c, u)  
}


MMCC_InitPn_def <- function(x)
{
  pn <- numeric()

  u <- x$lambda / x$mu
  ro <- u / x$c
  
  prod <- 1
  acum <- 1

  i <- 1
  pn[i] <- prod

  while (i <= x$c-1)
  {
    prod <- prod * u/i
    acum <- acum + prod
    pn[i+1] <- prod    
    i <- i + 1
  }

  prod <- prod * ro
  pn[i+1] <- prod # i has the value c

  p0 <- 1 / (acum + prod)
  pn <- p0 * pn
}


QueueingModel.i_MMCC <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMCC(x, ...)

  Pn <- MMCC_InitPn(x)
  
  Lq <- 0
  Wq <- 0
  Wqq <- NA
  Lqq <- NA

  aux <- x$lambda / x$mu
  one_minus_b_erlang <- 1 - B_erlang(x$c, aux)
  L <- aux * one_minus_b_erlang
  Throughput <- x$lambda * one_minus_b_erlang
  RO <- Throughput / (x$mu * x$c)
  W <- 1 / x$mu 

  Qn <- Pn[1:x$c]/one_minus_b_erlang

  FW <- function(t){
    exp(x$mu)
  }

  FWq <- function(t){0}

  # variances
  VN <- L - (aux * (1 - one_minus_b_erlang) * (x$c - L))
  VT <- 1/(x$mu^2)

  VNq <- 0
  VTq <- 0
 
  # The result
  res <- list(
    Inputs = x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, VN = VN, W = W, VT = VT, Lqq = Lqq, Wqq = Wqq, Pn = Pn, Qn = Qn, FW = FW, FWq = FWq
  )
  
  class(res) <- "o_MMCC"
  res

} 

Inputs.o_MMCC     <- function(x, ...) { x$Inputs }
L.o_MMCC          <- function(x, ...) { x$L }
VN.o_MMCC         <- function(x, ...) { x$VN }
W.o_MMCC          <- function(x, ...) { x$W }
VT.o_MMCC         <- function(x, ...) { x$VT }
RO.o_MMCC         <- function(x, ...) { x$RO }
Lq.o_MMCC         <- function(x, ...) { x$Lq }
VNq.o_MMCC        <- function(x, ...) { x$VNq }
Wq.o_MMCC         <- function(x, ...) { x$Wq }
VTq.o_MMCC        <- function(x, ...) { x$VTq }
Lqq.o_MMCC        <- function(x, ...) { x$Lqq }
Wqq.o_MMCC        <- function(x, ...) { x$Wqq }
Pn.o_MMCC         <- function(x, ...) { x$Pn }
Qn.o_MMCC         <- function(x, ...) { x$Qn }
Throughput.o_MMCC <- function(x, ...) { x$Throughput }

Report.o_MMCC <- function(x, ...)
{
  reportAux(x)  
}


summary.o_MMCC <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMCC  <- function(x, ...)
{
  print_summary(x, ...)
}
