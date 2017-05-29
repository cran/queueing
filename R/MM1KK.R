###############################################################
###############################################################
## MODEL M/M/1/K/K - Finite Poblation.                       ##
###############################################################
###############################################################
NewInput.MM1KK <- function(lambda=0, mu=0, k=1, method=3)
{
  res <- list(lambda = lambda, mu = mu, k = k, method = method)
  class(res) <- "i_MM1KK"
  res
}

CheckInput.i_MM1KK <- function(x, ...)
{

  MM1KK_class     <- "The class of the object x has to be M/M/1/K/K (i_MM1KK)"
  MM1KK_anomalous <- "Some value of lambda, mu or k is anomalous. Check the values."
  MM1KK_method    <- "method variable has to be 0 to be exact calculus, 1 to be aproximate calculus, 2 to use Jain's Method or 3 to use Poisson truncated distribution"


 if (class(x) != "i_MM1KK")
   stop(MM1KK_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$k))
   stop(MM1KK_anomalous)

 if (x$lambda < 0)
   stop(ALL_lambda_zpositive)

 if (x$mu <= 0)
 	 stop(ALL_mu_positive)

 if (x$k < 1)
   stop(ALL_k_warning)

 if (!is.wholenumber(x$k))
		stop(ALL_k_integer)

 if (x$method != 0 && x$method != 1 && x$method != 2 && x$method != 3)
   stop(MM1KK_method)
}


MM1KK_InitPn_Aprox_Aux <- function(n, lambda, mu, c, k, m)
{
  (lfactorial(k) - lfactorial(k-n)) + (n * log(lambda/mu))
}


MM1KK_InitPn_Aprox <- function(x)
{
  ProbFactCalculus(x$lambda, x$mu, 1, x$k, x$k, x$k,
    MM1KK_InitPn_Aprox_Aux, MM1KK_InitPn_Aprox_Aux, MM1KK_InitPn_Aprox_Aux)
}


MM1KK_InitPn_Exact <- function(x)
{
  pn <- c(0:x$k)

  z <- x$mu / x$lambda
	u <- x$lambda / x$mu
  
  pn[1] <- B_erlang(x$k, z)
  
  totu <- 1
	totk <- 1

  i <- 2
  while (i <= (x$k + 1))
  {
    totu <- totu * u
		totk <- totk * ((x$k + 1) - i + 1)
    pn[i] <- pn[1] * totu * totk 
		i <- i + 1
  }	

  pn
}


MM1KK_method2_Aux <- function(x, i)
{
  r <- x$lambda/x$mu

  if (i == 0)
  {
    (x$k - i) * r / (i+1)
  }
  else
  {
    (x$k - i) * r
  }
}


MM1KK_method2_Prod <- function(x,n)
{
  prod <- 1
  
  for (i in 0:(n-1))
  {
    prod <- prod * MM1KK_method2_Aux(x, i)
  }

  prod

}


MM1KK_method2_Prob <- function(x)
{

  pn <- c()
  
  sumAux <- 1

  for (i in (1:x$k))
  {
    sumAux <- sumAux + MM1KK_method2_Prod(x, i)
  }

  pn[1] <- 1/sumAux

  for (i in 2:(x$k+1))
  {
    pn[i] <- MM1KK_method2_Aux(x, i-2) * pn[i-1]
  }

  pn
}


MM1KK_method3_Prob <- function(x)
{
  z <- x$mu/x$lambda

  funMethod3 <- function(n){ dpois(x$k-n, z)/ppois(x$k, z) }

  pn <- sapply(0:x$k, funMethod3)
  pn
}



MM1KK_InitPn <- function(x)
{
  if (x$method == 0)
    pn <- MM1KK_InitPn_Exact(x)
  else if (x$method == 1)
    pn <- MM1KK_InitPn_Aprox(x)
  else if (x$method == 2)
    pn <- MM1KK_method2_Prob(x)
  else
    pn <- MM1KK_method3_Prob(x)

  pn
}


QueueingModel.i_MM1KK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MM1KK(x, ...)
  z <- x$mu / x$lambda

  Pn <- MM1KK_InitPn(x)

  RO <- 1 - Pn[1]
  Throughput <- x$mu * RO
    
  L <- x$k - (Throughput / x$lambda)
  W <- (x$k / Throughput) - ( 1 / x$lambda)
  Wq <- W - (1 / x$mu)
  Lq <- Throughput * Wq
  WWs <- (x$k / RO) - z
  SP <- 1 + z

  QnAux <- function(n){ Pn[n] * (x$k - (n-1)) / (x$k - L) }
  Qn <- sapply(1:x$k, QnAux)

  if (x$k == 1)
  {
    Wqq <- NA
    Lqq <- NA
  }
  else
  {
    Wqq <- Wq / (1 - Qn[1])
    Lqq <- Wqq * x$mu    
  }

  #Wqq <- Wq / RO
  
  FW <- function(t){
    aux <- function(i, t) { Qn[i] * ppois(i-1, x$mu * t) }
    1 - sum(sapply(seq(1, x$k, 1), aux, t))
  }

  if (x$k == 1)
    FWq <- function(t){ 0 }
  else
  {
    FWq <- function(t){
      aux <- function(i, t) { Qn[i+1] * ppois(i-1, x$mu * t) }
      1 - sum(sapply(seq(1, x$k-1, 1), aux, t))
    }
  }
 
  # variances
  VN  <- sum( (0:x$k)^2 * Pn ) - (L^2)  

  xFWc  <- Vectorize(function(t){t * (1 - FW(t))})
  xFWqc <- Vectorize(function(t){t * (1 - FWq(t))})

  FWInt <- integrate(xFWc, 0, Inf)

  if (FWInt$message == "OK")
    VT <- (2 * FWInt$value) - (W^2) 
  else
    VT <- NA

  if (x$k == 1)
  {
    VNq <- 0
    VTq <- 0
  }
  else
  {
    VNq <- sum( c(0, 0, 1:(x$k-1))^2 * Pn ) - (Lq^2)
    FWqInt <- integrate(xFWqc, 0, Inf)
    
    if (FWqInt$message == "OK")
      VTq <- (2 * FWqInt$value) - (Wq^2) 
    else
      VTq <- NA
  }

  # The result
  res <- list(
    Inputs=x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, VN = VN, W = W, VT = VT, Lqq = Lqq, Wqq = Wqq, WWs = WWs, SP = SP, Pn = Pn, Qn = Qn, FW = FW, FWq = FWq
  )
  
  class(res) <- "o_MM1KK"
  res

} 

Inputs.o_MM1KK     <- function(x, ...) { x$Inputs }
RO.o_MM1KK         <- function(x, ...) { x$RO }
Lq.o_MM1KK         <- function(x, ...) { x$Lq }
VNq.o_MM1KK        <- function(x, ...) { x$VNq }
Wq.o_MM1KK         <- function(x, ...) { x$Wq }
VTq.o_MM1KK        <- function(x, ...) { x$VTq }
L.o_MM1KK          <- function(x, ...) { x$L }
VN.o_MM1KK         <- function(x, ...) { x$VN }
W.o_MM1KK          <- function(x, ...) { x$W }
VT.o_MM1KK         <- function(x, ...) { x$VT }
Lqq.o_MM1KK        <- function(x, ...) { x$Lqq }
Wqq.o_MM1KK        <- function(x, ...) { x$Wqq }
WWs.o_MM1KK        <- function(x, ...) { x$WWs }
SP.o_MM1KK         <- function(x, ...) { x$SP }
Pn.o_MM1KK         <- function(x, ...) { x$Pn }
Qn.o_MM1KK         <- function(x, ...) { x$Qn }
Throughput.o_MM1KK <- function(x, ...) { x$Throughput }


Report.o_MM1KK <- function(x, ...)
{ 
  reportAux(x)  
}


summary.o_MM1KK <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MM1KK  <- function(x, ...)
{
  print_summary(x, ...)
}
