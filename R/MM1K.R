############################################################
############################################################
## MODEL M/M/1/K - Capacity limited of the system         ##
############################################################
############################################################
NewInput.MM1K <- function(lambda=0, mu=0, k=1)
{
  res <- list(lambda = lambda, mu = mu, k = k)
  class(res) <- "i_MM1K"
  res
}

CheckInput.i_MM1K <- function(x, ...)
{
  MM1K_k_one <- "k must be equal or greater than one"
  MM1K_class <- "the class of the object x has to be M/M/1/K (i_MM1K)"
  MM1K_anomalous <- "Some value of lambda, mu, or k is anomalous. Check the values."

  if (class(x) != "i_MM1K")
   	stop(MM1K_class)	

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$k))
    stop(MM1K_anomalous)

  if (x$mu <= 0)
 	  stop(ALL_mu_positive)

  if (x$lambda < 0)
	  stop(ALL_lambda_zpositive)

  if (x$k < 1)
	  stop(ALL_k_warning)

  if (!is.wholenumber(x$k))
    stop(ALL_k_integer)
 
}


MM1K_InitPn <- function(x)
{

  # pn <- numeric()
  
  if (x$lambda == x$mu)
  {
    # pn[1:(x$k+1)] <- 1 / (x$k + 1)
    pn <- rep(1 / (x$k + 1), x$k+1)
  }
  else if (x$lambda < x$mu)
  {
    one_minus_u <- 1 - ( x$lambda / x$mu ) 
    pn <- dgeom(0:x$k, one_minus_u)/pgeom(x$k, one_minus_u)
  }
  else # x$lambda > x$mu
  {
    pow <- function(e, b, k){k * (b^e)}
    u <- x$lambda / x$mu
    aux <- (1 - u) / (1 - (u^(x$k+1)))
    pn <- sapply(0:x$k, pow, u, aux)
  }
	pn
}


MM1K_L <- function(x)
{
 if (x$lambda == x$mu) ( x$k / 2 )
 else
 {
		u <- x$lambda / x$mu
    u_up_k <- u^(x$k)
    u_up_k_plus_1 <- u_up_k * u
    numerator <- x$lambda * (1 - ((x$k + 1) * u_up_k) + (x$k * u_up_k_plus_1))
    denominator <- (x$mu - x$lambda) * (1 - u_up_k_plus_1)
		numerator / denominator
 }
}

QueueingModel.i_MM1K <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MM1K(x, ...)

  Pn <- MM1K_InitPn(x)

  RO <- 1 - Pn[1]

  L <- MM1K_L(x)
  Lq <- L - RO
  Throughput <- x$lambda * (1 - Pn[x$k+1])
  W <- L / Throughput
  Wq <- Lq / Throughput
  #Wqq <- Wq / RO
  Qn <- Pn[seq(1, x$k, 1)] / (1 - Pn[x$k+1])
  
  if (x$k == 1)
  {
    Lqq <- NA
    Wqq <- NA
  }
  else
  {
    Wqq <- Wq / (1 - Qn[1])
    Lqq <- Wqq * x$mu    
  }
     
  FW <- function(t){
    aux <- function(i) { Qn[i] * ppois(i-1, x$mu * t) }
    1 - sum(sapply(seq(1, x$k, 1), aux))
  }

  if (x$k == 1)
    FWq <- function(t){0}
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

  FWInt  <- integrate(xFWc, 0, Inf)

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
    L = L, VN = VN, W = W, VTq = VTq, Wqq = Wqq, Lqq = Lqq, Pn = Pn, Qn = Qn, FW = FW, FWq = FWq)

  class(res) <- "o_MM1K"
  res
} 

Inputs.o_MM1K     <- function(x, ...) { x$Inputs }
L.o_MM1K          <- function(x, ...) { x$L }
VN.o_MM1K         <- function(x, ...) { x$VN }
W.o_MM1K          <- function(x, ...) { x$W }
VT.o_MM1K         <- function(x, ...) { x$VT }
RO.o_MM1K         <- function(x, ...) { x$RO }
Lq.o_MM1K         <- function(x, ...) { x$Lq }
VNq.o_MM1K        <- function(x, ...) { x$VNq }
Lqq.o_MM1K        <- function(x, ...) { x$Lqq }
Wq.o_MM1K         <- function(x, ...) { x$Wq }
VTq.o_MM1K        <- function(x, ...) { x$VTq }
Wqq.o_MM1K        <- function(x, ...) { x$Wqq }
Pn.o_MM1K         <- function(x, ...) { x$Pn }
Qn.o_MM1K         <- function(x, ...) { x$Qn }
Throughput.o_MM1K <- function(x, ...) { x$Throughput }


Report.o_MM1K <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MM1K <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MM1K  <- function(x, ...)
{
  print_summary(x, ...)
}
