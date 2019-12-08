###############################################################
###############################################################
## MODEL M/M/c/K - Capacity limited of the system, c servers.##
###############################################################
###############################################################
NewInput.MMCK <- function(lambda=0, mu=0, c=1, k=1)
{
  res <- list(lambda = lambda, mu = mu, c = c, k = k)
  class(res) <- "i_MMCK"
  res
}

CheckInput.i_MMCK <- function(x, ...)
{
  MMCK_class <- "the class of the object x has to be M/M/C/K (i_MMCK)"
  MMCK_anomalous <- "Some value of lambda, mu, c or k is anomalous. Check the values."

 if (!inherits(x, "i_MMCK"))
   	stop(MMCK_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$k)
  )
    stop(MMCK_anomalous)

 if (x$lambda < 0)
	stop(ALL_lambda_zpositive)

 if (x$mu <= 0)
 	stop(ALL_mu_positive)

 if (x$c < 1)
	stop(ALL_c_warning)

 if (!is.wholenumber(x$c))
   stop(ALL_c_integer)

 if (x$k < 1)
	stop(ALL_k_warning)

 if (!is.wholenumber(x$k))
   stop(ALL_k_integer)

 if (x$k < x$c)
	stop(ALL_k_c)
}


MMCK_InitPn <- function(x)
{
  
  pn <- numeric()

  u <- x$lambda / x$mu
  ro <- u / x$c

  prod <- 1
	acum <- 1
  aux <- 1

  i <- 1
  pn[i] <- prod # in the final, to multiply by p0

	while (i <= x$c-1)
  {
		prod <- prod * u/i
   	acum <- acum + prod
    pn[i+1] <- prod
    i <- i + 1
  }  

  prod <- prod * ro # this is the case of i = c
  pn[x$c+1] <- prod

  if (ro == 1)
		p0 <- 1 / (acum + (prod * (x$k - x$c + 1)))
	else
		p0 <- 1 / (acum + (prod * ((1 - ro^(x$k - x$c + 1)) / (1 - ro))))

  # from c+1 to k
  i <- x$c + 1

  while (i <= x$k)
  {
    prod <- prod * u/x$c
    pn[i+1] <- prod
    i <- i + 1
  }

  p0 * pn
}


QueueingModel.i_MMCK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMCK(x, ...)


  Pn <- MMCK_InitPn(x)
 
 	aux <- x$lambda / (x$c * x$mu)
	queue_max_length <- x$k - x$c

  Lq <-
    if (aux == 1)
		  Pn[x$c+1] * queue_max_length * (queue_max_length + 1) / 2
	  else
	  {
      one_minus_aux <- 1 - aux
		  aux_up_queue_max_length <- aux^queue_max_length
		  aux_up_queue_max_length_plus_1 <- aux_up_queue_max_length * aux
		  tmp1 <- 1 - aux_up_queue_max_length_plus_1 - ((queue_max_length+1) * aux_up_queue_max_length * one_minus_aux)
		  tmp2 <- one_minus_aux^2
		  Pn[x$c+1] * aux * (tmp1 / tmp2)
	  }

  Throughput <- x$lambda * (1 - Pn[x$k+1])

  L <- Lq + (Throughput / x$mu)

  RO <- Throughput / (x$mu * x$c)
  W <- L / Throughput
  Wq <- W - (1/x$mu)
  Qn <- Pn[seq(1, x$k, 1)] / (1 - Pn[x$k+1])

  if (x$c == x$k)
  {
    Wqq <- NA
    Lqq <- NA
  }
  else
  {
    Wqq <- Wq / (1 - sum(Qn[1:x$c]))
    Lqq <- Wqq * (x$c * x$mu)
  }

  auxFwq <- function(n, t) { Qn[n] * ppois(n-x$c-1, x$c * x$mu * t) }

  if (x$c == x$k)
    FWq <- 0
  else
    FWq <- function(t){
      1 - sum(sapply(seq(x$c+1, x$k, 1), auxFwq, t))
    }

  # variances
  VN  <- sum( (0:x$k)^2 * Pn ) - (L^2)

  if (x$c == x$k)
    VNq <- 0
  else
    VNq <- sum( ( c( rep(0, x$c+1), 1:(x$k-x$c) )^2 * Pn) - (Lq^2) )

  xFWqc <- function(t){Vectorize(t * (1 - FWq(t)))}

  if (x$c == x$k)
    VTq <- 0
  else
  {
    FWqInt  <- integrate(xFWqc, 0, Inf)
    
    if (FWqInt$message == "OK")
      VTq <- (2 * FWqInt$value) - (Wq^2) 
    else
      VTq <- NA
  }
  
  # The result
  res <- list(
    Inputs = x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, VN = VN, W = W, Lqq = Lqq, Wqq = Wqq, Pn = Pn, Qn = Qn, FWq = FWq
  )

  class(res) <- "o_MMCK"
  res

} 

Inputs.o_MMCK     <- function(x, ...) { x$Inputs }
L.o_MMCK          <- function(x, ...) { x$L }
VN.o_MMCK         <- function(x, ...) { x$VN }
W.o_MMCK          <- function(x, ...) { x$W }
RO.o_MMCK         <- function(x, ...) { x$RO }
Lq.o_MMCK         <- function(x, ...) { x$Lq }
VNq.o_MMCK        <- function(x, ...) { x$VNq }
Lqq.o_MMCK        <- function(x, ...) { x$Lqq }
Wq.o_MMCK         <- function(x, ...) { x$Wq }
VTq.o_MMCK        <- function(x, ...) { x$VTq }
Wqq.o_MMCK        <- function(x, ...) { x$Wqq }
Pn.o_MMCK         <- function(x, ...) { x$Pn }
Qn.o_MMCK         <- function(x, ...) { x$Qn }
Throughput.o_MMCK <- function(x, ...) { x$Throughput }


Report.o_MMCK <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MMCK <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMCK  <- function(x, ...)
{
  print_summary(x, ...)
}
