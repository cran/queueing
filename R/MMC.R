############################################################
## Model M/M/C
############################################################
NewInput.MMC <- function(lambda=0, mu=0, c=1, n=0, method=0)
{
  res <- list(lambda = lambda, mu = mu, c = c, n = n, method = method)
  class(res) <- "i_MMC"
  res
}

CheckInput.i_MMC <- function(x, ...)
{
	MMC_r_c_warning <- "( lambda/(mu*c) ) has to be less than one!!"
  MMC_class <- "the class of the object x has to be M/M/C (i_MMC)"
  MMC_anomalous <- "Some value of lambda, mu, c or n is anomalous. Check the values."
  MMC_method <- "method variable has to be 0 to be exact calculus, 1 to be aproximate calculus"

  if (class(x) != "i_MMC")
   	stop(MMC_class)

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$n)
  )
    stop(MMC_anomalous)    

  r <- x$lambda / x$mu  

	if (x$c < 1)
    stop(ALL_c_warning)

	if (x$lambda < 0)
		stop(ALL_lambda_zpositive)

	if (x$mu <= 0)
		stop(ALL_mu_positive)
  
  if (r >= x$c)
  {
    ro <- r/x$c
    cat(paste("Throughput is: ", x$mu * x$c, "\n", sep=""))
    cat(paste("Utilization exceeds 100% use!!: ", ro * 100, "%\n", sep=""))
    stop(MMC_r_c_warning)
  }

  if (!is.wholenumber(x$c))
		stop(ALL_c_integer)

  if (!is.wholenumber(x$n))
		stop(ALL_n_integer)

  if (x$method != 0 && x$method != 1)
    stop(MMC_method)

}


MMC_InitPn_Exact <- function(x)
{
  r <- x$lambda / x$mu
  ro <- r / x$c
  one_minus_ro <- 1 - ro
    
  prod <- 1
  acum <- prod

  if (x$n > x$c)
    pn <- rep(0, x$n)
  else
    pn <- rep(0, x$c)

  i <- 1
	pn[i] <- prod

  while ( i <= (x$c - 1) )
  {
   	prod <- prod * r/i
   	acum <- acum + prod
    pn[i+1] <- prod
   	i <- i + 1
  }

  prod <- prod * r/x$c
  pn[x$c+1] <- prod
    
  p0 <- 1 / (acum + (prod / one_minus_ro))
      
  if (x$n > x$c)
  {
   	for (j in (x$c+1):x$n)
		{
			prod <- prod * r/x$c
      pn[j+1] <- prod
		}
  }    

  # Now, calculate the complete probabilities
  pn <- p0 * pn
  
  # Return the number of elements requested
  pn[1:(x$n+1)]

}

MMC_InitPn_Aprox_AuxToC <- function(n, lambda, mu, c, k, m)
{
  (n * log(lambda/mu)) - lfactorial(n)
}


MMC_InitPn_Aprox_AfterC <- function(n, lambda, mu, c, k, m)
{
  (n * log(lambda/mu)) - lfactorial(c) - (n - c) * log(c)
}


MMC_InitPn_Aprox <- function(x)
{
  (ProbFactCalculus(
      x$lambda, x$mu, x$c, max(x$c, x$n), max(x$c, x$n), max(x$c, x$n),
      MMC_InitPn_Aprox_AuxToC, MMC_InitPn_Aprox_AfterC, MMC_InitPn_Aprox_AfterC
  ))[1:(x$n+1)]
}


MMC_InitPn <- function(x)
{   
  if (x$method == 0)
    MMC_InitPn_Exact(x)
  else # method == 1
     MMC_InitPn_Aprox(x)
}


QueueingModel.i_MMC <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMC(x, ...)

  r <- x$lambda / x$mu
  RO <- r / x$c
  one_minus_ro <- 1 - RO
  inverse_lambda <- 1 / x$lambda
  cErlang <- C_erlang(x$c, r)

  Throughput <- x$lambda

  if (x$n < 0)
    Pn <- numeric()
  else
    Pn <- MMC_InitPn(x)

  Lq <- (cErlang * RO) / (one_minus_ro)
  Wq <- Lq * inverse_lambda
  L <- Lq + r  
  W <- L * inverse_lambda
  Wqq <- 1 / (x$c * one_minus_ro * x$mu)
  Lqq <- Wqq * x$mu * x$c

  FWq <- function(t)
  {
    1 - ( cErlang * exp( (-1) * (1 - RO) * x$c * x$mu * t ) )  
  }

  FW <- function(t)
  {
      
    if (r == (x$c - 1))
    {
      res <- 1 - ( (1 + cErlang * x$mu * t) * exp(-x$mu * t) )
    }
    else
    {
      aux1 <- ( r - x$c + 1 - cErlang ) * exp(-x$mu * t)
      aux2 <- cErlang * exp( (-1) * (1 - RO) * x$c * x$mu * t )
      aux <-  (aux1 + aux2)/( x$c - 1 - r )
      res <- 1 + aux
    }
    res
  }

  square_one_minus_ro <- one_minus_ro^2
  square_mu <- x$mu^2
  square_c <- x$c^2
  square_w <- W^2

  VNq <- ( RO * cErlang * (1 + RO - (RO * cErlang)) ) / square_one_minus_ro
  VTq <- ( (2 - cErlang) * cErlang) / (square_mu * square_c * square_one_minus_ro)
  VN  <- VNq + (r * (1 + cErlang)) 

  if ( r == (x$c - 1) )
    VT <- ((2 * (2 * cErlang + 1)) / square_mu) - square_w
  else
    VT <- ( ((2 * cErlang * (1 - (square_c * square_one_minus_ro))) / ((r + 1 -x$c) * square_c * square_one_minus_ro * square_mu)) + (2 / square_mu) ) - square_w

 res <- list(
    Inputs = x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, VN = VN, W = W, VT = VT, Wqq = Wqq, Lqq = Lqq, Pn = Pn, Qn = Pn, FW = FW, FWq = FWq
  )

  class(res) <- "o_MMC"
  res
}


Inputs.o_MMC     <- function(x, ...) { x$Inputs }
RO.o_MMC         <- function(x, ...) { x$RO }
Lq.o_MMC         <- function(x, ...) { x$Lq }
VNq.o_MMC        <- function(x, ...) { x$VNq }
Wq.o_MMC         <- function(x, ...) { x$Wq }
VTq.o_MMC        <- function(x, ...) { x$VTq }
L.o_MMC          <- function(x, ...) { x$L }
VN.o_MMC         <- function(x, ...) { x$VN }
W.o_MMC          <- function(x, ...) { x$W }
VT.o_MMC         <- function(x, ...) { x$VT }
Lqq.o_MMC        <- function(x, ...) { x$Lqq }
Wqq.o_MMC        <- function(x, ...) { x$Wqq } 
Pn.o_MMC         <- function(x, ...) { x$Pn }
Qn.o_MMC         <- function(x, ...) { x$Qn }
Throughput.o_MMC <- function(x, ...) { x$Throughput }


Report.o_MMC <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MMC <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMC  <- function(x, ...)
{
  print_summary(x, ...)
}
