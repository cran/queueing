###############################################################
###############################################################
## MODEL M/M/c/K/K - Finite Plobation, c servers        		 ##
###############################################################
###############################################################
NewInput.MMCKK <- function(lambda=0, mu=0, c=1, k=1, method=0)
{
  res <- list(lambda = lambda, mu = mu, c = c, k = k, method = method)
  class(res) <- "i_MMCKK"
  res
}


CheckInput.i_MMCKK <- function(x, ...)
{
  MMCKK_class <- "The class of the object x has to be M/M/c/K/K (i_MMCKK)"
  MMCKK_anomalous <- "Some value of lambda, mu, c or k is anomalous. Check the values."
  MMCKK_method <- "method variable has to be 0 to be exact calculus, 1 to be aproximate calculus"


  if (class(x) != "i_MMCKK")
    stop(MMCKK_class)

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$k)
  )
    stop(MMCKK_anomalous)

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

  if (x$method != 0 && x$method != 1 && x$method != 2)
    stop(MMCKK_method)

}


MMCKK_InitPn_Aprox_AuxC <- function(n, lambda, mu, c, k, m)
{
  (lfactorial(k) - lfactorial(k-n) - lfactorial(n)) + (n * log(lambda/mu))
}


MMCKK_InitPn_Aprox_AuxK <- function(n, lambda, mu, c, k, m)
{
  toC <- MMCKK_InitPn_Aprox_AuxC(n, lambda, mu, c, k, m)
  toK <- lfactorial(n) - lfactorial(c) - (n - c) * log(c)
  toC + toK
}


MMCKK_InitPn_Aprox <- function(x)
{
  ProbFactCalculus(
    x$lambda, x$mu, x$c, x$k, x$k, x$k, MMCKK_InitPn_Aprox_AuxC, MMCKK_InitPn_Aprox_AuxK, MMCKK_InitPn_Aprox_AuxK
  )
}


MMCKK_method2_Aux <- function(x, i)
{
  r <- x$lambda/x$mu

  if (i <= x$c-1)
  {
    (x$k - i) * r / (i+1)
  }
  else
  {
    (x$k - i) * r / x$c
  }
}


MMCKK_method2_Prod <- function(x,n)
{
  prod <- 1
  
  for (i in 0:(n-1))
  {
    prod <- prod * MMCKK_method2_Aux(x, i)
  }

  prod

}


MMCKK_method2_Prob <- function(x)
{

  pn <- c()
  
  sumAux <- 1

  for (i in (1:x$k))
  {
    sumAux <- sumAux + MMCKK_method2_Prod(x, i)
  }

  pn[1] <- 1/sumAux

  for (i in 2:(x$k+1))
  {
    pn[i] <- MMCKK_method2_Aux(x, i-2) * pn[i-1]
  }

  pn
}



MMCKK_InitPn_Exact <- function(x)
{
		pn <- c(0:x$k)
		fn <- c(0:x$k)

		u <- x$lambda / x$mu
		totu <- 1
		totfact <- 1
		factn <- 1
		factc <- 0
		totaux <- 1
		potc <- 1
    sumpn <- 0
		
		pn[1] <- totu * totfact
    fn[1] <- totfact
    sumpn <- sumpn + pn[1]

		i <- 1
		while (i <= x$k)
		{
			totu <- totu * u
		  factn <- factn * i
			# Factorial calculus
		  if (i <= x$k/2)
			{
				totfact <- totfact * (x$k - i + 1) / i
				fn[i+1] <- totfact
      }
			else
				totfact <- fn[x$k - i + 1]
				
			if (i == x$c) factc <- factn

			if (i > x$c)
			{
				potc <- potc * x$c
				totaux <- factn / (factc * potc)
		    pn[i+1] <- totfact * totu * totaux
        sumpn <- sumpn + pn[i+1]
			}
			else
      {
        pn[i+1] <- totfact * totu
        sumpn <- sumpn + pn[i+1]
      }
			i <- i + 1
		}
    pn/sumpn
}


MMCKK_InitPn <- function(x)
{
  # check if c=k so the distribution is a binomial
  if (x$c == x$k)
  {
    u <- x$lambda / x$mu
    prob <- u / (u + 1)
    pn <- sapply(seq(0, x$k, 1), function(i){dbinom(i, x$k, prob)})
  }
  else
  {
    if (x$method == 0)
      pn <- MMCKK_InitPn_Exact(x)
    else
    {
      if (x$method == 1)
        pn <- MMCKK_InitPn_Aprox(x)
      else
        pn <- MMCKK_method2_Prob(x)
    }
  }
  
  pn
}


QueueingModel.i_MMCKK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMCKK(x, ...)

  Pn <- MMCKK_InitPn(x)

  # To control the cases where the probabilties doesn't make sense, that it is going to be saturation
  if ( (x$method == 1 && sum(Pn) == 0) || (x$method == 0 && sum(is.nan(Pn)) != 0) )
  {
    #print(paste("sum(Pn):", sum(Pn)))
    RO <- 1
    Throughput <- (x$c * x$mu)
    L <- (x$k - (Throughput/x$lambda))
    
    if (L <= 0)
    {
      W <- NA
      L <- NA 
      Wq <- NA   
      Lq <- NA
      Wqq <- NA
      Lqq <- NA
    }
    else
    {
      W <- L/Throughput
      Wq <- W - (1/x$mu)
      Lq <- Throughput * Wq
      Wqq <- NA
      Lqq <- NA
    }
    VN <- NA
    VNq <- NA
    VTq <- NA  
  }
  else
  {
    k_per_pk <- c(0:x$k) * Pn[1:(x$k+1)]
    sum_pn_0_c_minus_1 <- sum(Pn[1:x$c])
    L <- sum(k_per_pk)
    Lq <- L - x$c - sum(k_per_pk[1:x$c]) + (x$c * sum_pn_0_c_minus_1)
    Throughput <- x$lambda * (x$k - L)
    W <- L / Throughput
    RO <-  Throughput / (x$c * x$mu)
    Wq <- Lq / Throughput

    QnAux <- function(n){ Pn[n] * (x$k - (n-1)) / (x$k - L) }
    Qn <- sapply(1:x$k, QnAux)

    if (x$k == x$c)
    {
      Wqq <- NA
      Lqq <- NA
    }
    else
    {
      Wqq <- Wq / (1 - sum(Qn[1:x$c]))
      Lqq <- Wqq * x$c * x$mu   
    }

    if (x$c == x$k)
      FWq <- function(t){0}
    else
    {
      FWq <- function(t){
        aux <- function(n) { Qn[n+x$c] * ppois(n-1, x$c * x$mu * t) }
        1 - sum(sapply(seq(1, x$k-x$c, 1), aux))
      }
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

    #Wqq <- Wq / (1-sum_pn_0_c_minus_1) 
  }    
  
  # dist <- function(n) { ppois(n, x$c * x$mu * t) }
  
  # FW <- function(t){    
  #   aux <- function(n) { Qn[n+x$c-1] * dist(n-1) }
  #   1 - sum(sapply(seq(1, x$k-x$c+1, 1), aux))
  # }


  # The result
  res <- list(
    Inputs=x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, VTq = VTq, Throughput = Throughput,
    L = L, VN = VN, W = W, Lqq = Lqq, Wqq = Wqq,
    Pn = Pn, Qn = Qn, FWq = FWq
  )
  
  class(res) <- "o_MMCKK"
  res

} 

Inputs.o_MMCKK     <- function(x, ...) { x$Inputs }
L.o_MMCKK          <- function(x, ...) { x$L }
VN.o_MMCKK         <- function(x, ...) { x$VN }
Lq.o_MMCKK         <- function(x, ...) { x$Lq }
VNq.o_MMCKK        <- function(x, ...) { x$VNq }
Lqq.o_MMCKK        <- function(x, ...) { x$Lqq }
Throughput.o_MMCKK <- function(x, ...) { x$Throughput }
W.o_MMCKK          <- function(x, ...) { x$W }
RO.o_MMCKK         <- function(x, ...) { x$RO }
Wq.o_MMCKK         <- function(x, ...) { x$Wq }
VTq.o_MMCKK        <- function(x, ...) { x$VTq }
Wqq.o_MMCKK        <- function(x, ...) { x$Wqq }
Pn.o_MMCKK         <- function(x, ...) { x$Pn }
Qn.o_MMCKK         <- function(x, ...) { x$Qn }

Report.o_MMCKK <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MMCKK <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMCKK  <- function(x, ...)
{
  print_summary(x, ...)
}
