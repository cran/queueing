###############################################################
###############################################################
## MODEL M/M/c/K/m - Finite Poblation, c servers, system     ##
## capacity lesser or equal than the poblation        		 	 ##
###############################################################
###############################################################
NewInput.MMCKM <- function(lambda=0, mu=0, c=1, k=1, m=1, method=0)
{
  res <- list(lambda = lambda, mu = mu, c = c, k = k, m = m, method = method)
  class(res) <- "i_MMCKM"
  res
}


CheckInput.i_MMCKM <- function(x, ...)
{

  MMCKM_m_warning <- "m has to be at least one"
  MMCKM_m_k <- "k must be equal or lesser than the poblation m"
  MMCKM_class <- "The class of the object x has to be M/M/c/K/m (i_MMCKM)"
  MMCKM_anomalous <- "Some value of lambda, mu, c, k or m is anomalous. Check the values."
  MMCKM_method <- "method variable has to be 0 to be exact calculus, 1 to be aproximate calculus"
  MMCKM_m_integer <- "the poblation (m) must be an integer number"

 if (class(x) != "i_MMCKM")
  stop(MMCKM_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$k) || is.anomalous(x$m)
  )
    stop(MMCKM_anomalous)

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

 if (x$m < 1)
 	 stop(MMCKM_m_warning)

 if (!is.wholenumber(x$m))
   stop(MMCKM_m_integer)

 if (x$k < x$c)
	 stop(ALL_k_c)
 
 if (x$m < x$k)
	 stop(MMCKM_m_k)

 if (x$method != 0 && x$method != 1)
   stop(MMCKM_method)

}


MMCKM_InitPn_Exact <- function(x)
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
		
		pn[1] <- totu * totfact
    fn[1] <- totfact
    sum <- pn[1]

		i <- 1
		while (i <= x$k)
		{
			totu <- totu * u
		  factn <- factn * i
			# Factorial calculus
		  if (i <= x$m/2)
			{
				totfact <- totfact * (x$m - i + 1) / i
				fn[i+1] <- totfact
      }
			else
				totfact <- fn[x$m - i + 1]
		  	
			if (i == x$c) factc <- factn

			if (i > x$c)
			{
				potc <- potc * x$c
				totaux <- factn / (factc * potc)
		    pn[i+1] <- totfact * totu * totaux
        sum <- sum + pn[i+1]
			}
			else
      {
        pn[i+1] <- totfact * totu
        sum <- sum + pn[i+1]
      }		

			i <- i + 1
		}
    pn/sum
}


MMCKM_InitPn_Aprox_AuxC <- function(n, lambda, mu, c, k, m)
{
  (lfactorial(m) - lfactorial(m-n) - lfactorial(n)) + (n * log(lambda/mu))
}


MMCKM_InitPn_Aprox_AuxK <- function(n, lambda, mu, c, k, m)
{
  toC <- MMCKM_InitPn_Aprox_AuxC(n, lambda, mu, c, k, m)
  toK <- lfactorial(n) - lfactorial(c) - (n - c) * log(c)
  toC + toK
}


MMCKM_InitPn_Aprox <- function(x)
{
  ProbFactCalculus(
    x$lambda, x$mu, x$c, x$k, x$m, x$k, MMCKM_InitPn_Aprox_AuxC, MMCKM_InitPn_Aprox_AuxK, MMCKM_InitPn_Aprox_AuxK
  )
}


MMCKM_InitPn <- function(x)
{
  if (x$method == 0)
    pn <- MMCKM_InitPn_Exact(x)
  else
    pn <- MMCKM_InitPn_Aprox(x)

  pn
}


QueueingModel.i_MMCKM <- function(x, ...)
{
 CheckInput.i_MMCKM(x, ...)
 Pn <- MMCKM_InitPn(x)

 # To control the cases where the probabilties doesn't make sense, that it is going to be saturation
 if ( (x$method == 1 && sum(Pn) == 0) || (x$method == 0 && sum(is.nan(Pn)) != 0) )
 {
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
 }
 else
 {
   i_per_pn_i <- (0:x$k) * Pn[1:(x$k+1)]
   sum_pn_0_c_minus_1 <- sum(Pn[1:x$c])

   L <- sum(i_per_pn_i)
   Lq <- L - x$c - sum(i_per_pn_i[1:x$c]) + (x$c * sum_pn_0_c_minus_1)
   
   if (x$k < x$m)
     Throughput <- x$mu * (L - Lq)
   else #x$k == x$m
     Throughput <- x$lambda * (x$m - L)

   W <- L / Throughput
   Wq <- Lq / Throughput 
   RO <- Throughput / (x$c * x$mu)

   QnAux <- function(n){ Pn[n] * (x$m - (n-1)) / ( (x$m - L) - ( (x$m - x$k) * Pn[x$k+1] ) ) }
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

   # variances
   VN  <- sum( (0:x$k)^2 * Pn ) - (L^2)

   if (x$c == x$k)
     VNq <- 0
   else
     VNq <- sum( ( c( rep(0, x$c+1), 1:(x$k-x$c) )^2 * Pn) - (Lq^2) ) 
 }


  # FW <- function(t){
  #  aux <- function(n) { Qn[n+x$c-1] * dist(n-1) }
  #  1 - sum(sapply(seq(1, x$k-x$c+1, 1), aux))
  # }

  # FWq <- function(t){
  #   aux <- function(n) { Qn[n+x$c] * dist(n-1) }
   
  #   if (x$c == x$k)
  #     0
  #   else
  #     1 - sum(sapply(seq(1, x$k-x$c, 1), aux))
  # }


  # The result
  res <- list(Inputs=x, RO = RO, Lq = Lq, VNq = VNq, Wq = Wq, Throughput = Throughput,
     L = L, VN = VN, W = W, Lqq = Lqq, Wqq = Wqq, Pn = Pn, Qn = Qn)

  class(res) <- "o_MMCKM"
  res
} 

Inputs.o_MMCKM     <- function(x, ...) { x$Inputs }
L.o_MMCKM          <- function(x, ...) { x$L }
VN.o_MMCKM         <- function(x, ...) { x$VN }
Lq.o_MMCKM         <- function(x, ...) { x$Lq }
VNq.o_MMCKM        <- function(x, ...) { x$VNq }
Lqq.o_MMCKM        <- function(x, ...) { x$Lqq }
Throughput.o_MMCKM <- function(x, ...) { x$Throughput }
W.o_MMCKM          <- function(x, ...) { x$W }
RO.o_MMCKM         <- function(x, ...) { x$RO }
Wq.o_MMCKM         <- function(x, ...) { x$Wq }
Wqq.o_MMCKM        <- function(x, ...) { x$Wqq }
Pn.o_MMCKM         <- function(x, ...) { x$Pn }
Qn.o_MMCKM         <- function(x, ...) { x$Qn }


Report.o_MMCKM    <- function(x, ...)
{ 
  reportAux(x)
}


summary.o_MMCKM <- function(object, ...)
{ 
  aux <- list(el=CompareQueueingModels(object))
  class(aux) <- "summary.o_MM1"
  aux
}


print.summary.o_MMCKM  <- function(x, ...)
{
  print_summary(x, ...)
}
