############################################################
############################################################
## GENERAL BIRTH AND DEATH PROCESS
############################################################
############################################################
NewInput.BnD <- function(lambda=NULL, mu=NULL)
{
  res <- list(lambda = lambda, mu = mu)
  class(res) <- "i_BnD"
  res
}

CheckInput.i_BnD <- function(x, ...)
{
  BnD_class <- "the class of the object x has to be Birth and Death (i_BnD)"
  BnD_anomalous <- "Some value of lambda or mu is anomalous. Check the values."
  BnD_length_differ <- "The lengths of vectors lambda and mu differs. Check that the same number of elements are included"

 if (!inherits(x, "i_BnD"))
  stop(BnD_class)
  
 nlambda <- length(x$lambda)
 nmu     <- length(x$mu)

 if (nlambda != nmu)
   stop(BnD_length_differ)
 
 for (i in 1:nlambda)
 {
   if (is.anomalous(x$lambda[i]) || is.anomalous(x$mu[i]))
     stop(BnD_anomalous)
   
   if (x$mu[i] <= 0)
     stop(ALL_mu_positive)
   
   if (x$lambda[i] <= 0)
     stop(ALL_lambda_positive)
 }
}


doModelBnD <- function(lambda, mu)
{
  n <- length(lambda)
  
  prob <- rep(0, n+1)
  prob[1] <- 1
  
  for (i in 2:(n+1))
    prob[i] <- prob[i-1] * (lambda[i-1] / mu[i-1])

  
  p0 <- 1/sum(prob)
  prob <- prob * p0
  
  return(prob)
}

QueueingModel.i_BnD <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_BnD(x, ...)
  
  Pn <- doModelBnD(x$lambda, x$mu)
  len_Pn <- length(Pn)
  L <- sum((0:(len_Pn - 1)) * Pn)

  res <- list(
    Inputs = x, Pn = Pn, L = L 
  )

  class(res) <- "o_BnD"
  res
} 

Pn.o_BnD         <- function(x, ...) { x$Pn }
L.o_BnD          <- function(x, ...) { x$L }
Inputs.o_BnD     <- function(x, ...) { x$Inputs }

Report.o_BnD <- function(x, ...)
{ 
  reportBnD(x)
}

summary.o_BnD <- function(object, ...)
{ 
  aux <- list(el=summaryBnD(object))
  class(aux) <- "summary.o_BnD"
  aux
}


print.summary.o_BnD  <- function(x, ...)
{
  print_summary(x, ...)
}
