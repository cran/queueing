############################################################
## Auxiliary functions
############################################################
is.anomalous <- function(x)
{
  # is.nan(x) doesn't work for lists
  #is.null(x) || is.na(x) || is.nan(x)
  is.null(x) || is.na(x)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


print_summary  <- function(x, ...)
{
  print(x$el, ...)
  invisible(x)
}


C_erlang2 <- function(c, r)
{
  ro <- r / c

  totr <- 1
  totn <- 1
  total <- totr / totn

	i <- 1
	while (i <= c-1)
  {
		totr <- totr * r
		totn <- totn * i
		total <- total + (totr / totn)
		i <- i + 1
  }

  totr <- totr * r
  totn <- totn * c
  numerator <- totr / totn
  denominator <- (1 - ro) * (total + (numerator / (1 - ro)))
  numerator / denominator  
}

nodes <- function(...)
{
  list(...)
}


checkNegative <- function(v)
{
  return(sum(v < 0) > 0)
}


checkNegativeOrZero <- function(v)
{
  return(sum(v <= 0) > 0)
}


checkAtLeastOne <- function(v)
{
  return(sum(v < 1) > 0)
}


checkAllZero <- function(v)
{
  return(sum(v <= 0) == length(v))
}


C_erlang3 <- function(c, r)
{
  b_result <- B_erlang(c, r)
  num <- c * b_result
  den <- c - (r * (1 - b_result))
  num / den    
}


# this saves one step of B_erlang, more efficient
C_erlang <- function(c=1, r=0)
{
  if (is.anomalous(c))
    stop("The parameter c is anomalous. Check it!")

  if (is.anomalous(r))
    stop("The parameter r is anomalous. Check it!")

  if (c<1)
    stop(ALL_c_warning)

  if (!is.wholenumber(c))
    stop(ALL_c_integer)

  b_result <- B_erlang(c-1, r)
  ( r * b_result ) / ( c - (r * (1 - b_result)) )   
}


# recursive version, in R has problems of stack overflow when c is large
B_erlang2 <- function(c, u)
{

	f <- function(c, u)
	{
		if (c == 0) 1
		else ( 1 + ( f(c-1, u) * (c/u) ) )
	}
	1 / f(c, u)
}


# definition version
B_erlang3 <- function(c, u)
{
  n_fact <- 1
  u_power <- 1
  tot <- u_power / n_fact

  i <- 1
  while (i <= c)
  {
	n_fact <- i * n_fact
    u_power <- u_power * u
    tot <- tot + (u_power / n_fact)
    i <- i + 1
  }
  
  (u_power / n_fact) / tot
	
}


B_erlang4 <- function(c=1, u=0)
{

  if (is.anomalous(c))
    stop("The parameter c is anomalous. Check it!")

  if (is.anomalous(u))
    stop("The parameter u is anomalous. Check it!")

  if (c<0)
    stop("The number of servers can not be less than zero!")

  tot <- 1
  aux <- 1 / u
  i <- 1

  while (i <= c)
  {
    tot <- 1 + (tot * aux)
    aux <- aux + (1 / u)
    i <- i + 1
  }

  1/tot
	
}

B_erlang5 <- function(c=1, u=0)
{

  if (is.anomalous(c))
    stop("The parameter c is anomalous. Check it!")

  if (is.anomalous(u))
    stop("The parameter u is anomalous. Check it!")

  if (c<0)
    stop("The number of servers can not be less than zero!")

  if (u<0)
    stop("The u parameter can not be negative!")

  if (!is.wholenumber(c))
    stop("The parameter c has to be an integer number")

  1/Reduce(function(i, j) 1 + (i * j), (1 / u) * seq(1, c), init=1)    
}


B_erlang_6 <- function(c=1, u=0)
{

  if (is.anomalous(c))
    stop("The parameter c is anomalous. Check it!")

  if (is.anomalous(u))
    stop("The parameter u is anomalous. Check it!")

  if (c<0)
    stop("The number of servers can not be less than zero!")

  if (u<0)
    stop("The u parameter can not be negative!")

  if (!is.wholenumber(c))
    stop("The parameter c has to be an integer number")

  if (c==0) return(1)
  if (u==0) return(0)

  val <- 1/u
  tot <- 1
  
  for (i in 1:c)
    tot <- 1 + (tot * i * val)

  1/tot    
}


B_erlang <- function(c=1, u=0)
{

  if (is.anomalous(c))
    stop("The parameter c is anomalous. Check it!")

  if (is.anomalous(u))
    stop("The parameter u is anomalous. Check it!")

  if (c<0)
    stop("The number of servers can not be less than zero!")

  if (u<0)
    stop("The u parameter can not be negative!")

  if (!is.wholenumber(c))
    stop("The parameter c has to be an integer number")

  if (c==0) return(1)
  if (u==0) return(0)

  exp(
    dgamma(u, shape=c+1, scale=1, log=TRUE) -
      pgamma(u, shape=c+1, scale=1, log.p=TRUE, lower.tail=FALSE)
  )    
}


Engset <- function(k=1, c=0, r=0)
{
  if (is.anomalous(c)) 
    stop("The parameter c is anomalous. Check it!")
  if (is.anomalous(r)) 
    stop("The parameter r is anomalous. Check it!")
  if (is.anomalous(k))
    stop("The parameter k is anomalous. Check it!")
  if (c < 0) 
    stop("The number of servers can not be less than zero!")
  if (r < 0) 
    stop("The r parameter can not be negative!")
  if (!is.wholenumber(c)) 
    stop("The parameter c has to be an integer number")
  if (!is.wholenumber(k)) 
    stop("The parameter k has to be an integer number")

  if (c > k)
    stop("c can not be greater than k")

  if (c == k)
    return(0)
  if (c == 0) 
    return(1)
  if (r == 0) 
    return(0)


  acum <- 1 
  for (i in (1:c)) acum <- ( (i * acum) / ((k - i + 1) * r) ) + 1
  return(1/acum)
}


Engset_def <- function(k=1, c=0, r=0)
{
  num <- choose(k-1, c) * (r^c)
  
  acum <- 0
  
  for (i in (0:c))
    acum <- acum + (choose(k-1, i) * (r^i))
  
  den <- acum
  num/den
}

engset_bin <- function(x, k=1, c=1, r=0)
{
  a <- r / (1 + r)
  num <- dbinom(x, k, a)
  den <- pbinom(c, k, a)
  num/den
}




ProbFactCalculus <- function(lambda, mu, c, k, m, limit, fAuxC, fAuxK, fAuxM)
{
  pn <- c(0:limit)

  pn[1] <- 0
  
  i <- 1
  while (i <= limit)
  {
    if (i <= c)
    {
      pn[i+1] <- fAuxC(i, lambda, mu, c, k, m)
      #print(paste(paste(paste("pn[", i), "]: "), pn[i+1]))      
    }
    else
    {
      if (i <= k)
      {
        pn[i+1] <- fAuxK(i, lambda, mu, c, k, m)
        #print(paste(paste(paste("pn[", i), "]: "), pn[i+1]))      
      }  
      else
      {
        pn[i+1] <- fAuxM(i, lambda, mu, c, k, m)
        #print(paste("pn[i+1]: ", pn[i+1]))
      }
    }
    i <- i + 1
  }

  p0 <- -log(sum(exp(pn)))

  pn <- exp(pn + p0)
  pn
}


tpoisson <- function(n, maximum, lambda)
{
  dpois(n, lambda)/ppois(maximum, lambda)
}


############################################################
############################################################
## FUNTION TO REPORT BASIC MARKOVIAN MODELS 
############################################################
############################################################

reportAux <- function(object)
{ 
  oclass <- class(object)

  Ls <- object$L - object$Lq

  if (oclass == "o_MM1")
  {
    cat("The inputs of the M/M/1 model are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", n: ", object$Inputs$n, "\n", sep=""))
    cat("\n")
    cat("The outputs of the M/M/1 model are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pn) of the n = ", object$Inputs$n, " clients in the system are:\n", sep=""))
    cat(object$Pn)
    cat("\n")
  }  
  else if (oclass == "o_MMCKM")
  {
    cat("The inputs of the model M/M/c/K/m are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", c: ", object$Inputs$c, ", k: ", object$Inputs$k, " ,m: ", object$Inputs$m, ", method: ", object$Inputs$method, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/c/K/m are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
    cat(paste("The mean think time is : ", 1/object$Inputs$lambda, "\n", sep=""))
  }
  else if (oclass == "o_MMC")
  {
    method <- if (object$Inputs$method == 0) "Exact" else "Aprox"

    cat("The inputs of the model M/M/c are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", c: ", object$Inputs$c, ", n: ", object$Inputs$n, ", method: ", method, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/c are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pn) of the n = ", object$Inputs$n, " clients in the system are:\n", sep=""))
    cat(object$Pn)
    cat("\n")
  }
  else if (oclass == "o_MM1KK")
  {
    cat("The inputs of the model M/M/1/K/K are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", k: ", object$Inputs$k, ", method: ", object$Inputs$method, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/1/K/K are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
    cat(paste("The mean think time is : ", 1/object$Inputs$lambda, "\n", sep=""))
  }
  else if (oclass == "o_MMCKK")
  {
    cat("The inputs of the model M/M/c/K/K are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", c: ", object$Inputs$c, ", k: ", object$Inputs$k, ", method: ", object$Inputs$method, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/c/K/K are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
    cat(paste("The mean think time is : ", 1/object$Inputs$lambda, "\n", sep=""))
  }
  else if (oclass == "o_MMCKM")
  {
    cat("The inputs of the model M/M/c/K/m are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", c: ", object$Inputs$c, ", k: ", object$Inputs$k, " ,m: ", object$Inputs$m, ", method: ", object$Inputs$method, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/c/K/m are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
    cat(paste("The mean think time is : ", 1/object$Inputs$lambda, "\n", sep=""))
  }
  else if (oclass == "o_MMInfKK")
  {
    cat("The inputs of the model M/M/Inf/K/K are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", k: ", object$Inputs$k, ", method: ", object$Inputs$method, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/Inf/K/K are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
    cat(paste("The mean think time is : ", 1/object$Inputs$lambda, "\n", sep=""))
  }
  else if (oclass == "o_MMInf")
  {
    cat("The inputs of the model M/M/Infinite are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", n: ", object$Inputs$n, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/Infinite are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pn) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")  
  }
  else if (oclass == "o_MM1K")
  {
    cat("The inputs of the model M/M/1/K are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", k: ", object$Inputs$k, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/1/K are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
    cat(paste("The probability (q0, q1, ..., qk-1) that a client that enters meets n clients in the system are:\n"))
    cat(object$Qn)
    cat("\n")
  }
  else if (oclass == "o_MMCK")
  {
    cat("The inputs of the model M/M/c/K are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", c: ", object$Inputs$c, ", k: ", object$Inputs$k, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/c/K are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pk) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
  }
  else if (oclass == "o_MMCC")
  { 
    cat("The inputs of the model M/M/c/c are:\n")
    cat(paste("lambda: ", object$Inputs$lambda, ", mu: ", object$Inputs$mu, ", c: ", object$Inputs$c, "\n", sep=""))
    cat("\n")
    cat("The outputs of the model M/M/c/c are:\n")
    cat("\n")
    cat(paste("The probability (p0, p1, ..., pc) of the clients in the system are:\n"))
    cat(object$Pn)
    cat("\n")
  }
 
  cat(paste("The traffic intensity is: ", Ls, "\n", sep=""))
  cat(paste("The server use is: ", object$RO, "\n", sep=""))
  cat(paste("The mean number of clients in the system is: ", object$L, "\n", sep=""))
  cat(paste("The mean number of clients in the queue is: ", object$Lq, "\n", sep=""))
  cat(paste("The mean number of clients in the server is: ", Ls, "\n", sep=""))
  cat(paste("The mean time spend in the system is: ", object$W, "\n", sep=""))
  cat(paste("The mean time spend in the queue is: ", object$Wq, "\n", sep=""))
  cat(paste("The mean time spend in the server is: ", object$W - object$Wq, "\n", sep=""))
  cat(paste("The mean time spend in the queue when there is queue is: ", object$Wqq, "\n", sep=""))
  cat(paste("The throughput is: ", object$Throughput, "\n", sep=""))

  if (oclass == "o_MM1KK")
  {
    cat(paste("The normalized average response time is: ", object$WWs, "\n", sep=""))
    cat(paste("The saturation point is: ", object$SP, "\n", sep=""))
  }  
}


############################################################
############################################################
## FUNTION TO REPORT SINGLE CLASS NETWORKS 
############################################################
############################################################

reportSingleClass <- function(object)
{
   
  classObject <- class(object)

  if (classObject == "o_OJN")
    cat("The inputs of the open Jackson network are:\n\n")
  else # has to be o_CJN
    cat("The inputs of the closed Jackson network are:\n\n")
  
  print(object$Inputs)
  cat("\n\n")

  if (classObject == "o_OJN")
    cat(paste("The outputs of the open Jackson network are:", "\n\n", sep=""))
  else # has to be o_CJN
    cat("The outputs of the closed Jackson network are:\n\n")
  
  cat("---------- Complete network -------------------------\n\n")
  cat(paste("The mean number of clients in the network is: ", object$L, "\n", sep=""))
  cat(paste("The mean time spend in the network is: ", object$W, "\n", sep=""))
  cat(paste("The throughput of the network is: ", object$Throughput, "\n", sep=""))
  cat("\n\n")


  cat("--------- Per node ---------------------------------\n\n")
 
  i <- 1
  while (i <= length(object$ROk))
  {
    cat(paste("The use of node ", i, " is: ", object$ROk[i], "\n", sep=""))
    cat(paste("The throughput of node ", i, " is: ", object$Throughputk[i], "\n", sep=""))
    cat(paste("The mean number of clients in node ", i, " is: ", object$Lk[i], "\n", sep=""))
    cat(paste("The mean time spend in node ", i, " is: ", object$Wk[i], "\n", sep=""))

    if (classObject == "o_OJN")
    {
      cat(paste("The probability (p0, p1, ..., pn) or visit ratio of node ", i, " is: ", "\n", sep=""))
      print(object$Pn[[i]])  
    }
   
    cat("\n\n")
    i <- i + 1
  }
}


############################################################
############################################################
## FUNTION TO REPORT MULTIPLE CLASS NETWORKS 
############################################################
############################################################

reportMultiClass <- function(object)
{
  if (class(object) != "o_MCON" && class(object) != "o_MCCN" && class(object) != "o_MCMN")
    stop("Incorrect class")

  if (class(object) == "o_MCON")
    netType <- "open"
  else if (class(object) == "o_MCCN")
    netType <- "closed"
  else
    netType <- "mixed"

  cat(paste("The inputs of the multiclass ", netType, " network are:", "\n", sep=""))
  print(object$Inputs)
  cat("\n\n")
  cat(paste("The outputs of the multiclass ", netType, " network are:", sep=""))
  cat("\n\n")

  cat("---------- Complete network -------------------------\n\n")
  cat(paste("The mean number of clients in the network is: ", object$L, "\n", sep=""))
  cat(paste("The mean time spend in the network is: ", object$W, "\n", sep=""))
  cat(paste("The throughput of the network is: ", object$Throughput, "\n", sep=""))
  cat("\n\n")
  
  cat("---------- Per Class -------------------------\n\n")

  for (i in (1:object$Inputs$classes))
  {
    cat(paste("The mean number of class ", i, " clients in the network is: ", object$Lc[i], "\n", sep=""))
    cat(paste("The mean time spend in the network per class ", i ," is: ", object$Wc[i], "\n", sep=""))
    cat(paste("The throughput of class " , i ," of the network is: ", object$Throughputc[i], "\n", sep=""))
    cat("\n\n")    
  }
  
  cat("--------- Per node ---------------------------------\n\n")
  
  for (i in (1:object$Inputs$nodes))
  {
    cat(paste("The use of node ", i, " is: ", object$ROk[i], "\n", sep=""))
    cat(paste("The mean number of clients in node ", i, " is: ", object$Lk[i], "\n", sep=""))
    cat(paste("The mean time spend in node ", i, " is: ", object$Wk[i], "\n", sep=""))
    cat(paste("The throughput of node ", i, " is: ", object$Throughputk[i], "\n", sep=""))
    cat("\n\n")
  }

  cat("--------- Per class and node -----------------------\n\n")

  for (i in (1:object$Inputs$classes))
  {
    for (j in (1:object$Inputs$nodes))
    {
      cat(paste("The class ", i, " use of node ", j, " is: ", object$ROck[i, j], "\n", sep=""))
      cat(paste("The mean number of class ", i, " clients in node ", j, " is: ", object$Lck[i, j], "\n", sep=""))
      cat(paste("The mean time spend by class ", i, " in node ", j, " is: ", object$Wck[i, j], "\n", sep=""))
      cat(paste("The throughput of class " , i, " in node ", j, " is: ", object$Throughputck[i, j], "\n", sep=""))
      cat("\n\n")
    }
  }  
}



############################################################
############################################################
## FUNTION TO TABULATE DIFFERENT MODELS 
############################################################
############################################################

CompareQueueingModels2 <- function(models)
{
  num_elems <- length(models)

  # Check that every object has the correct class
  for (i in 1:num_elems)
  {
    classm <- class(models[[i]])
    if (!(classm == "o_MM1" || classm == "o_MMC" || classm == "o_MM1K" || classm == "o_MMCK" ||
      classm == "o_MMCC" || classm == "o_MMInf" || classm == "o_MMInfKK" || classm == "o_MM1KK" ||
      classm == "o_MMCKK" || classm == "o_MMCKM"))
    stop("Function called with incorrect class")
  }

  #build the table
  lambda      <- c()
  mu          <- c()
  c           <- c()
  k           <- c()
  m           <- c()  
  RO          <- c()
  Lq          <- c()
  Wq          <- c()
  Throughput  <- c()
  L           <- c()
  W           <- c()
  Wqq         <- c()
  Lqq         <- c()
  P0          <- c()

  for (i in (1:num_elems))
  {
    mod <- models[[i]]
    classm <- class(mod)   

    lambda <- c(lambda, mod$Inputs$lambda)
    mu <- c(mu, mod$Inputs$mu)

    if (classm == "o_MMInf" || classm == "o_MMInfKK")
      c <- c(c, NA)
    else if (classm == "o_MM1" || classm == "o_MM1K" || classm == "o_MM1KK")
      c <- c(c, 1)
    else
      c <- c(c, mod$Inputs$c)

    if (classm == "o_MMInf" || classm == "o_MM1" || classm == "o_MMC")
      k <- c(k, NA)
    else if (classm == "o_MMCC")
      k <- c(k, mod$Inputs$c)
    else
      k <- c(k, mod$Inputs$k)

    if (classm != "o_MMCKM")
      m <- c(m, NA)
    else
      m <- c(m, mod$Inputs$m)

    if (classm == "o_MMInf" || classm == "o_MM1" || classm == "o_MMC")
    {
      if (mod$Inputs$n >= 0)
        P0 <- c(P0, mod$Pn[1])
      else
        P0 <- c(P0, NA)
    }
    else
      P0 <- c(P0, mod$Pn[1])
   

    RO          <- c(RO, mod$RO)
    Lq          <- c(Lq, mod$Lq)
    Wq          <- c(Wq, mod$Wq)
    Throughput  <- c(Throughput, mod$Throughput)
    L           <- c(L, mod$L)
    W           <- c(W, mod$W)
    Lqq         <- c(Lqq, mod$Lqq)
    Wqq         <- c(Wqq, mod$Wqq)
  }
  
  data.frame(
    lambda=lambda, mu=mu, c=c, k=k, m=m, RO=RO, P0 = P0, Lq=Lq, Wq=Wq, X=Throughput,
    L=L, W=W, Wqq=Wqq, Lqq=Lqq
  )
}


CompareQueueingModels <- function(model, ...)
{
  models <- c(list(model), list(...))
  CompareQueueingModels2(models)    
}



############################################################
############################################################
## FUNTION TO SUMMARIZE SINGLE CLASS NETWORKS 
############################################################
############################################################

summarySingleClass <- function(object)
{

  classObject <- class(object)

  if (class(object) != "o_OJN" && class(object) != "o_CJN")
    stop("Incorrect class")

  #build the table

  # complete network
  Throughput  <- c()
  L           <- c()
  W           <- c()

  # Node
  ROk          <- c()
  Throughputk  <- c()
  Lk           <- c()
  Wk           <- c()

  rnames       <- c()

  # Complete Network
  L <- c(L, object$L)  
  W <- c(W, object$W)
  Throughput <- c(Throughput, object$Throughput)

  ROk          <- c(ROk, NA)
  Throughputk  <- c(Throughputk, NA)
  Lk           <- c(Lk, NA)
  Wk           <- c(Wk, NA)

  rnames       <- c(rnames, "Net")   
 
  i <- 1
  while (i <= length(object$ROk))
  {
    L            <- c(L, NA)  
    W            <- c(W, NA)
    Throughput   <- c(Throughput, NA)
    ROk          <- c(ROk, object$ROk[i])
    Throughputk  <- c(Throughputk, object$Throughputk[i])
    Lk           <- c(Lk, object$Lk[i])
    Wk           <- c(Wk, object$Wk[i])

    rnames       <- c(rnames, paste("Nd", as.character(i), sep=""))

    i <- i + 1
  }


  res <- 
    data.frame(
      L  = L,  W  = W,  X  = Throughput,
      Lk = Lk, Wk = Wk, Xk = Throughputk, ROk = ROk 
    )

  # Rowname change
  rownames(res) <- rnames
  res

}


############################################################
############################################################
## FUNTION TO SUMMARIZE MULTIPLE CLASS NETWORKS 
############################################################
############################################################


summaryMultiClass <- function(object)
{
  if (class(object) != "o_MCON" && class(object) != "o_MCCN" && class(object) != "o_MCMN")
    stop("Incorrect class")


  #build the table

  # complete network
  L            <- c()
  W            <- c()
  Throughput   <- c()
  

  # Class
  Lc           <- c()
  Wc           <- c()
  Throughputc  <- c()

  # Node
  ROk          <- c()
  Lk           <- c()
  Wk           <- c()
  Throughputk  <- c()
  
  # ClassNode
  ROck         <- c()
  Lck          <- c()
  Wck          <- c()
  Throughputck <- c()

  rnames       <- c() 


  # Complete Network
  L <- c(L, object$L)  
  W <- c(W, object$W)
  Throughput <- c(Throughput, object$Throughput)

  Lc          <- c(Lc, NA)
  Wc          <- c(Wc, NA)
  Throughputc <- c(Throughputc, NA)

  ROk          <- c(ROk, NA)
  Lk           <- c(Lk, NA)
  Wk           <- c(Wk, NA)
  Throughputk  <- c(Throughputk, NA)
  
  ROck         <- c(ROck, NA)
  Lck          <- c(Lck, NA)
  Wck          <- c(Wck, NA)
  Throughputck <- c(Throughputck, NA)

  rnames       <- c(rnames, "Net")  

 
  # Classes
  for (i in (1:object$Inputs$classes))
  {
    L            <- c(L, NA)
    W            <- c(W, NA)
    Throughput   <- c(Throughput, NA)
    
    Lc           <- c(Lc, object$Lc[i])
    Wc           <- c(Wc, object$Wc[i])
    Throughputc  <- c(Throughputc, object$Throughputc[i])
    
    ROk          <- c(ROk, NA)
    Lk           <- c(Lk, NA)
    Wk           <- c(Wk, NA)
    Throughputk  <- c(Throughputk, NA)

    ROck         <- c(ROck, NA)
    Lck          <- c(Lck, NA)
    Wck          <- c(Wck, NA)
    Throughputck <- c(Throughputck, NA)

    rnames       <- c(rnames, paste("Cl", as.character(i), sep=""))

  }

  # Nodes
  for (i in (1:object$Inputs$nodes))
  {
    L            <- c(L, NA)
    W            <- c(W, NA)
    Throughput   <- c(Throughput, NA)
    
    Lc           <- c(Lc, NA)
    Wc           <- c(Wc, NA)
    Throughputc  <- c(Throughputc, NA)
    
    ROk          <- c(ROk, object$ROk[i])
    Lk           <- c(Lk, object$Lk[i])
    Wk           <- c(Wk, object$Wk[i])
    Throughputk  <- c(Throughputk, object$Throughputk[i])

    ROck         <- c(ROck, NA)
    Lck          <- c(Lck, NA)
    Wck          <- c(Wck, NA)
    Throughputck <- c(Throughputck, NA)

    rnames       <- c(rnames, paste("Nd", as.character(i), sep=""))

  }

  #Per class and node
  for (i in (1:object$Inputs$classes))
  {
    for (j in (1:object$Inputs$nodes))
    {
      
      L            <- c(L, NA)
      W            <- c(W, NA)
      Throughput   <- c(Throughput, NA)
    
      Lc           <- c(Lc, NA)
      Wc           <- c(Wc, NA)
      Throughputc  <- c(Throughputc, NA)
    
      ROk          <- c(ROk, NA)      
      Lk           <- c(Lk, NA)
      Wk           <- c(Wk, NA)
      Throughputk  <- c(Throughputk, NA)

      ROck         <- c(ROck, object$ROck[i, j])
      Lck          <- c(Lck, object$Lck[i, j])
      Wck          <- c(Wck, object$Wck[i, j])
      Throughputck <- c(Throughputck, object$Throughputck[i, j])

      rnames       <- c(rnames, paste("CN", as.character(i), as.character(j), sep=""))
    }
  }


  res <-
    data.frame(
      L   = L,    W   = W,    X   = Throughput,
      Lc  = Lc,   Wc  = Wc,   Xc  = Throughputc,
      Lk  = Lk,   Wk  = Wk,   Xk  = Throughputk,   ROk  = ROk,
      Lck = Lck,  Wck = Wck,  Xck = Throughputck,  ROck = ROck
    )
  
  # Rowname change
  rownames(res) <- rnames
  res  
  
}
