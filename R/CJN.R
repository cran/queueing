#######################################################################################
## Closed Jackson Network
#######################################################################################

NewInput2.CJN <- function(prob=NULL, n=0, z=0, operational=FALSE, method=0, tol=0.001, nodes)
{
  nds <- list(prob=prob, n=n, z=z, operational=operational, method=method, tol=tol, nodes=nodes)
  class(nds) <- "i_CJN"
  nds
}

NewInput.CJN <- function(prob=NULL, n=0, z=0, operational=FALSE, method=0, tol=0.001, ...)
{
  NewInput2.CJN(prob=prob, n=n, z=z, operational=operational, method=method, tol=tol, nodes=nodes(...))
}


NewInput3.CJN <- function(n, z, numNodes, vType, vVisit, vService, vChannel, method=0, tol=0.001)
{
  prob        <- vVisit
  operational <- TRUE
  method      <- method
  tol         <- tol

  nodes <- list()

  # Build each node
  for (i in 1:numNodes)
  {
    if (vType[i] == "Q")
    {
      if (vChannel[i] > 1)
        nodes <- c(nodes, list(NewInput.MMC(0, 1/vService[i], vChannel[i])))
      else
        nodes <- c(nodes, list(NewInput.MM1(0, 1/vService[i])))
    }
    else
      nodes <- c(nodes, list(NewInput.MMInf(0, 1/vService[i])))
  }
  
  NewInput2.CJN(prob, n, z, operational, method, tol, nodes)
}


CheckInput.i_CJN <- function(x, ...)
{
 x_class_CJN <- "x has to be of class i_CJN (Closed Jackson Network)" 
 x_anomalous <- "x has some anomalous value. Check the value(s)."
 row_distinct_col <- "x$prob (matrix class) has the number of rows distinct of the number of columns."
 row_distinct_nodes <- "x$prob (matrix class) has distinct number of rows that the number of nodes in x$nodes."
 visit_ratios_wrong <- "x$prob contains a different number of visit ratios than x$nodes."
 n_greater_zero <- "n has to be greater than zero"
 CJN_operational_logical <- "The x$operational parameter has to be of class logical (TRUE or FALSE)"
 CJN_method_values <- "The x$method has to be 0 (exact) or 1 (aprox)"
 CJN_tol_value <- "The x$tol has to be positive"

 if (
   is.anomalous(x$prob) || is.anomalous(x$nodes) || is.anomalous(x$n) ||
   is.anomalous(x$z) || is.anomalous(x$operational) || is.anomalous(x$method) ||
   is.anomalous(x$tol)
 )
    stop(x_anomalous)

 if (class(x) != "i_CJN")
   stop(x_class_CJN) 

 if (x$n <= 0)
   stop(n_greater_zero)

 num_nodes <- length(x$nodes)

 if (x$method != 0 && x$method != 1)
   stop(CJN_method_values)

 if (!(x$tol > 0))
   stop(CJN_tol_value)

 is_prob_a_matrix <- (class(x$prob) == "matrix")

 if (is_prob_a_matrix)
 {
   if (nrow(x$prob) != ncol(x$prob))
     stop(row_distinct_col)

   if (nrow(x$prob) != num_nodes)
     stop(row_distinct_nodes)
 }
 else
 {
   if (length(x$prob) != num_nodes)
     stop(visit_ratios_wrong)
 }

 if (class(x$operational) != "logical")
   stop(CJN_operational_logical) 

 i <- 1
 while (i <= num_nodes)
 {
   n = x$nodes[[i]]

   if (x$method == 0)
   {
     if (class(n) != "i_MM1" && class(n) != "i_MMC" && class(n) != "i_MMInf")
       stop(paste(paste("Node ", i), "is not of class i_MM1 or i_MMC or i_MMInf!!"))
   }
   else
   {
     if (class(n) != "i_MM1" && class(n) != "i_MMInf")
       stop(paste(paste("Node ", i), "is not of class i_MM1 or i_MMInf!!"))
   }

   CheckInput(n)
     
   i <- i + 1
 }
}


QueueingModel.i_CJN <- function(x, ...)
{
  CheckInput(x)

  if (x$method == 0)
    QueueingModelExact(x, ...)
  else
    QueueingModelApprox(x, ...)  
}


QueueingModelExact <- function(x, ...)
{
  num_nodes <- length(x$nodes)

  Throughputn <- rep(0, x$n)

  if (class(x$prob) == "matrix")
  {
    ident <- diag(dim(x$prob)[1])
    const <- matrix(data=1, nrow=dim(x$prob)[1], ncol=1)
    all1 <- matrix(data=1, ncol=dim(x$prob)[2], nrow=dim(x$prob)[1])
    prob_est <- t(solve(t(x$prob + all1 - ident), const))
  }
  else
  {
    if (x$operational) #Visit ratios as counts of repetitions has been given
    {
      prob_est <- rep(1, num_nodes)
      #We have to "correct" the mu values
      for (i in (1:num_nodes))
        x$nodes[[i]]$mu <- x$nodes[[i]]$mu / x$prob[i]     
    }
    else
      prob_est <- x$prob
  }
  
    
  #print(paste("prob_est: ", prob_est))
  
  # create the list to hold the prob
  mclass <- list()

  k <- 1
  while (k <= num_nodes)
  {
    if (class(x$nodes[[k]]) == "i_MMC" && x$nodes[[k]]$c > 1)
      mclass <- c(mclass, list(array(0, dim=c(x$nodes[[k]]$c, x$n))))
    k <- k + 1
  }

  # array initialization
  Wk <- numeric()
  Lk <- rep(times=num_nodes, 0)

  alfa <- function(i, c)
  {
    if (i <= c)
      res <- i
    else
      res <- c

    res
  }

  CalcProb <- function
    (n, c, mu, prob_est_mmcnode, thro, probC)
  {
    if (n == 1)
    {
      probC[, n] <- 0
      probC[1, n] <- 1
    }
    else #n>=2
    {
      sum <- 0
      j <- c
      while (j>=2)
      {
        #print(paste("j: ", j))
        #print(paste("n: ", n))
        #print(paste("probC[j-1, n-1]:", probC[j-1, n-1]))
        #print(paste("probC[j, n]:", probC[j, n]))
        probC[j, n] <- 
          ( (prob_est_mmcnode * thro) / (mu * alfa(j-1, c)) ) * probC[j-1, n-1]
        #print(paste("iterating inside, prob cond de ", j, " | ", n, " es: ", probC[j, n]))
        sum <- sum + ( (c - (j-1)) * probC[j, n] ) 
        j <- j - 1
      }

      sum <- (1/c) * (sum + (prob_est_mmcnode * thro / mu)) 
      probC[1, n] <- 1 - sum
      #print(paste("inside, prob cond de 1", " | ", n, " es: ", probC[1, n]))
    }
    probC
  }

  #print("Se compila la funcion, se va a entrar en el bucle")

  i <- 1
  while (i <= x$n)
  {  
    # change before putting the x$z  
    #tmp <- 0
    tmp <- x$z
    k <- 1
    num_mmc <- 0
    while (k <= num_nodes)
    {
      if (class(x$nodes[[k]]) == "i_MMInf")
        Wk[k] <-  1/x$nodes[[k]]$mu
      else
      {
        if (class(x$nodes[[k]]) == "i_MM1" ||
            (class(x$nodes[[k]]) == "i_MMC" && x$nodes[[k]]$c == 1)
           )
          {
            #print("entrando en node mm1")
            Wk[k] <- (1 + Lk[k]) / x$nodes[[k]]$mu
          }
        else
        {
          num_mmc <- num_mmc + 1
          #print(paste("num_mmc:", num_mmc))
          
          #calculate probabilities
          mclass[[num_mmc]] <- CalcProb(i, x$nodes[[k]]$c,
            x$nodes[[k]]$mu, prob_est[k], Throughput, mclass[[num_mmc]]
          )

          aux1 <- 1/(x$nodes[[k]]$mu * x$nodes[[k]]$c)
          sum <- 1 + Lk[k]
  
          z <- 1
          while (z <= x$nodes[[k]]$c - 1)
          {
            sum <- sum + 
             ( (x$nodes[[k]]$c -(z-1) -1) * (mclass[[num_mmc]][z, i]) )
            z <- z + 1
          }
          Wk[k] <- aux1 * sum
        }
      }
              
      tmp <- tmp + (prob_est[k] * Wk[k])
      #print(paste("wi[", k, "]: ", wi[k]))
      #print(paste("tmp: ", tmp))

      k <- k + 1
    }

    Throughput <- i / tmp
    #print(paste("throughput: ", throughput))
    
    # The vector with the n values is updated
    Throughputn[i] <- Throughput

    k <- 1
    while (k <= num_nodes)
    {
      Lk[k] <- Throughput * prob_est[k] * Wk[k]
      #print(paste("li[", k, "]: ", li[k]))
      k <- k + 1
    }

    i <- i + 1
  }

  Throughputk <- prob_est * Throughput
  
  ROk <- rep(0, num_nodes)

  k <- 1
  while (k <= num_nodes)
  {
    if (class(x$nodes[[k]]) == "i_MMInf")
      ROk[k] <- Lk[k]
    else if (class(x$nodes[[k]]) == "i_MM1")
      ROk[k] <- Throughputk[k] * (1/x$nodes[[k]]$mu)
    else #class i_MMC  
      ROk[k] <- Throughputk[k] * (1/(x$nodes[[k]]$mu * x$nodes[[k]]$c))
    k <- k + 1
  }

  W <- (x$n / Throughput) - x$z
  L <- x$n - (Throughput * x$z)

  if (x$operational)
  {
    for (i in (1:num_nodes))
    {
      x$nodes[[i]]$mu <- x$nodes[[i]]$mu * x$prob[i]
      Throughputk[i] <- Throughputk[i] * x$prob[i]
    }
  }

  res <-
    list(
      Inputs = x,
      Throughput = Throughput,
      L = L,
      W = W,
      ROk = ROk,
      Throughputk = Throughputk,
      Wk = Wk,
      Lk = Lk,
      Throughputn = Throughputn
    )
  
  class(res) <- "o_CJN"
  res
}


QueueingModelApprox <- function(x, ...)
{
  num_nodes <- length(x$nodes)

  Throughputn <- numeric()

  if (class(x$prob) == "matrix")
  {
    ident <- diag(dim(x$prob)[1])
    const <- matrix(data=1, nrow=dim(x$prob)[1], ncol=1)
    all1 <- matrix(data=1, ncol=dim(x$prob)[2], nrow=dim(x$prob)[1])
    prob_est <- t(solve(t(x$prob + all1 - ident), const))
  }
  else
  {
    if (x$operational) #Visit ratios as counts of repetitions has been given
    {
      prob_est <- rep(1, num_nodes)
      #We have to "correct" the mu values
      for (i in (1:num_nodes))
        x$nodes[[i]]$mu <- x$nodes[[i]]$mu / x$prob[i]     
    }
    else
      prob_est <- x$prob
  }
  
  # array initialization
  Wk <- rep(0, num_nodes)
  Lk <- rep(x$n / num_nodes, num_nodes)
  Throughputk <- rep(0, num_nodes)
  Throughput <- 0
  
  finIter <- FALSE
  numIter <- 1

  while (!finIter)
  {
    Arrk <- ((x$n-1)/x$n) * Lk
  
    acum <- 0
    for (i in (1:num_nodes))
    {
      if (class(x$nodes[[i]]) == "i_MMInf")
        Wk[i] <- 1/x$nodes[[i]]$mu
      else
        Wk[i] <- 1/x$nodes[[i]]$mu * (1 + Arrk[i])

      acum <- acum + Wk[i] 
    }

    Throughput <- x$n / (x$z + acum)
    Throughputn[numIter] <- Throughput

    LkAux <- Throughput * Wk

    if (sum(abs((LkAux - Lk)) < x$tol) == num_nodes)
      finIter <- TRUE
    else
      {
        Lk <- LkAux
        numIter <- numIter + 1
      }
    }
  
  Throughputk <- prob_est * Throughput
  
  ROk <- rep(0, num_nodes)

  for (i in (1:num_nodes))
  {
    if (class(x$nodes[[i]]) == "i_MMInf")
      ROk[i] <- Lk[i] 
    else
      ROk[i] <- Throughputk[i] * (1/x$nodes[[i]]$mu)
  }

  W <- (x$n / Throughput) - x$z
  L <- x$n - (Throughput * x$z)

  if (x$operational)
  {
    for (i in (1:num_nodes))
    {
      x$nodes[[i]]$mu <- x$nodes[[i]]$mu * x$prob[i]
      Throughputk[i] <- Throughputk[i] * x$prob[i]
    }
  }

  res <-
    list(
      Inputs = x,
      Throughput = Throughput,
      L = L,
      W = W,
      ROk = ROk,
      Throughputk = Throughputk,
      Wk = Wk,
      Lk = Lk,
      Throughputn = Throughputn
    )
  
  class(res) <- "o_CJN"
  res
}



Inputs.o_CJN      <- function(x, ...) { x$Inputs }
Throughput.o_CJN  <- function(x, ...) { x$Throughput }
L.o_CJN           <- function(x, ...) { x$L }
W.o_CJN           <- function(x, ...) { x$W }
ROk.o_CJN         <- function(x, ...) { x$ROk }
Throughputk.o_CJN <- function(x, ...) { x$Throughputk }
Lk.o_CJN          <- function(x, ...) { x$Lk }
Wk.o_CJN          <- function(x, ...) { x$Wk }
Throughputn.o_CJN <- function(x, ...) { x$Throughputn }

Report.o_CJN <- function(x, ...)
{
  reportSingleClass(x)
}


summary.o_CJN <- function(object, ...)
{ 
  aux <- list(el=summarySingleClass(object))
  class(aux) <- "summary.o_CJN"
  aux
}


print.summary.o_CJN  <- function(x, ...)
{
  print_summary(x, ...)
}
