###############################################################
###############################################################
## Open Jackson networks
###############################################################
###############################################################

clambda <- function(x)
{
  res <- numeric()
  i <- 1
  while (i <= length(x))
  {
    res[i] <- x[[i]]$lambda
    i <- i + 1
  }
  cbind(res)
}


newNodes <- function(rawNodes, arrivals)
{
  res <- list()
  i <- 1
  while (i <= length(rawNodes))
  {
    rawNode = rawNodes[[i]]
    
    if (inherits(rawNode, "i_MM1"))
      res[[i]] <- NewInput.MM1(lambda = arrivals[i], mu = rawNode$mu, n = rawNode$n)
    else if (inherits(rawNode, "i_MMC"))
      res[[i]] <- NewInput.MMC(lambda = arrivals[i], mu = rawNode$mu, c = rawNode$c, n = rawNode$n)
    else if (inherits(rawNode, "i_MMInf"))
      res[[i]] <- NewInput.MMInf(lambda = arrivals[i], mu = rawNode$mu, n = rawNode$n)
    else 
      stop(paste(paste("Node ", i), "is not of class i_MM1, i_MMC or i_MMInf !!"))

    i <- i + 1
  }
  res
}


doModel <- function(x, newNodes, tLambda)
{
  prob <- numeric()
  Lk <- numeric()
  Wk <- numeric()
  ROk <- numeric()
  Throughputk <- numeric()
  is_prob_a_matrix <- (inherits(x$prob, "matrix"))
  totalL <- 0
  
  i <- 1
  while (i <= length(newNodes))
  {
    aux <- QueueingModel(newNodes[[i]])
    prob[i] <- Pn(aux)
 
    if (is_prob_a_matrix)
      Wk[i] <- W(aux) * (Throughput(aux) / tLambda)    
      # old code wrong: it doesn't sum W: Wk[i] <- W(aux)
    else 
      Wk[i] <- W(aux) * x$prob[i]

    auxL <- L(aux)
    Lk[i] <- auxL
    ROk[i] <- RO(aux)
    Throughputk[i] <- Throughput(aux)
    totalL <- totalL + auxL
    i <- i + 1
  }

  W <- totalL/tLambda
  Throughput <- sum(tLambda)

  res <-
    list(
      Inputs = x,
      Throughput = Throughput,
      L = totalL,
      W = W,
      ROk = ROk,
      Throughputk = Throughputk,      
      Lk = Lk,      
      Wk = Wk,
      Pn = prob
    )

  class(res) <- "o_OJN"
  res

}



CheckInput.i_OJN <- function(x, ...)
{
 x_class_OJN <- "x has to be of class i_OJN (Open Jackson Network)" 
 x_anomalous <- "x has some anomalous value. Check the value(s)."
 row_distinct_col <- "x$prob has the number of rows distinct of the number of columns."
 row_distinct_nodes <- "x$prob has distinct number of rows that the number of nodes in x$nodes."
 visit_ratios_wrong <- "x$prob contains a different number of visit ratios than x$nodes."
 prob_zero <- "If neither a routing x$prob is given nor a visit ratio vector, x$prob should be 0"
 all_lambda_equals <- "if visit ratios are given, all nodes must have the same lambda (the sum of all external arrivals)"
 
 is_prob_a_matrix <- (inherits(x$prob, "matrix"))

 if (is.anomalous(x$prob) || is.anomalous(x$nodes))
    stop(x_anomalous)

 if (!inherits(x, "i_OJN"))
   stop(x_class_OJN) 

 num_nodes <- length(x$nodes)

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

 i <- 1
 while (i <= num_nodes)
 {
   n = x$nodes[[i]]

   if (!inherits(n, c("i_MM1", "i_MMC", "i_MMInf")))
     stop(paste(paste("Node ", i), "is not of class i_MM1, i_MMC or i_MMInf!!"))
   
   if (!is_prob_a_matrix && (x$nodes[[i]]$lambda != x$nodes[[1]]$lambda))
     stop(all_lambda_equals)

   CheckInput(n)

   i <- i + 1
 }

}


QueueingModel.i_OJN <- function(x, ...)
{
  CheckInput(x)
  
  if (inherits(x$prob, "matrix"))
  {
    vlambda <- -clambda(x$nodes)
    tProb <- t(x$prob) 
    sol <- solve(tProb - diag(nrow=nrow(tProb)), vlambda)
    newNd <- newNodes(x$nodes, sol)  
    model <- doModel(x, newNd, -sum(vlambda))
  }
  else
  {
    lambda <- x$nodes[[1]]$lambda
    arrivals <- x$prob * lambda
    newNd <- newNodes(x$nodes, arrivals)
    model <- doModel(x, newNd, lambda)
  }
  
  model

}


NewInput.OJN <- function(prob=NULL, ...)
{
  NewInput2.OJN(prob, nodes(...))
}


NewInput2.OJN <- function(prob=NULL, nodes)
{
  nds <- list(prob=prob, nodes=nodes)
  class(nds) <- "i_OJN"
  nds
}



NewInput3.OJN <- function(vLambda, numNodes, vType, vVisit, vService, vChannel)
{

  nodes <- list()

  # Build each node
  for (i in 1:numNodes)
  {
    if (vType[i] == "Q")
    {
      if (vChannel[i] > 1)
        nodes <- c(nodes, list(NewInput.MMC(vLambda[i], 1/vService[i], vChannel[i])))
      else
        nodes <- c(nodes, list(NewInput.MM1(vLambda[i], 1/vService[i])))
    }
    else
      nodes <- c(nodes, list(NewInput.MMInf(vLambda[i], 1/vService[i])))
  }
  
  NewInput2.OJN(vVisit, nodes)
}

Inputs.o_OJN      <- function(x, ...) { x$Inputs }
Throughput.o_OJN  <- function(x, ...) { x$Throughput }
L.o_OJN           <- function(x, ...) { x$L }
W.o_OJN           <- function(x, ...) { x$W }
ROk.o_OJN         <- function(x, ...) { x$ROk }
Throughputk.o_OJN <- function(x, ...) { x$Throughputk }
Lk.o_OJN          <- function(x, ...) { x$Lk }
Wk.o_OJN          <- function(x, ...) { x$Wk }
Pn.o_OJN          <- function(x, ...) { x$Pn }


Report.o_OJN <- function(x, ...)
{
  reportSingleClass(x)  
}


summary.o_OJN <- function(object, ...)
{ 
  aux <- list(el=summarySingleClass(object))
  class(aux) <- "summary.o_OJN"
  aux
}


print.summary.o_OJN  <- function(x, ...)
{
  print_summary(x, ...)
}
