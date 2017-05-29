#######################################################################################
## MultiClass Open Network
#######################################################################################

NewInput.MCON <- function(classes, vLambda, nodes, vType, vVisit, vService)
{
  nds <- list(classes=classes, vLambda=vLambda, nodes=nodes, vType=vType, vVisit=vVisit, vService=vService)
  class(nds) <- "i_MCON"
  nds
}


CheckInput.i_MCON <- function(x, ...)
{

  MCON_vLambda_negatives <- "Some lambda has a negative value. Lambda has to be zero or positive"
  MCON_vService_negatives <- "Some service time is negative. Service time has to be zero or negative"
  MCON_lenght_vType_nodes <- "The lenght of Vtype vector doesn't coincide with nodes"
  x_class_MCON <- "The class of x has to be i_MCON"
  x_anomalous <- "Some parameter has a anomalous value" 
  MCON_dimension_visit_service <- "The matrix vVisit and the matrix vService has to have the same dimension"
  MCON_vVisit_negatives <- "Some visit has a negative value. Visits has to be zero or positive"
  MCON_vVisit_class_matrix <- "vVisit has to be of class matrix"
  MCON_vService_class_matrix <- "vService has to be of class matrix"
  MCON_dim_vVisit_nodes_vLambda <- "The dimension of the vVisit matrix doesn't coincide with the dimension of vLambda and nodes"
  MCON_vType_wrong <- "The types for the nodes has to be \"Q\" or \"D\""
  MCON_vlambda_classes_wrong <- "The number of elements of the vector vLambda has to be equal to classes"


  if (
    is.anomalous(x$vLambda) || is.anomalous(x$nodes) || is.anomalous(x$vType) ||
    is.anomalous(x$vVisit) || is.anomalous(x$vService)
  )
    stop(x_anomalous)

  if (class(x) != "i_MCON")
    stop(x_class_MCON)

  # Check negatives in parameters
  if (checkNegative(x$vLambda))
    stop(MCON_vLambda_negatives)

  if (checkNegative(x$vVisit))
    stop(MCON_vVisit_negatives)

  if (checkNegative(x$vService))
    stop(MCON_vService_negatives)

  if (x$classes != length(x$vLambda)) 
    stop(MCON_vlambda_classes_wrong)

  if (length(x$vType) != x$nodes)
    stop(MCON_lenght_vType_nodes)

  dimVisit <- dim(x$vVisit)


  if (sum(dimVisit == dim(x$vService)) != 2)
    stop(MCON_dimension_visit_service)

  if (class(x$vVisit) != "matrix")
    stop(MCON_vVisit_class_matrix)

  if (class(x$vService) != "matrix")
    stop(MCON_vService_class_matrix)

  if (sum(dimVisit == c(x$classes, x$nodes)) != 2)
    stop(MCON_dim_vVisit_nodes_vLambda)

  #vService has to has at least one element positive
  j <- 1
  while (j <= x$nodes)
  {
    if (sum(x$vService[, j]) <= 0)
      stop("At least some service time has to be greater than zero at each node")

    j <- j + 1
  }

  #vVisit has to has at least one element positive
  j <- 1
  while (j <= x$nodes)
  {
    if (sum(x$vVisit[, j]) <= 0)
      stop("At least some visit has to be greater than zero at each node")

    j <- j + 1
  }

  i <- 1
  while (i <= x$nodes)
  {
    if (x$vType[i] != "Q" && x$vType[i] != "D")
      stop(MCON_vType_wrong)
     
    ro_aux <- sum(x$vLambda * x$vVisit[, i] * (x$vService[, i])) 

    if (x$vType[i] == "Q" && ro_aux >= 1 )
      stop(paste("The processing capacity of node ", i, " is saturated. The utilization is: ", ro_aux * 100, "%", sep=""))

    i <- i + 1
  }

  
}

QueueingModel.i_MCON <- function(x, ...)
{
  CheckInput(x)

  Throughputck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  ROck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Wck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Lck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Wc <- rep(0, x$classes)
  Lc <- rep(0, x$classes)
  Throughputc <- x$vLambda
  Throughput <- sum(x$vLambda)

  Wk <- rep(0, x$nodes)
  Lk <- rep(0, x$nodes)
  Throughputk <- rep(0, x$nodes)
  ROk <- rep(0, x$nodes)

  for (nd in (1:x$nodes))
  {
    Sk <- x$vService[, nd]
    Throughputck[, nd] <- x$vLambda * x$vVisit[, nd]
    ROck[, nd] <- Throughputck[, nd] * Sk
    inf_i <- 1 - sum(ROck[, nd])

    if (x$vType[nd] == "Q")
    {
      Wck[, nd] <- (x$vVisit[, nd] * Sk)/inf_i
      Lck[, nd] <- ROck[, nd]/inf_i
    }
    else
    {
      Wck[, nd] <- (x$vVisit[, nd] * Sk)
      Lck[, nd] <- ROck[, nd]
    }
  }

  # values for class
  W <- 0
  for (cla in (1:x$classes))
  {
    Wc[cla] <- sum(Wck[cla, ])
    Lc[cla] <- sum(Lck[cla, ])
    W <- W + (Wc[cla] * Throughputc[cla])
  }

  W <- W / Throughput

  for (nd in (1:x$nodes))
  {
    Lk[nd] <- sum(Lck[, nd])
    Throughputk[nd] <- sum(Throughputck[, nd])
    ROk[nd] <- sum(ROck[, nd])
    Wk[nd] <- sum(Wck[, nd] * Throughputc)
  }
  
  Wk <- Wk / Throughput

  L <- sum(Lc)

  res <-
    list(
      Inputs=x,
      W=W,
      Throughput=Throughput,
      L=L,
      Wc=Wc,
      Throughputc=Throughputc,      
      Lc=Lc,
      ROk=ROk,
      Wk=Wk,
      Throughputk=Throughputk,      
      Lk=Lk,
      ROck=ROck,
      Wck=Wck,
      Throughputck=Throughputck,     
      Lck=Lck
    )

  class(res) <- "o_MCON"
  res    

}

Inputs.o_MCON       <- function(x, ...) { x$Inputs }
W.o_MCON            <- function(x, ...) { x$W }
L.o_MCON            <- function(x, ...) { x$L }
Throughput.o_MCON   <- function(x, ...) { x$Throughput }
Wc.o_MCON           <- function(x, ...) { x$Wc }
Lc.o_MCON           <- function(x, ...) { x$Lc }
Throughputc.o_MCON  <- function(x, ...) { x$Throughputc }
ROk.o_MCON          <- function(x, ...) { x$ROk }
Wk.o_MCON           <- function(x, ...) { x$Wk }
Lk.o_MCON           <- function(x, ...) { x$Lk }
Throughputk.o_MCON  <- function(x, ...) { x$Throughputk }
ROck.o_MCON         <- function(x, ...) { x$ROck }
Wck.o_MCON          <- function(x, ...) { x$Wck }
Lck.o_MCON          <- function(x, ...) { x$Lck }
Throughputck.o_MCON <- function(x, ...) { x$Throughputck }


Report.o_MCON <- function(x, ...)
{   
  reportMultiClass(x)
}


summary.o_MCON <- function(object, ...)
{ 
  aux <- list(el=summaryMultiClass(object))
  class(aux) <- "summary.o_MCON"
  aux
}


print.summary.o_MCON  <- function(x, ...)
{
  print_summary(x, ...)
}
