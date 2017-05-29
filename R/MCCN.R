#######################################################################################
## MultiClass Closed Network
#######################################################################################

NewInput.MCCN <- function(classes, vNumber, vThink, nodes, vType, vVisit, vService, method=1, tol=0.01)
{
  nds <- list(classes=classes, vNumber=vNumber, vThink=vThink, nodes=nodes, vType=vType, vVisit=vVisit, vService=vService, method=method, tol=tol)

  class(nds) <- "i_MCCN"
  nds
}


CheckInput.i_MCCN <- function(x, ...)
{
  MCCN_x_class <- "The class of x has to be i_MCCN"
  MCCN_x_anomalous <- "Some parameter has a anomalous value"
  

  MCCN_vNumber_negatives_or_zero <- "Some class number of clients are zero or negative. The number of clientes has to be positive"
  
  MCCN_vService_negatives <- "Some service time is negative"
  MCCN_vThink_negatives <- "Some think time is negative"

  MMCN_classes_at_least_one <- "The number of clasess has to be one or greater"
  MMCN_nodes_at_least_one <- "The number of nodes has to be one or greater"
  MMCN_vVisit_negatives <- "Some visit is negative"

  MCCN_vNumber_classes <- "The length of vector vNumber doesn't coincide with classes"
  MCCN_lenght_vType_nodes <- "The lenght of Vtype vector doesn't coincide with nodes"
  MCCN_length_vNumber_vThink <- "The length of the vector vNumber does not coincide with lenght of vector vThink"  


  MCCN_dimension_visit_service <- "The matrix vVisit and the matrix vService has to have the same dimension"
  
  MCCN_vVisit_class_matrix <- "vVisit has to be of class matrix"
  MCCN_vService_class_matrix <- "vService has to be of class matrix"

  MCCN_dim_vVisit_nodes_vNumber <- "The dimension of the vVisit matrix doesn't coincide with the dimension of vNumber and nodes"
  MCCN_vType_wrong <- "The types for the nodes has to be \"Q\" or \"D\""
  MCCN_method_values <- "The x$method has to be 0 (exact) or 1 (aprox)"
  MCCN_tol_value <- "The x$tol has to be positive"

  if (
    is.anomalous(x$classes) || is.anomalous(x$vNumber) || is.anomalous(x$vThink) ||
    is.anomalous(x$nodes) || is.anomalous(x$vType) || is.anomalous(x$vVisit) || is.anomalous(x$vService) ||
    is.anomalous(x$method) || is.anomalous(x$tol)
  )
    stop(MCCN_x_anomalous)

  if (class(x) != "i_MCCN")
    stop(MCCN_x_class)

  # Check negatives, zero or one in parameters
  if (checkNegativeOrZero(x$vNumber))
    stop(MCCN_vNumber_negatives_or_zero)

  if (checkNegative(x$vService))
    stop(MCCN_vService_negatives)

  if (checkNegative(x$vThink))
    stop(MCCN_vThink_negatives)

  if (checkAtLeastOne(x$classes))
    stop(MMCN_classes_at_least_one)

  if (checkAtLeastOne(x$nodes))
    stop(MMCN_nodes_at_least_one)

  if (checkNegative(x$vVisit))
    stop(MMCN_vVisit_negatives) 

  if (x$method != 0 && x$method != 1)
    stop(MCCN_method_values)

  if (!(x$tol > 0))
    stop(MCCN_tol_value)

  # dimension and lengths
  if (length(x$vType) != x$nodes)
    stop(MCCN_lenght_vType_nodes)

  if (length(x$vNumber) != x$classes)
    stop(MCCN_vNumber_classes)

  if (length(x$vNumber) != length(x$vThink))
    stop(MCCN_length_vNumber_vThink)

  if (sum(dim(x$vVisit) == dim(x$vService)) != 2)
    stop(MCCN_dimension_visit_service)

  # classes
  if (class(x$vVisit) != "matrix")
    stop(MCCN_vVisit_class_matrix)

  if (class(x$vService) != "matrix")
    stop(MCCN_vService_class_matrix)

  if (sum(dim(x$vVisit) == c(length(x$vNumber), x$nodes)) != 2)
    stop(MCCN_dim_vVisit_nodes_vNumber)

  # vType has the correct types
  i <- 1
  while (i <= x$nodes)
  {
    if (x$vType[i] != "Q" && x$vType[i] != "D")
      stop(MCCN_vType_wrong)

    i <- i + 1
  }
  
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

}


QueueingModel.i_MCCN <- function(x, ...)
{
  CheckInput(x)

  if (x$method == 0)
    QueueingModelMCCNExact(x, ...)
  else
    QueueingModelMCCNApprox(x, ...)  
}


QueueingModelMCCNExact <- function(x, ...)
#check <- function(x, ...)
{
  Throughputc <- rep(0, x$classes)
  Throughputck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  ROck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Wck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Lck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)  
  
  ROk <- rep(0, x$nodes)
  Lk <- rep(0, x$nodes)
  Throughputk <- rep(0, x$nodes)
  Wk <- rep(0, x$nodes)

  LkAux <- list()

  # from the last to the first of x$vNumber
  i <- x$classes
  while (i > 0)
  {
    LkAux <- c(LkAux, list(seq(1, x$vNumber[i] + 1)))  
    i <- i - 1
  }

  # Generate the matrix of elements to iterate
  # the elements has to be reversed when iterate over it
  iter <- as.matrix(expand.grid(LkAux))
  
  numRowsIter <- nrow(iter)

  orderCol <- rep(0, numRowsIter)

  # calculate the order num of each row
  for (i in 1:numRowsIter)
    orderCol[i] <- sum(iter[i, ])

  #attach it to the iter matrix
  iter <- cbind(iter, orderCol)

  # Order ir by the OrderCol column
  iter <- iter[order(orderCol), ]

  #orderCol has to be replaced to the column already ordered
  orderCol <- iter[, ncol(iter)]

  LkV <- array(0, dim=c(x$vNumber + 1, x$nodes))
  Throughputcn <- array(0, dim=c(x$vNumber + 1, x$classes))
  
  for (n in (3:(sum(x$vNumber+1))))
  {
    
    #Get the population to iterate
    population <- iter[orderCol == n, ]

    cardin <- nrow(population)
    if (is.null(cardin))
      cardin <- 1

    #Now, getting each elem
    for (subpopulation in 1:cardin)
    {
      if (cardin == 1)
        elem <- as.vector(population)
      else
        elem <- as.vector(population[subpopulation, ])
   
      ## Fixed elem to remove the order value
      elem <- rev(elem[-length(elem)])
      #revElem <- rev(elem)   
   
      for (cla in (1:x$classes))
      {
        vaux <- elem
        vaux[cla] <- max(elem[cla]-1, 1)
   
        for (server in (1:x$nodes))
        {
          demand <- x$vVisit[cla, server] * (x$vService[cla, server])

          if (x$vType[server] == "D")
            Wck[cla, server] <- demand
          else        
            Wck[cla, server] <- demand * (1 + LkV[array(data=c(vaux, server), dim=c(1, x$classes+1))])
        }
      }

      for (cla in (1:x$classes))
      {
        Throughputc[cla] <- (elem[cla]-1) / (x$vThink[cla] + sum(Wck[cla, ]))
        Throughputcn[array(data=c(elem, cla), dim=c(1, x$classes+1))] <- Throughputc[cla]
      }
 
      for (server in (1:x$nodes))
        LkV[array(data=c(elem, server), dim=c(1, x$classes+1))] <- sum(Throughputc * Wck[, server])
    }
  }

  Wc <- (x$vNumber/Throughputc) - x$vThink
  Lc <- x$vNumber - (Throughputc * x$vThink)

  for (i in (1:x$nodes))
  {
    ThAux <- Throughputc * x$vVisit[, i]
    Throughputck[ , i] <- ThAux
    ROck[, i] <- ThAux * (x$vService[, i])
    Lck[, i] <- Throughputc * Wck[, i]
  }

  Throughput <- 0
  L <- 0

  for (i in (1:x$nodes))
  {
    ROk[i] <- sum(ROck[, i])
    Lk[i] <- sum(Lck[, i])
    L <- L + Lk[i]
    Throughputk[i] <- sum(Throughputck[, i])
    Throughput <- Throughput + Throughputk[i]
    Wk[i] <- sum(Wck[, i] * Throughputc)
  }
  
  Wk <- Wk / Throughput 

  W <- (sum(Wc * Throughputc)) / Throughput

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
      Lck=Lck,
      Throughputcn=Throughputcn
    )

  class(res) <- "o_MCCN"
  res    

}


#############################################################################
#############################################################################
QueueingModelMCCNApprox <- function(x, ...)
{
  Throughputc <- rep(0, x$classes)
  Throughputck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  ROck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Wck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)
  Lck <- matrix(data=0, nrow=x$classes, ncol=x$nodes)  
  Ack <- matrix(data=0, nrow=x$classes, ncol=x$nodes)

  ROk <- rep(0, x$nodes)
  Lk <- rep(0, x$nodes)
  Throughputk <- rep(0, x$nodes)
  Wk <- rep(0, x$nodes)

  Throughputcn <- array(0, dim=c(x$vNumber + 1, x$classes))

  for (i in (1:x$classes))
    Lck[i, ] <- x$vNumber[i]/x$nodes

  # To ensure that at least one iteration is done, other option is to use repeat and break
  AuxLck <- Lck + 2*x$tol

  while (sum( (AuxLck - Lck >= x$tol) > 0))
  {

    for (i in (1:x$classes))
      Ack[i, ] <- ((x$vNumber[i] - 1)/x$vNumber[i]) * Lck[i, ] + colSums(rbind(Lck[-i, ], 0))

    for (i in (1:x$classes))
    {
      if (x$vType[i] == "D")
        Wck[i, ] <- x$vVisit[i, ] * x$vService[i, ]
      else
        Wck[i, ] <- (x$vVisit[i, ] * x$vService[i, ]) * (1 + Ack[i, ])
    }

    AuxLck <- Lck

    for (i in (1:x$classes))
    {
      Throughputc[i] <- x$vNumber[i] / (x$vThink[i] + sum(Wck[i, ]))
      Lck[i, ] <- Throughputc[i] * Wck[i, ]
    }
    
    #print("LckF")
    #print(Lck)

  }
  
  Wc <- (x$vNumber/Throughputc) - x$vThink  
  Lc <- Wc * Throughputc

  for (i in (1:x$nodes))
  {
    Throughputck[, i] <- Throughputc * x$vVisit[, i]
    ROck[, i] <- Throughputck[, i] * x$vService[, i]    
  }

  Throughput <- 0
  L <- 0

  for (i in (1:x$nodes))
  {
    ROk[i] <- sum(ROck[, i])
    Lk[i] <- sum(Lck[, i])
    L <- L + Lk[i]
    Throughputk[i] <- sum(Throughputck[, i])
    Throughput <- Throughput + Throughputk[i]
    Wk[i] <- sum(Wck[, i] * Throughputc)
  }

  Wk <- Wk / Throughput 

  W <- (sum(Wc * Throughputc)) / Throughput

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
      Lck=Lck,
      Throughputcn=Throughputcn
    )

  class(res) <- "o_MCCN"
  res    

}


Inputs.o_MCCN       <- function(x, ...) { x$Inputs }
L.o_MCCN            <- function(x, ...) { x$L }
W.o_MCCN            <- function(x, ...) { x$W }
Throughput.o_MCCN   <- function(x, ...) { x$Throughput }
Lc.o_MCCN           <- function(x, ...) { x$Lc }
Wc.o_MCCN           <- function(x, ...) { x$Wc }
Throughputc.o_MCCN  <- function(x, ...) { x$Throughputc }
ROk.o_MCCN          <- function(x, ...) { x$ROk }
Lk.o_MCCN           <- function(x, ...) { x$Lk }
Wk.o_MCCN           <- function(x, ...) { x$Wk }
Throughputk.o_MCCN  <- function(x, ...) { x$Throughputk }
ROck.o_MCCN         <- function(x, ...) { x$ROck }
Lck.o_MCCN          <- function(x, ...) { x$Lck }
Wck.o_MCCN          <- function(x, ...) { x$Wck }
Throughputck.o_MCCN <- function(x, ...) { x$Throughputck }
Throughputcn.o_MCCN <- function(x, ...) { x$Throughputcn }

Report.o_MCCN <- function(x, ...)
{   
  reportMultiClass(x)  
}


summary.o_MCCN <- function(object, ...)
{ 
  aux <- list(el=summaryMultiClass(object))
  class(aux) <- "summary.o_MCCN"
  aux
}


print.summary.o_MCCN  <- function(x, ...)
{
  print_summary(x, ...)
}
