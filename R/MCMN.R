#######################################################################################
## Multi Class Mixed Queueing Networks
#######################################################################################

# Open Class has to be defined first
NewInput.MCMN <- function(classes, vLambda, vNumber, vThink, nodes, vType, vVisit, vService, method=0, tol=0.01)
{
  nds <- list(classes=classes, vLambda=vLambda, vNumber=vNumber, vThink=vThink, nodes=nodes, vType=vType, vVisit=vVisit, vService=vService, method=method, tol=tol)
  class(nds) <- "i_MCMN"
  nds
}

# To check the params, the functions of the others classes are to be reused
CheckInput.i_MCMN <- function(x, ...)
{

  # Ckeck the class
  if (!inherits(x, "i_MCMN"))
    stop("The class of the object has to be i_MCMN")

  numOC <- length(x$vLambda)
  numCC <- length(x$vNumber)

  if (numOC < 1)
    stop("The number of Open Classes has to be at least one")

  if (numCC < 1)
    stop("The number of Closed Classes has to be at least one")
   
  if ((numCC + numOC) != x$classes)
    stop("The number of classes declared does not coincide with length of the lambda vector and the vNumber vector")

  openMod <- NewInput.MCON(numOC, x$vLambda, x$nodes, x$vType, x$vVisit[1:numOC, ], x$vService[1:numOC, ])

  closedMod <- NewInput.MCCN(numCC, x$vNumber, x$vThink, x$nodes,
     x$vType, x$vVisit[(numOC+1):x$classes, ], x$vService[(numOC+1):x$classes, ], x$method, x$tol)

  # Check each one of the models
  CheckInput(openMod)
  CheckInput(closedMod)
}


QueueingModel.i_MCMN <- function(x, ...)
{
  CheckInput(x)

  numOC <- length(x$vLambda)
  numCC <- length(x$vNumber)
  
  openInput <- NewInput.MCON(numOC, x$vLambda, x$nodes, x$vType, x$vVisit[1:numOC, ], x$vService[1:numOC, ])

  #Solve the open model
  openModel <- QueueingModel(openInput)

  openInflated <- 1 - ROk(openModel)

  closedInput <- NewInput.MCCN(numCC, x$vNumber, x$vThink, x$nodes,
     x$vType, x$vVisit[(numOC+1):x$classes, ],
     (x$vService[(numOC+1):x$classes, ])/matrix(openInflated, nrow=numCC, ncol=x$nodes, byrow=TRUE),
    x$method, x$tol)

  #Solve the closed model
  closedModel <- QueueingModel(closedInput)

  closedROck <- Throughputck(closedModel) * x$vVisit[(numOC+1):x$classes, ] * x$vService[(numOC+1):x$classes, ]

  openWck <- Wck(openModel) * (1 + t(array(Lk(closedModel), dim=c(x$nodes, numCC))))
  openLck <- Throughputck(openModel) * openWck 
  
  
  # Build the complete result
  Wck <- rbind(openWck, Wck(closedModel))
  Lck <- rbind(openLck, Lck(closedModel))
  Throughputck <- rbind(Throughputck(openModel), Throughputck(closedModel))
  ROck <- rbind(ROck(openModel), closedROck)

  Lc <- rowSums(Lck)
  Throughputc <- rowSums(Throughputck)
  ROc <- rowSums(ROck)

  Lk <- colSums(Lck)
  Throughputk <- colSums(Throughputck)
  ROk <- colSums(ROck)
  
  L <- sum(Lc)
  Throughput <- sum(Throughputc)
  W <- L / Throughput

  Wk <- colSums(Wck * array(Throughputc, dim=c(x$classes, x$nodes)))/Throughput
  Wc <- rowSums(Wck * matrix(data=Throughputk, nrow=x$classes, ncol=x$nodes, byrow = TRUE))/Throughput
  
  # W <- (Wc * Throughputc)/Throughput <-- to check that values are correct

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

  class(res) <- "o_MCMN"  
  res
}


Inputs.o_MCMN       <- function(x, ...) { x$Inputs }
L.o_MCMN            <- function(x, ...) { x$L }
W.o_MCMN            <- function(x, ...) { x$W }
Throughput.o_MCMN   <- function(x, ...) { x$Throughput }
Lc.o_MCMN           <- function(x, ...) { x$Lc }
Wc.o_MCMN           <- function(x, ...) { x$Wc }
Throughputc.o_MCMN  <- function(x, ...) { x$Throughputc }
ROk.o_MCMN          <- function(x, ...) { x$ROk }
Lk.o_MCMN           <- function(x, ...) { x$Lk }
Wk.o_MCMN           <- function(x, ...) { x$Wk }
Throughputk.o_MCMN  <- function(x, ...) { x$Throughputk }
ROck.o_MCMN         <- function(x, ...) { x$ROck }
Lck.o_MCMN          <- function(x, ...) { x$Lck }
Wck.o_MCMN          <- function(x, ...) { x$Wck }
Throughputck.o_MCMN <- function(x, ...) { x$Throughputck }

Report.o_MCMN <- function(x, ...)
{   
  reportMultiClass(x)  
}


summary.o_MCMN <- function(object, ...)
{ 
  aux <- list(el=summaryMultiClass(object))
  class(aux) <- "summary.o_MCMN"
  aux
}


print.summary.o_MCMN  <- function(x, ...)
{
  print_summary(x, ...)
}
