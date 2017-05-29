############################################################
## class methods
############################################################

CheckInput     <- function(x, ...) UseMethod("CheckInput")
QueueingModel  <- function(x, ...) UseMethod("QueueingModel")
Inputs         <- function(x, ...) UseMethod("Inputs")
RO             <- function(x, ...) UseMethod("RO")
Lq             <- function(x, ...) UseMethod("Lq")
VNq            <- function(x, ...) UseMethod("VNq")
Wq             <- function(x, ...) UseMethod("Wq")
VTq            <- function(x, ...) UseMethod("VTq")
L              <- function(x, ...) UseMethod("L")
VN             <- function(x, ...) UseMethod("VN")
W              <- function(x, ...) UseMethod("W")
VT             <- function(x, ...) UseMethod("VT")
Wqq            <- function(x, ...) UseMethod("Wqq")
Pn             <- function(x, ...) UseMethod("Pn")
Qn             <- function(x, ...) UseMethod("Qn")
Lqq            <- function(x, ...) UseMethod("Lqq")
Throughput     <- function(x, ...) UseMethod("Throughput")
WWs            <- function(x, ...) UseMethod("WWs")
SP             <- function(x, ...) UseMethod("SP")
Throughputc    <- function(x, ...) UseMethod("Throughputc")
Throughputk    <- function(x, ...) UseMethod("Throughputk")
Throughputck   <- function(x, ...) UseMethod("Throughputck")
Throughputn    <- function(x, ...) UseMethod("Throughputn")
Throughputcn   <- function(x, ...) UseMethod("Throughputcn")
Lc             <- function(x, ...) UseMethod("Lc")
Lk             <- function(x, ...) UseMethod("Lk")
Lck            <- function(x, ...) UseMethod("Lck")
Wc             <- function(x, ...) UseMethod("Wc")
Wk             <- function(x, ...) UseMethod("Wk")
Wck            <- function(x, ...) UseMethod("Wck")
ROk            <- function(x, ...) UseMethod("ROk")
ROck           <- function(x, ...) UseMethod("ROck")
Report         <- function(x, ...) UseMethod("Report")


############################################################
## Error Messages
############################################################
ALL_mu_positive      <- "mu must be greater than zero"
ALL_lambda_zpositive <- "lambda must be equal or greater than zero"
ALL_n_integer        <- "the number of clients must be an integer number"
ALL_c_integer        <- "the number of servers (c) must be an integer number"
ALL_k_integer        <- "k must be a integer number"
ALL_c_warning        <- "c has to be at least one"
ALL_k_warning        <- "k has to be at least one"
ALL_k_c              <- "k must be equal or greater than the number of servers c"

