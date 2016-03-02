#' @title Quick and dirty stochastic generation of seasonal streamflow replicates for a single site.
#' @description Generates seasonal time series using a numerically-fitted PARMA(1,1) model. The function features automatic procedures for transforming seasonal sub-series to normal and deseasonalizing data.
#' @param Q             time series object with seasonal resolution (e.g., frequency = 2, 3, 4, 6 or 12 for monthly data). 
#' @param reps          integer. The number of replicates to be generated.  The default is 100.
#' @param repyrs        integer. The length of each replicate in years. The default is length(Q) / frequency(Q).
#' @param adjust        logical. If TRUE (the default) the final output time series X will be coerced for 0 <= X <= 1.2*max(Q). 
#' @return Returns a multi time series object containing synthetic streamflow replicates.
#' @examples 
#' Q <- resX$Q_Mm3
#' replicates <- dirtyreps(Q)
#' mean(replicates); mean(Q)
#' sd(replicates); sd(Q)
#' layout(1:4)
#' plot(Q, main = "Recorded series", ylim = c(0,1500))
#' plot(replicates[,1], main = "Replicate # 1", ylim = c(0,1500))
#' plot(replicates[,20], main = "Replicate # 20", ylim = c(0,1500))
#' plot(replicates[,99], main = "Replicate # 99", ylim = c(0,1500))
#' 
#' @import stats
#' @export
dirtyreps <- function(Q, reps = 100, repyrs = length(Q) / frequency (Q), adjust = TRUE){
  
  #if(min(Q) < 0){
  #  stop("Q cannot be log-transformed to normal; contains negative or zero values")
  #}
  
  frq <- frequency(Q)
  len <- length(Q)
  cyc <- cycle(Q)
  
  if (start(Q)[2] != 1){
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
  }
  if(end(Q)[2] != frq){
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
  }

  # CREATE FUNCTION TO ESTIMATE PARAMETER OF THE LOG TRANSFORM
  getlogtransparam <- function(x){
    a <- (min(x) * max(x) - median(x) ^ 2) / (min(x) + max(x) - 2 * median(x))
    if(a >= min(x)) a <- min(x) - 1
    return(a)
  }
  
  # CREATE FUNCTION TO LOG-TRANSFORM SERIES
  getQtrans <- function(x){
    a <- (min(x) * max(x) - median(x) ^ 2) / (min(x) + max(x) - 2 * median(x))
    if(a >= min(x)) a <- min(x) - 1
    xT <- log(x - a)
    return(xT)
  }
  
  # CREATE FUNCTIONS FOR PARMA FITTING
  get_params <- function(x){
    model <- lm(x[,1] ~ 0 + x[,2] + x[,4])
    ph <- model[[1]][[1]] 
    th <- model[[1]][[2]]
    return(c(ph,th))
  }
  get_resids <- function(x){
    model <- lm(x[,1] ~ 0 + x[,2] + x[,4])
    return(model[[2]])
  }
  
  
  # TRANSFORM THE SEASONAL SUBSERIES INDIVIDUALLY
  Qtr <- ts(as.vector(t(matrix(unlist(tapply(Q, cyc, getQtrans)),ncol = frq))),
               start = start(Q), frequency = frq)
  
  ## REMOVE MEAN AND DIVIDE BY STANDARD DEVIATION TO GET WHITE NOISE
  Qtr_ds <- ((Qtr - rep(tapply(Qtr, cyc, mean),len / frq)) / 
               rep(tapply(Qtr, cyc, sd), len / frq))
  
  ## INITIALISE RESIDUALS FOR LEAST SQUARES PARMA FITTING
  M <- list()
  Y <- Qtr_ds[cyc == 1][2:(len / frq)]
  X <- Qtr_ds[cyc == frq][1:((len / frq) - 1)]
  ar1_model <- lm(Y ~ 0 + X)
  E <- as.numeric(ar1_model[[2]])
  M[[1]] <- cbind(Y, X, E, vector("numeric", length(X)))
  for(v in 2:frq){
    Y <- Qtr_ds[cyc == v]
    X <- Qtr_ds[cyc == (v - 1)]
    ar1_model <- lm(Y ~ 0 + X)
    E <- as.numeric(ar1_model[[2]])
    M[[v]] <- cbind(Y, X, E, vector("numeric", length(X)))
  }
  paramsx <- matrix(0, frq, 2)
  
  #ITERATE
  n = 0
  repeat{
    n = n + 1
    for(i in 1:frq){
      if(i==1){
        M[[i]][,4] <- head(M[[frq]][,3], -1)
      }else if(i==2){
        M[[i]][,4] <- c(0 , M[[i-1]][,3])  
      }else{
        M[[i]][,4] <- M[[i-1]][,3]
      }
    }
    params <- matrix(unlist(lapply(M, get_params)),ncol=2,byrow=TRUE)
    resids <- lapply(M, get_resids)
    for (i in 1:frq){
      M[[i]][,3] <- resids[[i]]
    }
    if (sum(params - paramsx > 0.01) == 0 | n > 50) break
    paramsx <- params
  }
  
  # GENERATE WHITE NOISE FROM THE FITTED PARMA MODEL
  params <- cbind(params, unlist(lapply(resids,sd)))
  colnames(params) <- c("phi", "theta", "st_dev")
  lenX <- repyrs * frq * reps
  Xtr_ds <- vector("numeric", lenX)
  period <- rep(1:frq, lenX/frq)
  error <- rnorm(1, 0, params[1, 3])
  Xtr_ds[1] <- params[1, 1] *  tail(Qtr_ds,1) + params[1, 2] * tail(resids[[frq]],1) + error
  for (t in 2:length(Xtr_ds)){
    preverror <- error
    error <- rnorm(1, 0, params[period[t], 3])
    Xtr_ds[t] <- params[period[t], 1] * Xtr_ds[t - 1] + params[period[t], 2] * preverror + error
  }
  
  ## ADD SEASONS BACK
  Xtr <- as.vector(Xtr_ds * rep(tapply(Qtr, cyc, sd), lenX / frq)
                   + rep(tapply(Qtr, cyc, mean), lenX / frq))
  
  ## REVERSE THE LOG TRANSFORM
  a <- tapply(Q, cyc, getlogtransparam)
  mins <- tapply(Q, cyc, min)
  a[which(a>=mins)] <- mins[which(a>=mins)] - 1
  a <- rep(a, length(Xtr) / frq)
  X <- ts(matrix((exp(Xtr)) + a, ncol = reps), frequency = frq)
  
  if (adjust){
    X[which(X < 0)] <- 0
    X[which(X > 1.2 * max(Q))] <- 1.2 * max(Q)
  }
  
  return(X)
}