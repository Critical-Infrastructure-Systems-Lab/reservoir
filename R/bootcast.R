#' @title Develop seasonal streamflow forecasts for an observed period ("hindcasts") using the k nearest neighbour (kNN) bootstrap in forecast mode.
#' @description For seasonal streamflow "hindcasts" based solely on persistence.
#' @param Q              time series object - seasonal streamflow rate or streamflow totals.
#' @param H              integer giving the number of seasons in the forecast horizon (e.g., H = 3 means a three month forecast if the frequency of Q is 12).
#' @param start_yr       start year (if omitted, forecast will be for one year)
#' @param end_yr         end_year (if omitted, forecast will be for last year of input)
#' @param k              integer. The k parameter of the kNN Bootstrap (i.e., number of nearest neighbors from which to sample). If left blank k = n ^ 0.5., where n is the number of years in the input data.
#' @param d              integer. The d parameter of the kNN Bootstrap (i.e., number of previous time periods to inform the model). If left blank d = 1. 
#' @param sampling_mode  used to define period of observed streamflow data (Q) from which to sample the forecasts. "past" will use the period of data prior to the forecast start year. "all" will use all data within Q. "adapt" (the default) updates the sample each year to all years prior to forecast release.
#' @param plot           logical. If TRUE (the default) the function will return a plot of the forecasts against observed flow during the period start_yr : end_yr.
#' @return Returns the critical drawdown period of the reservoir, giving the shortest time from full to empty by default.
#' @examples # prepare three month ahead forecasts of resX$Q_Mm3 for the period 1995 to 2000...
#' Q <- resX$Q_Mm3
#' Q_fcast <- bootcast(Q, H = 3, start_yr = 1995, end_yr = 2000)
#' @import stats
#' @export
bootcast <- function(Q, H, start_yr, end_yr, k, d, sampling_mode, plot){
  
  if (is.ts(Q)==FALSE) stop("Q must be seasonal time series object")
  
  frq <- frequency(Q)
  len <- length(Q)
  cyc <- cycle(Q)
  
  if (missing(H)) stop("Please provide forecast horizon H")
  if (H < 1 || H > frq || H%%1 != 0) stop(paste0("H must be an integer between 1 and ", frq))
  if (missing(start_yr)) stop("Please provide start_yr")
  if (start(Q)[1] >= start_yr) stop("start_yr must be >= 1 yr ahead of the start year of Q")
  if (missing(end_yr)){
    end_yr <- end(Q)[1]
    warning(paste0("No forecast end year entered... end_yr taken as ", end_yr))
  }
  if (missing(k)) k <- round((len/frq) ^ 0.5)
  if (missing(d)) d <- 1
  if (missing(sampling_mode)) sampling_mode = "adapt"
  if (missing(plot)) plot <- TRUE
  if (d > frq) stop("d must be <= frequency(Q)")
  
  # DEFINE FUNCTION FOR EUCLIDEAN DISTANCE CALCULATION
  getEucDist <- function(D_t, D_i) {
    EucDist <- sqrt(sum((D_t - D_i)^2))
    return(EucDist)
  }
  
  # DEFINE PROBABILITY MASS FUNCTION FOR THE KERNAL
  kernal.weights <- (1 / seq(1:k)) / (sum(1 / seq(1:k)))
  
  # set up X variable, which defines data set to be sampled from...
  if (sampling_mode == "all"){
    Qs <- Q
  } else {
    Qs <- window(Q, end = c(start_yr - 1, frq))
  }
  
  cycQs <- cycle(Qs)
  
  X <- ts(0, start=c(start_yr,1), end = c(end_yr,12), frequency = frq)
  
  Q_Qs_X <- cbind(Q, Qs, X)
  forecasts <- matrix(NA, nrow = length(Q), ncol = H)
  
  for (t in which(Q_Qs_X[,3] == 0)){
    
    D_i <- Qs[(t - d):(t - 1)]
    locs <- which(cycQs==cyc[t])
    locs<- locs[locs - d > 0]                # Remove instances for which insufficient training data exist
    locs <- locs[locs + H <= length(Qs)]      # Remove instances for which insufficient samples exist
    locs <- locs[locs != t]                   # Remove current instance
    D_locs <- matrix(nrow=length(locs), ncol = d)
    F_locs <- matrix(nrow=length(locs), ncol = H)
    for(i in 1:length(locs)){
      D_locs[i,] <- (locs[i] - d):(locs[i] - 1)  
      F_locs[i,] <- locs[i]:(locs[i] + H - 1)
    }
    D_t <- matrix(Qs[D_locs], ncol = d)
    F_t <- matrix(Qs[F_locs], ncol = H)
    distances <- apply(D_t, 1, getEucDist, D_i = D_i)
    
    x_pot <- F_t[match(1:k, rank(distances, ties.method = "random")),]
    if (k == 1){
      forecasts[t,] <- x_pot
    } else {
      if (H == 1){
        forecasts[t] <- x_pot[sample(1:k, 1, prob = kernal.weights)]
      } else {
        forecasts[t,] <- x_pot[sample(1:k, 1, prob = kernal.weights),]
      }
    }
    if (sampling_mode == "adapt") Qs <- window(Q, end = time(Q)[t])
  }
  
  if (plot==TRUE){
    if (H == 1){
      forecasts <- ts(forecasts[,1], start = start(Q), frequency = frq)
      plot(window(Q, start = c(start_yr, 1)), col = "darkgrey", ylim = c(0, max(Q)), ylab = "Q")
      lines(window(forecasts, start = c(start_yr, 1)), col = "red")
    } else {
      plot(window(Q, start = c(start_yr,1)), col = "darkgrey", ylim = c(0,max(Q)), ylab = "Q")
      for(i in 1:len){
        lines(ts(forecasts[i,], start = time(Q)[i],frequency=frq), col = "red")
      }
    }
  }

  return(matrix(forecasts[is.na(forecasts)==FALSE], ncol = H))
  
}