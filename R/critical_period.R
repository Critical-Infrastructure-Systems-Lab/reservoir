#' @title Critical (drawdown) period.
#' @description For computing the critical period of a reservoir from its storage time series. The critical period is defined as the length of time for the reservoir to go from full to empty (without spilling in between). Input storage should fill and empty at least once for correct calculation of critical period.
#' @param x     vector or time series object giving the reservoir storage behaviour. The sequence should contain at least one zero values for the critical period to be estimated.
#' @param report  character string giving the critical period to report... "shortest" (default) returns shortest of all critical periods (i.e., the most severe drought); "average" returns the mean of all critical periods; "longest" returns the longest critical drawdown period; "all" returns all critical drawdown periods in the time series.   
#' @return Returns the critical drawdown period of the reservoir, giving the shortest time from full to empty by default.
#' @examples storage_behavior <- simRes(Q=resX$Q_Mm3, target = 0.9*mean(resX$Q_Mm3),
#'                                      capacity=1000, plot=FALSE)$storage
#' plot(storage_behavior)
#' 
#' ## longest critical period
#' critPeriod(storage_behavior) # or critPeriod(x, "longest")
#' 
#' ## shortest critical period
#' critPeriod(storage_behavior, "shortest")
#' 
#' ## average critical period
#' critPeriod(storage_behavior, "average")
#' 
#' ## all critical periods
#' critPeriod(storage_behavior, "all")
#' @import stats
#' @export
critPeriod <- function(x, report) {
  if (missing(report)) report <- "longest"
  if (is.ts(x) == FALSE && is.vector(x) == FALSE)
    stop("x must be time series or vector object")
  
  if(min(x)>0) warning("Cannot calculate actual critical period for a reservoir that never empties.")
  
  if (length(which(x == min(x))) == 1){
    hit_bottom_periods <- which(x == min(x))
  } else {
    hit_bottom_periods <- which(x == min(x))[-which(c(NA, diff(which(x == min(x)))) == 1)]
  }
  #hit_bottom_periods <- which(x == min(x))[-which(c(NA, diff(which(x == min(x)))) == 1)]
  CDPs <- vector("numeric", length(hit_bottom_periods))
  for (i in 1:length(hit_bottom_periods)){
    CDPs[i] <- hit_bottom_periods[i] - max(which(x == max(x))[which(x == max(x)) < hit_bottom_periods[i]])
  }
  
  if(length(CDPs > 1)){
    CDPs <- CDPs[!c(FALSE, diff(CDPs)==diff(hit_bottom_periods))]
  }
  
  if(report=="longest") CDP <- max(CDPs)
  if(report=="average") CDP <- round(mean(CDPs))
  if(report=="shortest") CDP <- min(CDPs)
  if(report=="all") CDP <- CDPs
  return(CDP)
}