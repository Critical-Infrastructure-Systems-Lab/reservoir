#' @title Reliability, resilience, and vulnerability analysis for water supply reservoirs
#' @description Computes time-based, annual, and volumetric reliability, as well as resilience and dimensionless vulnerability for a single reservoir.
#' @param X                  time series object representing the release time series of a reservoir.
#' @param target            numerical constant (or a vector of length = length(X)). If omitted then the trigger constant will be < max(X).
#' @return Returns reliability, resilience and vulnerability metrics based on supply deficits.
#' @examples # Compare reliability, resilience and vulnerability for two operating policies (SOP and SDP).
#' @import stats
#' @export
rrv <- function(X, target) {
  
  if(is.ts(X) == FALSE) stop("X must be a time series object... use ts(X, ...)")
  frq <- frequency(X)
  len <- length(X)
  
  if(length(target) == 1) {
    target <- rep(target, len)
  }
  if(length(target) != len) {
    stop(paste0("target must be a scalar constant or vector or time series of length = ",
                len, " for this input time series, X"))
  }
  
  deficit <- ts(round(1 - (X / target),5), start = start(X), frequency = frq)
  rel_ann <- sum(aggregate(deficit, FUN = mean) <= 0) /
    length(aggregate(deficit, FUN = mean))
  rel_time <- sum(deficit <= 0) / length(deficit)
  rel_vol <- sum(X) / sum(target)
  fail.periods <- which(deficit > 0)
  if (length(fail.periods) == 0) {
    resilience <- NA
    vulnerability <- NA
  } else {
    if (length(fail.periods) == 1) {
      resilience <- 1
      vulnerability <- max(deficit)
    } else {
      resilience <- (sum(diff(which(deficit > 0)) > 1) + 1) / (length(which(deficit > 0)))
      fail.refs <- vector("numeric", length = length(fail.periods))
      fail.refs[1] <- 1
      for (j in 2:length(fail.periods)) {
        if (fail.periods[j] > (fail.periods[j - 1] + 1)) {
          fail.refs[j] <- fail.refs[j - 1] + 1
        } else {
          fail.refs[j] <- fail.refs[j - 1]
        }
      }
      n.events <- max(fail.refs)
      event.starts <- by(fail.periods, fail.refs, FUN = min)
      event.ends <- by(fail.periods, fail.refs, FUN = max)
      max.deficits <- vector("numeric", length = n.events)
      for (k in 1:n.events) {
        max.deficits[k] <- max(deficit[event.starts[k]:event.ends[k]])
      }
      vulnerability <- mean(max.deficits)
    }
  }
  
  #===============================================================================
  
  results <- c(rel_ann, rel_time, rel_vol, resilience, vulnerability)
  names(results) <- c("annual_reliability",
                      "time_based_reliability", "volumetric_reliability",
                      "resilience", "vulnerability")
  
  return(results)
}
