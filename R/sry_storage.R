#' @title Storage-Reliability-Yield (SRY) relationships: Storage computation
#' @description Returns the required storage for given inflow time series, yield, and target time-based reliability. Assumes standard operating policy. Storage is computed iteratively using the bi-section method.
#' @param Q               time series or vector. The net inflows to the reservoir. This must be a time series object or vector of the net inflow volumes.
#' @param yield           numerical.  (must be same volumetric unit as Q and R).
#' @param reliability     numerical. The required time-based reliability.
#' @param profile         a vector of factors with length = frequency(Q). Represents within-year demand profile. Defaults to constant release if left blank.
#' @param plot            logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial       numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param max.iterations  Maximum number of iterations for yield computation.
#' @param double_cycle    logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series)
#' @return Returns the required storage capacity necessary to supply specified yield with specified reliability.
#' @examples # Determine the required storage for 95 % reliability and yield equal to 80 % of the mean inflow.
#' storage(ResX_inflow.ts, yield = 0.8*mean(ResX_inflow.ts), reliability = 0.95)
#' @import stats
#' @export
storage <- function(Q, yield, reliability, profile = rep(1, frequency(Q)),
                    plot = TRUE, S_initial = 1, max.iterations = 50, double_cycle = FALSE) {
    if (length(profile) != frequency(Q))
        stop("profile must have length equal to the time series frequency - e.g., (0.8, 1.2, 1.2, 0.8) for quarterly data")
    if (reliability < 0 || reliability > 1)
        stop("Reliability must be between 0 and 1")
    if (length(yield) > 1)
        stop("yield must be a scalar, not a vector")
    if (double_cycle) {
        Q <- ts(c(Q, Q), start = start(Q), frequency = frequency(Q))
    }
    if (reliability == 1){
      reliability <- 0.9999 # reliability of 1 allows convergence to stop at over-large storage
    }
  
    S <- vector("numeric", length = length(Q) + 1)
    R <- vector("numeric", length = length(Q))
    min_storage <- 0
    max_storage <- 100 * mean(Q) * frequency(Q)  #Provides upper bound of 100 years' storage
    i <- 0
    repeat {
        i <- i + 1
        mid_storage <- (min_storage + max_storage) / 2
        S[1] <- mid_storage * S_initial
        R.target <- rep(yield, length(Q)) * rep(profile, length(Q) / frequency(Q))
        for (t in 1:length(Q)) {
            x <- Q[t] - R.target[t] + S[t]
            if (x < 0) {
                S[t + 1] <- 0
                R[t] <- Q[t] + S[t]
            } else {
                if (x > mid_storage) {
                  S[t + 1] <- mid_storage
                  R[t] <- R.target[t]
                } else {
                  S[t + 1] <- x
                  R[t] <- R.target[t]
                }
            }
        }
        
        reliability_0 <- 1 - (sum((R / R.target) < 1) / length(Q))
        if (reliability_0 > reliability) {
            max_storage <- mid_storage
        }
        if (reliability_0 < reliability) {
            min_storage <- mid_storage
        }
        
        if (round(reliability_0, 3) == reliability) {
            break
        }
        if (i >= max.iterations) {
            break
        }
        
    }
    
    if (plot) {
        plot(ts(S[2:length(S)], start = start(Q),
                frequency = frequency(Q)), ylab = "Storage",main = (paste0("Storage = ", 
            round(mid_storage, 3), " at ", round(reliability_0, 3) * 100, " % reliability")))
        plot(ts(R, start = start(Q), frequency = frequency(Q)),
             ylim = c(0, max(R)), ylab = "Water supplied")
    }
    message(paste0("Converged after ", i,
                 " iterations and time-based reliability equal to ", reliability_0))
    return(mid_storage)
    } 
