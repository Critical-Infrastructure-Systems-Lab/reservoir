#' @title Storage-Reliability-Yield (SRY) relationships: Yield computation
#' @description Returns the yield for given inflow time series, reservoir capacity, and required time-based reliability. Assumes standard operating policy. Yield is computed iteratively using the bi-section method.
#' @param Q               a time series or vector  of net inflows to the reservoir.
#' @param capacity        numerical.  (must be same volumetric unit as Q and R).
#' @param reliability     numerical.  (must be same volumetric unit as Q and R).
#' @param profile         a vector of factors with length = frequency(Q). Represents within-year demand profile. Defaults to constant release if left blank.
#' @param plot            logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial       numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param max.iterations  Maximum number of iterations for yield computation.
#' @param double_cycle  logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series)
#' @return Returns the storage behaviour time series for the no-failure (Rippl) reservoir given net inflows Q and target release R.
#' @examples # Compute yield for 0.95 reliability
#' yield_ResX <- yield(ResX_inflow.ts, capacity = 100000, reliability = 0.95)
#' # Compute yield for quarterly time series with seasonal demand profile
#' 
#' quarterly.ts <- aggregate(ResX_inflow.ts, nfrequency = 4)
#' yield_ResX.quart <- yield(quarterly.ts,
#' capacity = 100000, reliability = 0.9, profile = c(0.8, 1.2, 1.2, 0.8))
#' @export
yield <- function(Q, capacity, reliability, profile = rep(1, frequency(Q)),
                  plot = TRUE, S_initial = 1, max.iterations = 50, double_cycle = FALSE) {
    
    if (length(profile) != frequency(Q)) 
        stop("profile must have length equal to the time series frequency - e.g., (0.8, 1.2, 1.2, 0.8) for quarterly data")
    if (reliability < 0 || reliability > 1) 
        stop("Reliability must be between 0 and 1")
  
    if (double_cycle) {
        Q <- ts(c(Q, Q), start = start(Q), frequency = frequency(Q))
    }
  
    S <- vector("numeric", length = length(Q) + 1)
    S[1] <- capacity * S_initial
    R <- vector("numeric", length = length(Q))
    min_yield <- 0
    max_yield <- 10 * mean(Q)
    i <- 0
    repeat {
        i <- i + 1
        mid_yield <- (min_yield + max_yield) / 2
        R.target <- rep(mid_yield, length(Q)) * rep(profile, length(Q) / frequency(Q))
        for (t in 1:length(Q)) {
            x <- Q[t] - R.target[t] + S[t]
            if (x < 0) {
                S[t + 1] <- 0
                R[t] <- Q[t] + S[t]
            } else {
                if (x > capacity) {
                  S[t + 1] <- capacity
                  R[t] <- R.target[t]
                } else {
                  S[t + 1] <- x
                  R[t] <- R.target[t]
                }
            }
        }
        reliability_0 <- 1 - (sum( (R / R.target) < 1) / length(Q))
        if (reliability_0 > reliability) {
            min_yield <- mid_yield
        }
        if (reliability_0 < reliability) {
            max_yield <- mid_yield
        }
        if (round(reliability_0, 3) == reliability) {
            break
        }
        if (i >= max.iterations) {
            break
        }
    }
    if (plot) {
        layout(1:2)
        plot(ts(S[2:length(S)], start = start(Q), frequency = frequency(Q)), ylab = "Storage")
    
    plot(ts(R, start = start(Q), frequency = frequency(Q)),
         ylim = c(0, max(R)), ylab = "Water supplied", main = (paste0("Yield = ", 
        round(mid_yield, 2), " at ", reliability_0 * 100, " % reliability")))
    }
    message(paste0("Converged after ", i,
                 " iterations and time-based reliability equal to ", reliability_0))
    return(mid_yield)
}
