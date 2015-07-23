#' @title Reliability, resilience, and vulnerability analysis
#' @description Computes time-based, annual, and volumetric reliability, as well as resilience and dimensionless vulnerability for a single reservoir.
#' @param Q                  time series or vector. The net inflows to the reservoir.
#' @param R_target           time series or vector. The target release. Must be the same length as Q.
#' @param capacity           numerical. The reservoir capacity. Should be same volumetric unit as Q and R.
#' @param double_cycle       logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series) 
#' @param plot               logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial          numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param policy             list. The output of the SDP function. The default is (NULL) is Standard Operating Policy.
#' @return Returns reliability, resilience and vulnerability metrics based on supply deficits.
#' @references McMahon, T.A., Adeloye, A.J., Zhou, S.L. (2006) Understanding performance measures of reservoirs, Journal of Hydrology 324 (359-382)
#' @examples # Determine the reliability, resilience and vulnerability for reservoir on Holland Creek
#' demand <- rep(0.8 * mean(ResX_inflow.ts), length = length(ResX_inflow.ts))
#' storage_cap <- 2*mean(aggregate(ResX_inflow.ts))  # 2 years' storage 
#' rrv(ResX_inflow.ts, R_target = demand, capacity = storage_cap)
#' @importFrom graphics layout
#' @import stats
#' @export
rrv <- function(Q, R_target, capacity, double_cycle = FALSE,
                plot = TRUE, S_initial = 1, policy = NULL) {
    

    frq <- frequency(Q) 
    if (double_cycle) {
        Q <- ts(c(Q, Q), start = start(Q), frequency = frq)
        R_target <- c(R_target, R_target)
    }
    R <- vector("numeric", length = length(Q))
    S <- vector("numeric", length = length(Q) + 1)
    S[1] <- capacity * S_initial
    
    # SIMULATION--------------------------------------------------------------------

    if (is.null(policy) == TRUE){
    for (t in 1:length(Q)) {
        x <- Q[t] - R_target[t] + S[t]
        if (x < 0) {
            S[t + 1] <- 0
            R[t] <- Q[t] + S[t]
        } else {
            if (x > capacity) {
                S[t + 1] <- capacity
                R[t] <- R_target[t]
            } else {
                S[t + 1] <- x
                R[t] <- R_target[t]
            }
        }
    }
    }
    
    if (is.null(policy) == FALSE){
      if (is.list(policy) == FALSE) stop("The policy must be the output returned by the SDP function")
      
      Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)
      S_states <- seq(from = 0, to = capacity, by = capacity / (length(policy$release_policy[,1,1]) - 1))
      
      Q.probs <- diff(policy$flow_disc)
      Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                           probs = policy$flow_disc[-1] - (Q.probs / 2))
      
      for (yr in 1:nrow(Q_month_mat)) {
        for (month in 1:frq) {
          t_index <- (frq * (yr - 1)) + month   
          S_state <- which.min(abs(S_states - S[t_index]))
          Qx <- Q_month_mat[yr, month]
          Q_class <- which.min(abs(as.vector(Q_class_med[,month] - Qx)))

          Rx <- R_target[t_index] * policy$release_policy[S_state,Q_class,month]
          
          if ( (S[t_index] - Rx + Qx) > capacity) {
            S[t_index + 1] <- capacity
            R[t_index] <- Rx
          }else{
            if ( (S[t_index] - Rx + Qx) < 0) {
              S[t_index + 1] <- 0
              R[t_index] <- S[t_index] + Qx
            }else{
              S[t_index + 1] <- S[t_index] - Rx + Qx
              R[t_index] <- Rx
            }
          }
        }
      }
    }

    #===============================================================================
    
    # COMPUTE METRICS FROM SIMULATION RESULTS---------------------------------------
    
    deficit <- ts(round(1 - (R / R_target),5), start = start(Q), frequency = frq)
    rel_ann <- sum(aggregate(deficit, FUN = mean) == 0) /
      length(aggregate(deficit, FUN = mean))
    rel_time <- sum(deficit == 0) / length(deficit)
    rel_vol <- sum(R) / sum(R_target)
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
    
    results <- list(rel_ann, rel_time, rel_vol, resilience, vulnerability)
    names(results) <- c("annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability")
    if (plot) {
        layout(1:3)
        plot(Q, ylab = "inflow")
        plot(ts(S[2:length(S)], start = start(Q), frequency = frq), ylab = "storage", ylim = c(0, capacity))
        plot(ts(R, start = start(Q), frequency = frq), ylab = "release", ylim = c(0, max(R_target)))
    }
    
    return(results)
}
