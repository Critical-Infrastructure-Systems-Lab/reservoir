#' @title Reliability, resilience, and vulnerability analysis for water supply reservoirs
#' @description Computes time-based, annual, and volumetric reliability, as well as resilience and dimensionless vulnerability for a single reservoir.
#' @param Q                  time series or vector. The net inflows to the reservoir.
#' @param target             numerical constant, or a time series or vector of the target releases. Must be the same length as Q is given as a vector or time series.
#' @param capacity           numerical. The reservoir capacity. Should be same volumetric unit as Q and R.
#' @param surface_area       numerical. The reservoir surface area at full capacity. Must be in square kilometers (km^2), or Mm^2.
#' @param max_depth          numerical. The maximum water depth of the reservoir at the dam at maximum capacity. Must be in meters. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only.
#' @param evap               vector or time series object of length Q, or a numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param double_cycle       logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series) 
#' @param plot               logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial          numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param policy             list. The output of the SDP function. If omitted, Standard Operating Policy is assumed.
#' @return Returns reliability, resilience and vulnerability metrics based on supply deficits.
#' @references McMahon, T.A., Adeloye, A.J., Zhou, S.L. (2006) Understanding performance measures of reservoirs, Journal of Hydrology 324 (359-382)
#' @examples # Compare reliability, resilience and vulnerability for two operating policies (SOP and SDP).
#' rrv(resX$Q_Mm3, capacity = 20*resX$cap_Mm3, target = 0.95 * mean(resX$Q_Mm3))
#' pol_Markov <- sdp_supply(resX$Q_Mm3, capacity = 20 * resX$cap_Mm3,
#' target = 0.95 * mean(resX$Q_Mm3), Markov = TRUE)
#' rrv(resX$Q_Mm3, capacity = 20*resX$cap_Mm3, target = 0.95 * mean(resX$Q_Mm3), policy = pol_Markov)
#' @import stats
#' @export
rrv <- function(Q, target, capacity, double_cycle = FALSE,
                surface_area, max_depth, evap,
                plot = TRUE, S_initial = 1, policy) {
    

    frq <- frequency(Q)
    if(length(target)==1){
      target <- rep(target, length(Q))
    }
    if (double_cycle) {
        Q <- ts(c(Q, Q), start = start(Q), frequency = frq)
        target <- c(target, target)
    }
    
    if (missing(evap)) {
      evap <- rep(0, length(Q))
    }
    if(length(evap) == 1) {
      evap <- rep(evap, length(Q))
    }
    if (length(evap) != length(Q)){
      stop("Evaporation must be either a vector (or time series) length Q, or a single numeric constant")
    }
    
    if (missing(surface_area)) {
      surface_area <- 0
    }
    
    if (missing(max_depth)){
      c <- sqrt(2) / 3 * (surface_area * 10 ^ 6) ^ (3/2) / (capacity * 10 ^ 6)
      GetLevel <- function(c, V){
        y <- (6 * V / (c ^ 2)) ^ (1 / 3)
        return(y)
      }
      GetArea <- function(c, V){
        Ay <- (((3 * c * V) / (sqrt(2))) ^ (2 / 3))
        return(Ay)
      }
    } else {
      c <- 2 * capacity / (max_depth * surface_area)
      GetLevel <- function(c, V){
        y <- max_depth * (V / (capacity * 10 ^ 6)) ^ (c / 2)
        return(y)
      }
      GetArea <- function(c, V){
        Ay <- ((2 * (capacity * 10 ^ 6)) / (c * max_depth * (V / (capacity * 10 ^ 6)) ^ (c / 2))) * ((V / (capacity * 10 ^ 6)) ^ (c / 2)) ^ (2 / c)
        Ay[which(is.nan(Ay) == TRUE)] <- 0
        return(Ay)
      }
    }
    
    GetEvap <- function(s, q, r, ev){
      e <- GetArea(c, V = s * 10 ^ 6) * ev / 10 ^ 6
      n <- 0
      repeat{
        n <- n + 1
        s_plus_1 <- max(min(s + q - r - e, capacity), 0)
        e_x <- GetArea(c, V = ((s + s_plus_1) / 2) * 10 ^ 6) * ev / 10 ^ 6
        if (abs(e_x - e) < 0.001 || n > 20){
          break
        } else {
          e <- e_x
        }
      }
      return(e)
    }
    
    S <- vector("numeric", length = length(Q) + 1); S[1] <- capacity * S_initial
    R <- vector("numeric", length = length(Q))
    E <- vector("numeric", length(Q))
    y <- vector("numeric", length(Q))
    Spill <- vector("numeric", length(Q))
    
    # SIMULATION--------------------------------------------------------------------

    if (missing(policy)){
      for (t in 1:length(Q)) {

        #x <- Q[t] - target[t] + S[t]
        
        E[t] <- GetEvap(s = S[t], q = Q[t], r = target[t], ev = evap[t])
        y[t] <- GetLevel(c, S[t] * 10 ^ 6)
        
        if ( (S[t] - target[t] + Q[t] - E[t]) > capacity) {
          S[t + 1] <- capacity
          Spill[t] <- S[t] - target[t] + Q[t] - capacity - E[t]
          R[t] <- target[t]
        } else {
          if (S[t] - target[t] + Q[t] - E[t] < 0) {
            S[t + 1] <- 0
            R[t] <- max(0, S[t] + Q[t] - E[t])
          } else {
            S[t + 1] <- S[t] +  Q[t] - target[t] - E[t]
            R[t] <- target[t]
          }
        }
        
        
        #if (x < 0) {
        #    S[t + 1] <- 0
        #    R[t] <- Q[t] + S[t]
        #} else {
        #    if (x > capacity) {
        #        S[t + 1] <- capacity
        #        R[t] <- target[t]
        #    } else {
         #       S[t + 1] <- x
        #        R[t] <- target[t]
         #   }
        #}
    }
    }
    
    if (!missing(policy)){
      if (is.list(policy) == FALSE) stop("The policy must be the full output returned by the sdp_supply function")
      
      Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)
      
      if (length(dim(policy$release_policy)) == 2){
        S_disc <- length(policy$release_policy[,1]) - 1
      } else {
        S_disc <- length(policy$release_policy[,1,1]) - 1
      }
      
      S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)
      Q.probs <- diff(policy$flow_disc)
      Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                           probs = policy$flow_disc[-1] - (Q.probs / 2))
      
      for (yr in 1:nrow(Q_month_mat)) {
        for (month in 1:frq) {
          t_index <- (frq * (yr - 1)) + month   
          S_state <- which.min(abs(S_states - S[t_index]))
          Qx <- Q_month_mat[yr, month]
          
          if (length(dim(policy$release_policy)) == 2){
            Rx <- target[t_index] * policy$release_policy[S_state, month]
          } else {
            Q_class <- which.min(abs(as.vector(Q_class_med[,month] - Qx)))
            Rx <- target[t_index] * policy$release_policy[S_state, Q_class, month]
          }
          
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
    
    deficit <- ts(round(1 - (R / target),5), start = start(Q), frequency = frq)
    rel_ann <- sum(aggregate(deficit, FUN = mean) == 0) /
      length(aggregate(deficit, FUN = mean))
    rel_time <- sum(deficit == 0) / length(deficit)
    rel_vol <- sum(R) / sum(target)
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
    
    results <- list(rel_ann, rel_time, rel_vol, resilience, vulnerability,
                    ts(S[1:length(S) - 1], start = start(Q), frequency = frq),
                    ts(R, start = start(Q), frequency = frq))
    names(results) <- c("annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability", "storage", "releases")
    if (plot) {
        plot(ts(S[1:length(S) - 1], start = start(Q), frequency = frq), ylab = "Storage", ylim = c(0, capacity))
        plot(ts(R, start = start(Q), frequency = frq), ylab = "Release", ylim = c(0, max(target)))
    }
    
    return(results)
}
