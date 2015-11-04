#' @title Stochastic Dynamic Programming (Simplified, multi-objective)
#' @description Derives the optimal release policy based on storage state and within-year period only. The objective is to minimise a penalty cost function based on water supply, spill, and water level.
#' @param Q             time series object. Net inflows to the reservoir.
#' @param capacity      numerical. The reservoir storage capacity (must be the same volumetric unit as Q and the target release).
#' @param target        numerical. The target release constant.
#' @param vol_targ      numerical. The target storage volume constant (as proportion of capacity).
#' @param R_max         numerical. The maximum controlled release.
#' @param weights       vector of length 3 indicating weighting to be applied to release, spill and water level objectives respectively.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param loss_exp      vector of length 3 indicating the exponents on release, spill and water level deviations from target. Default exponents are c(2,2,2).
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param tol           numerical. The tolerance for policy convergence. The default value is 0.990.
#' @param rep_rrv       logical. If TRUE then reliability, resilience and vulnerability metrics are computed and returned.
#' @return Returns a list that includes: the optimal policy as an array of release decisions dependent on storage state, month/season, and current-period inflow class; the Bellman cost function based on storage state, month/season, and inflow class; the optimized release and storage time series through the training inflow data; the flow discretization (which is required if the output is to be implemented in the rrv function); and, if requested, the reliability, resilience, and vulnerability of the system under the optimized policy. 
#' @references Loucks, D.P., van Beek, E., Stedinger, J.R., Dijkman, J.P.M. and Villars, M.T. (2005) Water resources systems planning and management: An introduction to methods, models and applications. Unesco publishing, Paris, France.
#' @seealso \code{\link{sdp}} for deterministic Dynamic Programming 
#' @examples \donttest{storage_cap <- 4 * mean(aggregate(ResX_inflow.ts)) # set storage ratio of 4 years
#' demand <- 0.8 * mean(ResX_inflow.ts) # set draft ratio of 0.8
#' optimal.releases <- sdp_simple_multi(ResX_inflow.ts, capacity = storage_cap,
#' target = demand, R_max = 2 * demand)
#' }
#' @import stats
#' @export
sdp_simple_multi <- function (Q, capacity, target, R_max, vol_targ = 0.75, weights = c(4, 2, 1),
                              S_disc = 1000, R_disc = 10,
                        Q_disc = c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0),
                        loss_exp = c(2,2,2), S_initial = 1, plot = TRUE, tol = 0.99, rep_rrv = FALSE){
  
  frq <- frequency(Q)
  if (is.ts(Q)==FALSE) stop("Q must be seasonal time series object with frequency of 12 or 4")
  if (frq != 12 && frq != 4) stop("Q must have frequency of 4 or 12")  
  if (start(Q)[2] != 1){
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
  }
  if(end(Q)[2] != frq){
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
  }
  
  Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)                                        
  #n_Qcl <- length(Q_disc) - 1
  Q.probs <- diff(Q_disc)
  Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                       probs = Q_disc[-1] - (Q.probs / 2))
  S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)                   
  R_disc_x <- seq(from = 0, to = R_max, by = R_max / R_disc)
  
  if (target %in% R_disc_x == FALSE) {
    warning("Warning: target not contained in R_disc")
  }
  
  Shell.array <- array(0,dim=c(length(S_states),length(R_disc_x),length(Q.probs)))
  R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))             
  #Q_class.mat <- matrix(nrow=length(Q_month_mat[,1]),ncol=frq)
  
  Cost_to_go <- vector("numeric",length=length(S_states))
  Results_mat <- matrix(0,nrow=length(S_states),ncol=12)
  R_policy <- matrix(0,nrow=length(S_states),ncol=12)
  Bellman <- R_policy
  R_policy_test <- R_policy
  message(paste0("policy converging... (>", tol,")"))
  
  # POLICY OPTIMIZATION----------------------------------------------------------------
  
  repeat{
    for (t in frq:1){
      R.cstr <- sweep(Shell.array, 3, Q_class_med[,t], "+") +
        sweep(Shell.array, 1, S_states, "+")                               
      R.star[which(R.star > R.cstr)] <- NaN                              
      Deficit.arr <- (target - R.star) / target                                           
      Deficit.arr[which(Deficit.arr < 0)] <- 0
      Cost_arr <- ( (abs(Deficit.arr)) ^ loss_exp[1])                          
      S.t_plus_1 <- R.cstr - R.star
      
      Spill_costs <- S.t_plus_1 - capacity
      Spill_costs[which(Spill_costs < 0)] <- 0
      Spill_costs <- (Spill_costs / quantile(Q, 0.95)) ^ loss_exp[2]
      
      S.t_plus_1[which(S.t_plus_1 > capacity)] <- capacity
      Vol_costs <- abs(((S.t_plus_1 - vol_targ * capacity) / (vol_targ * capacity))) ^ loss_exp[3]
      Implied_S_state <- round(1 + (S.t_plus_1 / capacity)
                               * (length(S_states) - 1))
      
      Cost_arr <- weights[1] * Cost_arr + weights[2] * Spill_costs  + weights[3] * Vol_costs
      
      Cost_to_go.arr <- array(Cost_to_go[Implied_S_state],
                              dim = c(length(S_states), length(R_disc_x) , length(Q.probs)))       
      
      Min_cost_arr <- Cost_arr + Cost_to_go.arr
      Min_cost_arr_weighted <- sweep(Min_cost_arr, 3, Q.probs, "*")
      Min_cost_expected <- apply(Min_cost_arr_weighted, c(1, 2), sum)
      
      Bellman[,t] <- Cost_to_go
      Cost_to_go <- apply(Min_cost_expected, 1, min, na.rm = TRUE)
      Results_mat[,t] <- Cost_to_go
      R_policy[,t] <- apply(Min_cost_expected, 1, which.min)
    }
    message(sum(R_policy == R_policy_test) /
              (frq * length(S_states)))   
    if (sum(R_policy == R_policy_test) /
        (frq * length(S_states)) > tol){
      break
    }
    
    R_policy_test <- R_policy
  }
  
  # ===================================================================================
  
  # POLICY SIMULATION------------------------------------------------------------------
  
  
  S <- vector("numeric",length(Q) + 1); S[1] <- S_initial * capacity    
  R_rec <- vector("numeric",length(Q))                      
  Spill <- vector("numeric",length(Q))
  
  for (yr in 1:nrow(Q_month_mat)) {
    for (month in 1:frq) {
      t_index <- (frq * (yr - 1)) + month   
      S_state <- which.min(abs(S_states - S[t_index]))
      Qx <- Q_month_mat[yr,month]
      #Q_class <- which.min(abs(as.vector(Q_class_med[,month] - Qx)))
      R <- R_disc_x[R_policy[S_state,month]]
      R_rec[t_index] <- R
      if ( (S[t_index] - R + Qx) > capacity) {
        S[t_index + 1] <- capacity
        Spill[t_index] <- S[t_index] - R + Qx - capacity
      }else{
        if ( (S[t_index] - R + Qx) < 0) {
          S[t_index + 1] <- 0
          R_rec[t_index] <- S[t_index] + Qx
        }else{
          S[t_index + 1] <- S[t_index] - R + Qx
        }
      }
    }
  }
  
  R_policy <- (R_policy - 1) / R_disc
  S <- ts(S[2:length(S)],start = start(Q),frequency = frq)
  R_rec <- ts(R_rec, start = start(Q), frequency = frq)
  Spill <- ts(Spill, start = start(Q), frequency = frq)
  
  if(plot) {
    layout(1:3)
    plot(S, ylab = "storage", ylim = c(0, capacity)); abline(h = vol_targ * capacity, lty = 2)
    plot(R_rec, ylab = "release", ylim = c(0, R_max)); abline(h = target, lty = 2)
    plot(Spill, ylab = "spill")
  }
  
  total_release_cost <- sum((R_rec/target)[which((R_rec/target) <  1)] ^ loss_exp[1])
  total_spill_cost <- sum((Spill / quantile(Q, 0.95)) ^ loss_exp[2])
  total_volume_cost <- sum(((S - vol_targ * capacity) / (vol_targ * capacity)) ^ loss_exp[3])
  total_weighted_cost <- weights[1] * total_release_cost + weights[2] * total_spill_cost + weights[3] * total_volume_cost 
  costs <- list(total_release_cost, total_spill_cost, total_volume_cost, total_weighted_cost)
  names(costs) <- c("total_release_cost", "total_spill_cost", "total_volume_cost", "total_weighted_cost")
  
  if (rep_rrv == TRUE){
    
    # COMPUTE RRV METRICS FROM SIMULATION RESULTS---------------------------------------
    
    deficit <- ts(round(1 - (R_rec / target),5), start = start(Q), frequency = frequency(Q))
    rel_ann <- sum(aggregate(deficit, FUN = mean) == 0) /
      length(aggregate(deficit, FUN = mean))
    rel_time <- sum(deficit == 0) / length(deficit)
    rel_vol <- sum(R_rec) / (target * length(deficit))
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
    
    results <- list(R_policy, Bellman, S, R_rec, Spill, rel_ann, rel_time, rel_vol, resilience, vulnerability, Q_disc, costs)
    names(results) <- c("release_policy", "Bellman", "storage", "releases", "spill", "annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability", "flow_disc", "costs")
    
    
    
  } else {
    results <- list(R_policy, Bellman, S, R_rec, Spill, Q_disc, costs)
    names(results) <- c("release_policy", "Bellman", "storage", "releases", "spill", "flow_disc", "total_costs")
  }
  
  return(results)
}