#' @title Stochastic Dynamic Programming (Simplified) for Hydropower Reservoir
#' @description Derives the optimal release policy based on storage state and within-year period only.
#' @param Q             time series object. Net inflows to the reservoir.
#' @param capacity      numerical. The reservoir storage capacity (must be the same volumetric unit as Q).
#' @param surface_area  numerical. The reservoir surface area
#' @param hydro_cap     numerical. The hydropower plant electric capacity (MW)
#' @param head          numerical. The maximum hydraulic head of the hydropower plant (m)
#' @param qmax          numerical. The maximum flow into the hydropower plant
#' @param depth         numerical. The reservoir maximum depth
#' @param dep_vol_curve string. The relationship between depth and volume of reservoir.
#' @param input_curve   data frame of 3 columns: depth, area, volume
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param tol           numerical. The tolerance for policy convergence. The default value is 0.990.
#' @param rep_rrv       logical. If TRUE then reliability, resilience and vulnerability metrics are computed and returned.
#' @return Returns a list that includes: the optimal policy as an array of release decisions dependent on storage state, month/season, and current-period inflow class; the Bellman cost function based on storage state, month/season, and inflow class; the optimized release and storage time series through the training inflow data; the flow discretization (which is required if the output is to be implemented in the rrv function); and, if requested, the reliability, resilience, and vulnerability of the system under the optimized policy. 
#' @references Loucks, D.P., van Beek, E., Stedinger, J.R., Dijkman, J.P.M. and Villars, M.T. (2005) Water resources systems planning and management: An introduction to methods, models and applications. Unesco publishing, Paris, France.
#' @seealso \code{\link{dp}} for deterministic Dynamic Programming 
#' @examples \donttest{storage_cap <- 4 * mean(aggregate(ResX_inflow.ts)) # set storage ratio of 4 years
#' demand <- 0.8 * mean(ResX_inflow.ts) # set draft ratio of 0.8
#' optimal.releases <- sdp_simple(ResX_inflow.ts, capacity = storage_cap, target = demand)
#' }
#' @import stats
#' @export
sdp_simple_hydro <- function (Q, capacity, capacity_live, surface_area, hydro_cap, head, qmax, depth,
                              dep_vol_curve = "l", input_curve, S_disc = 1000, R_disc = 10,
                              Q_disc = c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0),
                              S_initial = 1, plot = TRUE, tol = 0.99){
  
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
  
  if ( dep_vol_curve == "n") { 
    N <- 2 * capacity / (depth * surface_area)
  } else if ( dep_vol_curve == "f") {
    f <- sqrt(2) / 3 * (surface_area*10^6) ^ (3/2) / (capacity*10^6)
  }
  
  if (!missing(input_curve)){
    colnames(input_curve) <- c("depth", "area", "volume")
  }
  
  Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)                                        
  Q.probs <- diff(Q_disc)
  Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                       probs = Q_disc[-1] - (Q.probs / 2))
  
  S_states <- seq(from = capacity - capacity_live, to = capacity, by = capacity_live / S_disc)
  if ( !missing(qmax) ) {
    R_max <- qmax
  } else {
    R_max <- hydro_cap * 10^6 / (1000 * 9.81 * head) * (365/12 * 24 * 60 * 60) / 10^6 #10^6 m3/month
  }
  R_disc_x <- seq(from = 0, to = R_max, by = R_max / R_disc)
  
  Shell.array <- array(0,dim=c(length(S_states),length(R_disc_x),length(Q.probs)))
  R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))             
  
  Rev_to_go <- vector("numeric",length=length(S_states))
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
      R.star[which(R.star > (R.cstr - (capacity - capacity_live)))] <- NaN                              
      S.t_plus_1 <- R.cstr - R.star
      S.t_plus_1[which(S.t_plus_1 > capacity)] <- capacity
      
      if ( dep_vol_curve == "l"){
        H_arr <- (S.t_plus_1 + S_states) / 2 / surface_area
      } else if ( dep_vol_curve == "f") {
        H_arr <- ( 6 * (S.t_plus_1 + S_states) * 10^6 / 2 / f^2 ) ^ (1/3)
      } else if ( dep_vol_curve == "n") {
        H_arr <- depth * ( (S.t_plus_1 + S_states) / 2 / capacity ) ^ (N/2)
      } else if ( dep_vol_curve == "u") {
        H_arr <-  array(sapply( ((S.t_plus_1 + S_states) * 10^6 / 2), 
                                function (x) 
                                  ifelse(is.nan(x), NaN, input_curve$depth[ which.min(abs( input_curve$volume - x )) ]) ), 
                        dim=c(length(S_states),length(R_disc_x),length(Q.probs)))
      }
      
      Rev_arr <- R.star * H_arr
      Implied_S_state <- round(1 + ( ((S.t_plus_1 - (capacity - capacity_live)) / (capacity_live)) *
                                       (length(S_states) - 1)))
      Rev_to_go.arr <- array(Rev_to_go[Implied_S_state],
                             dim = c(length(S_states), length(R_disc_x) , length(Q.probs)))       
      
      Max_rev_arr <- Rev_arr + Rev_to_go.arr
      Max_rev_arr_weighted <- sweep(Max_rev_arr, 3, Q.probs, "*")
      Max_rev_expected <- apply(Max_rev_arr_weighted, c(1, 2), sum)
      
      Bellman[,t] <- Rev_to_go
      Rev_to_go <- apply(Max_rev_expected, 1, max, na.rm = TRUE)
      Results_mat[,t] <- Rev_to_go
      R_policy[,t] <- apply(Max_rev_expected, 1, which.max)
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
  
  
  S <- vector("numeric",length(Q) + 1); S[1] <- S_initial * capacity_live + (capacity - capacity_live)   
  R_rec <- vector("numeric",length(Q))                      
  for (yr in 1:nrow(Q_month_mat)) {
    for (month in 1:frq) {
      t_index <- (frq * (yr - 1)) + month   
      S_state <- which.min(abs(S_states - S[t_index]))
      Qx <- Q_month_mat[yr,month]
      R <- R_disc_x[R_policy[S_state,month]]
      R_rec[t_index] <- R
      if ( (S[t_index] - R + Qx) > capacity) {
        S[t_index + 1] <- capacity
      }else{
        if ( (S[t_index] - R + Qx) < (capacity - capacity_live)) {
          S[t_index + 1] <- (capacity - capacity_live)
          R_rec[t_index] <- S[t_index] + Qx - (capacity - capacity_live)
        }else{
          S[t_index + 1] <- S[t_index] - R + Qx
        }
      }
    }
  }
  
  R_policy <- (R_policy - 1) / R_disc
  S <- ts(S[2:length(S)],start = start(Q),frequency = frq)
  R_rec <- ts(R_rec, start = start(Q), frequency = frq)
  if(plot) {
    plot(S, ylab = "storage", ylim = c(0, capacity))
    plot(R_rec, ylab = "release/inflow", ylim = c(0, R_max))
    lines(Q, col="pink")
  }
  
  results <- list(R_policy, Bellman, S, R_rec, Q_disc)
  names(results) <- c("release_policy", "Bellman", "storage", "releases", "flow_disc")
  
  return(results)
}
