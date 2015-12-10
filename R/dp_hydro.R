#' @title Dynamic Programming for Hydropower Reservoir
#' @description Determines the optimal sequence of releases from the reservoir to minimise a penalty cost function based on water supply defict.
#' @param Q             vector or time series object. Net inflows to the reservoir.
#' @param capacity      numerical. The reservoir storage capacity (must be the same volumetric unit as Q).
#' @param capacity_live numerical. The reservoir live capacity (must be the same volumetric unit as Q).
#' @param surface_area  numerical. The reservoir surface area
#' @param hydro_cap     numerical. The hydropower plant electric capacity (MW)
#' @param head          numerical. The maximum hydraulic head of the hydropower plant (m)
#' @param qmax          numerical. The maximum flow into the hydropower plant
#' @param depth         numerical. The reservoir maximum depth
#' @param dep_vol_curve string. The relationship between depth and volume of reservoir.
#' @param input_curve   matrix of 3 columns: depth, area, volume
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @return Returns the time series of optimal releases and, if requested, the reliability, resilience and vulnerability of the system.
#' @references Loucks, D.P., van Beek, E., Stedinger, J.R., Dijkman, J.P.M. and Villars, M.T. (2005) Water resources systems planning and management: An introduction to methods, models and applications. Unesco publishing, Paris, France.
#' @examples \donttest{storage_cap <- 4 * mean(aggregate(ResX_inflow.ts)) # set storage ratio of 4 years
#' demand <- 0.8 * mean(ResX_inflow.ts) # set draft ratio of 0.8
#' optimal.releases <- dp(ResX_inflow.ts, capacity = storage_cap, target = demand)
#' }
#' @seealso \code{\link{sdp}} for Stochastic Dynamic Programming
#' @import stats 
#' @export
dp_hydro <- function(Q, capacity, capacity_live, surface_area, hydro_cap, head, qmax, depth,
                     dep_vol_curve = "l", input_curve, 
                     S_disc = 1000, R_disc = 10, S_initial = 1, plot = TRUE) {
  
  if (is.ts(Q) == FALSE && is.vector(Q) == FALSE) {
    stop("Q must be time series or vector object")
  }
  
  if ( dep_vol_curve == "n") { 
    N <- 2 * capacity / (depth * surface_area)
  } else if ( dep_vol_curve == "f") {
    f <- sqrt(2) / 3 * (surface_area*10^6) ^ (3/2) / (capacity*10^6)
  }
  
  if (!missing(input_curve)){
    colnames(input_curve) <- c("depth", "area", "volume")
  }
  
  S_states <- seq(from = capacity - capacity_live, to = capacity, by = capacity_live / S_disc)
  
  if ( !missing(qmax) ) {
    R_max <- qmax
  } else {
    R_max <- hydro_cap * 10^6 / (1000 * 9.81 * head) * (365/12 * 24 * 60 * 60) / 10^6 #10^6 m3/month
  }
  
  R_disc_x <- seq(from = 0, to = R_max, by = R_max / R_disc)
  State_mat <- matrix(0, nrow = length(S_states), ncol = length(R_disc_x))
  State_mat <- apply(State_mat, 2, "+", S_states)
  State_mat <- t(apply(State_mat, 1, "-", R_disc_x))
  Rev_to_go <- vector("numeric", length = length(S_states))
  Bellman <- matrix(0, nrow = length(S_states), ncol = length(Q))
  R_policy <- matrix(0, ncol = length(Q), nrow = length(S_states))
  
  # POLICY OPTIMIZATION----------------------------------------------------------------
  
  for (t in length(Q):1) {
    Release_mat <- matrix(rep(R_disc_x, length(S_states)),
                          ncol = length(R_disc_x), byrow = TRUE)
    Balance_mat <- State_mat + Q[t]
    Release_mat[which(Balance_mat < (capacity - capacity_live))] <- NaN
    Balance_mat[which(Balance_mat < (capacity - capacity_live))] <- NaN
    Balance_mat[which(Balance_mat > capacity)] <- capacity
    Implied_S_state <- round(1 + ( ((Balance_mat - (capacity - capacity_live)) / (capacity_live)) *
                                     (length(S_states) - 1)))
    if ( dep_vol_curve == "l"){
      H_mat <- (Balance_mat + S_states) / 2 / surface_area 
    } else if ( dep_vol_curve == "f") {
      H_mat <- ( 6 * (Balance_mat + S_states) * 10^6 / 2 / f^2 ) ^ (1/3)
    } else if ( dep_vol_curve == "n") {
      H_mat <- depth * ( (Balance_mat + S_states) / 2 / capacity ) ^ (N/2)
    } else if ( dep_vol_curve == "u") {
      H_mat <-  matrix(sapply( ((Balance_mat + S_states) * 10^6 / 2), 
                               function (x) 
                                 ifelse(is.nan(x), NaN, input_curve$depth[ which.min(abs( input_curve$volume - x )) ]) ), 
                       nrow = length(S_states), ncol = length(R_disc_x))
    }
    Rev_mat <- Release_mat * H_mat 
    Rev_mat2 <- Rev_mat + matrix(Rev_to_go[Implied_S_state],
                                 nrow = length(S_states))
    Rev_to_go <- apply(Rev_mat2, 1, max, na.rm = TRUE)
    Bellman[, t] <- Rev_to_go
    R_policy[, t] <- apply(Rev_mat2, 1, which.max)
  }
  # ===================================================================================
  
  # POLICY SIMULATION------------------------------------------------------------------
  
  S <- vector("numeric", length(Q) + 1)
  S[1] <- S_initial * capacity_live + (capacity - capacity_live)
  R <- vector("numeric", length(Q))
  for (t in 1:length(Q)) {
    S_state <- round(1 + ( ((S[t] - (capacity - capacity_live))/ capacity_live) *
                             (length(S_states) - 1)))
    R[t] <- R_disc_x[R_policy[S_state, t]]
    if ( (S[t] - R[t] + Q[t]) > capacity) {
      S[t + 1] <- capacity
    } else {
      S[t + 1] <- S[t] - R[t] + Q[t]
    }
  }
  S <- ts(S[2:length(S)], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  
  results <- list(S, R)
  names(results) <- c("storage", "releases")
  
  # ===================================================================================
  
  if (plot) {
    plot(S, ylab = "storage", ylim = c(0, capacity))
    plot(R, ylab = "release/inflow", ylim = c(0, R_max))
    lines(Q, col="pink")
  }
  return(results)
} 