#' @title Dynamic Programming for Hydropower Reservoir
#' @description Determines the optimal sequence of releases from the reservoir to minimise a penalty cost function based on water supply defict.
#' @param Q             time series object. Net inflows to the reservoir. Must be in volumetric units of Mm^3.
#' @param capacity      numerical. The total reservoir storage capacity (including unusable "dead" storage). Must be in Mm^3.
#' @param capacity_live numerical. The volume of usable water in the reservoir ("live capacity" or "active storage"). capacity_live <= capacity. Default capacity_live = capacity. Must be in Mm^3.
#' @param surface_area  numerical. The reservoir surface area at full capacity. Must be in square kilometers (km^2), or Mm^2.
#' @param installed_cap numerical. The hydropower plant electric capacity (MW).
#' @param efficiency    numerical. The hydropower plant efficiency. Default = 0.9.
#' @param head          numerical. The maximum hydraulic head of the hydropower plant (m).
#' @param qmax          numerical. The maximum flow into the hydropower plant
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
dp_hydro <- function(Q, capacity, capacity_live = capacity, surface_area, installed_cap, head, qmax,
                     efficiency = 0.9, S_disc = 1000, R_disc = 10, S_initial = 1, plot = TRUE) {
  
  
  if (is.ts(Q) == FALSE) {
    stop("Q must be time series object")
  }
  
  if (missing(head) && missing(qmax)) {
    stop("You must enter a value for either head or qmax")
  }
  
  if (!missing(head) && !missing(qmax) && missing(efficiency)) {
    efficiency <- installed_cap / (9.81 * 1000 * head * (qmax / ((365.25/frequency(Q)) * 24 * 60 * 60)))
    if (efficiency > 1) {
      warning("Check head, qmax and installed_cap: calculated efficiency exceeds 100 %")
    }
  }

  if (missing(head)) {
    head <- installed_cap / (efficiency * 9.81 * 1000 * (qmax / ((365.25/frequency(Q)) * 24 * 60 * 60)))
  }

  if (missing(qmax)) {
    qmax <- (installed_cap / (efficiency * 9.81 * 1000 * head)) * ((365.25/frequency(Q)) * 24 * 60 * 60)
  }
  
    f <- sqrt(2) / 3 * (surface_area * 10 ^ 6) ^ (3 / 2) / (capacity * 10 ^ 6)
    GetLevel <- function(f, V){
      y <- (6 * (V * 1000000) / (f ^ 2)) ^ (1 / 3)
      return(y)
    }
    GetArea <- function(f, V){
      Ay <- 0.5 * (6 * V * f) ^ (2 / 3)
      return(Ay)
    }
    
  ymax <- GetLevel(f, capacity)
  yconst <- head - ymax
  
  S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)
  R_disc_x <- seq(from = 0, to = qmax, by = qmax / R_disc)
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
    
    Release_mat[which(is.nan(Balance_mat[,1]))] <- 0           ## Correction to allow zero release under high evap
    Balance_mat[which(is.nan(Balance_mat[,1]))] <- capacity - capacity_live        ## Correction to allow zero release under high evap
    
    Balance_mat[which(Balance_mat > capacity)] <- capacity
    Implied_S_state <- round(1 + ((Balance_mat / capacity) * (length(S_states) - 1)))
    
    H_mat <- GetLevel(f, (Balance_mat + S_states) / 2) + yconst
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
  S[1] <- S_initial * capacity
  R <- vector("numeric", length(Q))
  Power <- vector("numeric", length(Q))
  for (t in 1:length(Q)) {
    S_state <- round(1 + ( (S[t] / capacity) * (length(S_states) - 1)))
    R[t] <- R_disc_x[R_policy[S_state, t]]
    if ( (S[t] - R[t] + Q[t]) > capacity) {
      S[t + 1] <- capacity
    } else {
      S[t + 1] <- max(0, S[t] - R[t] + Q[t])
    }
    Power[t] <- efficiency * 1000 * 9.81 * (GetLevel(f, mean(S[t:t+1])) + yconst) * R[t] / (365.25 / frequency(Q) * 24 * 60 * 60)
    
  }

  S <- ts(S[2:length(S)], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  Power <- ts(Power, start = start(Q), frequency = frequency(Q))
  Energy_MWh <- sum(Power * (365.25 / 12) * 24)
  results <- list(S, R, Power, Bellman, R_policy, Energy_MWh)
  names(results) <- c("Storage_Mm^3", "Release_Mm^3", "Power_MW", "Bellman", "R_policy", "Total_Energy_MWh")
  
  # ===================================================================================
  
  if (plot) {
    plot(S, ylab = "Storage", ylim = c(0, capacity))
    plot(R, ylab = "Release through turbine", ylim = c(0, qmax))
    plot(Power, ylab = "Power [MW]", ylim = c(0, installed_cap), main = paste0("Total output = ", round(Energy_MWh/1000000, 3), " TWh"))
  }
  return(results)
} 