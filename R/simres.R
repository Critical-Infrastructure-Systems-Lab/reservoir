#' @title Simulate a water supply reservoir with Standard Operating Policy
#' @description Simulates a reservoir for a given time series and assuming Standard Operating Policy (meet target at all times, unless constrained by available water in reservoir plus incoming flows). 
#' @param Q             vector or time series object. Net inflow totals to the reservoir. Mm^3 (Million cubic meters).
#' @param target        numerical constant, or a time series or vector of the target releases. Must be the same length as Q is given as a vector or time series. Mm^3 (Million cubic meters).
#' @param capacity      numerical. The reservoir capacity. Should be same volumetric unit as Q. Mm^3 (Million cubic meters).
#' @param surface_area  numerical. The reservoir surface area at full capacity. Must be in square kilometers (km^2), or Mm^2.
#' @param max_depth     numerical. The maximum water depth of the reservoir at the dam at maximum capacity. Must be in meters. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only.
#' @param evap          vector or time series object of length Q, or a numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param double_cycle  logical. If TRUE the Q and R time series will be replicated and placed end-to-end to double the simulation. Recommended if the critical period occurs at the end of the sequence.
#' @param plot          logical. If TRUE (the default) the storage and release time series are plotted.
#' @param policy        list. The output of the SDP function. If omitted, Standard Operating Policy is assumed.
#' @return Returns the no-fail storage capacity and corresponding storage behaviour time series.
#' @examples # simulate a reservoir assuming standard operating policy.
#' 
#' @import stats
#' @export
simSOP <- function(Q, target, capacity, surface_area, max_depth, evap,
                   double_cycle = FALSE, plot = TRUE, S_initial = 1, policy) {
  
  
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
  
  if (length(target) == 1){
    target <- rep(target, length(Q))
  }
  
  S <- vector("numeric",length(Q) + 1); S[1] <- S_initial * capacity    
  R <- vector("numeric",length(Q))
  E <- vector("numeric", length(Q))
  y <- vector("numeric", length(Q))
  Spill <- vector("numeric", length(Q))
  
  for (t in 1:length(Q)) {
    
    R[t] <- target[t]
    #E[t] <- GetArea(c, S[t] * 10 ^ 6) * evap[t] / 10 ^ 6
    E[t] <- GetEvap(s = S[t], q = Q[t], r = R[t], ev = evap[t])
    y[t] <- GetLevel(c, S[t] * 10 ^ 6)
    
    
    if ( (S[t] - R[t] + Q[t] - E[t]) > capacity) {
      S[t + 1] <- capacity
      Spill[t] <- S[t] - R[t] + Q[t] - capacity - E[t]
    } else {
      if (S[t] - R[t] + Q[t] - E[t] < 0) {
        S[t + 1] <- 0
        R[t] <- max(0, S[t] + Q[t] - E[t])
      } else {
        S[t + 1] <- S[t] +  Q[t] - R[t] - E[t]
      }
    }
  }
  
  S <- ts(S[1:length(S) - 1], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  E <- ts(E, start = start(Q), frequency = frequency(Q))
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Spill <- ts(Spill, start = start(Q), frequency = frequency(Q))
  
  results <- list(S, R, E, y, Spill)
  names(results) <- c("Storage", "Release", "Evaporation", "Water_level", "Spill")
  
  
  if (plot) {
    plot(R, ylab = "Release", ylim = c(0, max(target)))
    plot(S, ylab = "Storage", ylim = c(0, capacity))
    plot(Spill, ylab = "Spill")
  }
  
  return(results)

}













