#' @title Simulate water supply reservoir operations informed by seasonal forecasts.
#' @description For simulating a water supply reservoir operated with rolling horizon, adaptive control (Model Predictive Control).
#' @param Q             time series object. Observed reservoir inflow totals. Recommended units: Mm^3 (Million cubic meters).
#' @param forecast      matrix: N * H, where N is the number of forecast-issue periods and H is the forecast horizon (i.e., number of periods) .
#' @param start_yr      the start year of the forecast. If the 'Q' and 'forecast' parameters have the same start year then leave blank.
#' @param capacity      numerical. The reservoir storage capacity. Recommended units: Mm^3 (Million cubic meters).
#' @param target        numerical. The target release constant. Recommended units: Mm^3 (Million cubic meters).
#' @param surface_area  numerical. The reservoir water surface area at maximum capacity. Recommended units: km^2 (square kilometers).
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap          time series object of equal length to Q, vector of length frequency(Q), or numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param loss_exp      numeric. The exponent of the penalty cost function--i.e., Cost[t] <- ((target - release[t]) / target) ^ **loss_exp**). Default value is 2.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @return Returns a list of reservoir variables as time series for the forecast period. Also returns penalty cost during operating period and cost savings relative to operations without forecasts.
#' @examples Q <- resX$Q_Mm3
#' forecastQ <- bootcast(Q, start_yr = 1980, H = 3, plot = FALSE)
#' layout(1:3)
#' simQ <- simcast_supply(Q, resX$cap_Mm3, target = 0.3*mean(Q),
#' forecast = forecastQ, start_yr=1980, S_disc = 200)
#' @import stats
#' @importFrom graphics abline lines
#' @export
simcast_supply <- function(Q, forecast, start_yr, capacity, target, surface_area,
                               max_depth, evap, S_disc = 1000, R_disc = 10,
                               Q_disc = c(0, 0.2375, 0.4750, 0.7125, 0.95, 1),
                               loss_exp = 2, S_initial = 1, plot = TRUE){ # add rep_rrv
  
  frq <- frequency(Q)
  if (is.matrix(forecast)==TRUE){
    H <- H <- ncol(forecast)
  } else {
    H <- 1
  }
  if (start(Q)[2] != 1){
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
  }
  if(end(Q)[2] != frq){
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
  }
  if (missing(evap)) {
    evap <- ts(rep(0, length(Q)), start = start(Q), frequency = frq)
  }
  if(length(evap) == 1) {
    evap <- ts(rep(evap, length(Q)), start = start(Q), frequency = frq)
  }
  if (length(evap) != length(Q) && length(evap) != frq){
    stop("Evaporation must be either a time series of length Q, a vector of length frequency(Q), or a single numeric constant")
  }
  if (length(evap) == frq){
    evap <- ts(rep(evap, length(Q) / frq), start = start(Q), frequency = frq)
  } else {
    if(is.ts(evap)==FALSE) stop("Evaporation must be either a time series of length Q or a vector of length frequency(Q) for a seasonal evaporation profile")
    evap <- window(evap, start = start(Q), end = end(Q), frequency = frq)
  }
  if (missing(surface_area)) {
    surface_area <- 0
  }
  
  # SET UP (storage-depth-area relationships)----------------------------------------------------- 
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
  
  # SPLIT TIME SERIES INTO TRAINING AND FORECAST PERIODS IF START YEAR IS SPECIFIED
  if (missing(start_yr)){
    Qfc <- Q
    Qtr <- Q
    evap_fc <- evap
    evap_tr <- evap
  } else {
    Qfc <- window(Q, start = c(start_yr, 1), frequency = frq)
    Qtr <- window(Q, end = c(start_yr - 1, frq), frequency = frq)
    evap_fc <- window(evap, start = c(start_yr, 1), frequency = frq)
    evap_tr <- window(evap, end = c(start_yr - 1, frq), frequency = frq)
  }

  # GET COST-TO-GO FUNCTION FOR THE RESERVOIR
  message("Deriving the reservoir's cost-to-go function...")
  x <- reservoir::sdp_supply(Q = Qtr, capacity = capacity, target = target,
                               surface_area = surface_area, max_depth = max_depth, evap = evap_tr,
                               S_disc = S_disc, Q_disc = Q_disc, loss_exp = loss_exp,
                               S_initial = 1, plot = FALSE, tol = 0.999, Markov = FALSE)
  
  # SIMULATE FORECAST-INFORMED MODEL
  message("Beginning simuation...")
  S <- vector("numeric", length(Qfc)); S[1] <- S_initial * capacity
  R <- vector("numeric", length(Qfc))
  E <- vector("numeric", length(Qfc))
  y <- vector("numeric", length(Qfc))
  Sp <- vector("numeric", length(Qfc))
  Q_month_mat <- matrix(Qfc, byrow = TRUE, ncol = frq)
  evap_seas <- rep(as.vector(tapply(evap, cycle(evap), FUN = mean)), (length(Qfc) / 12) + 1)
  for(yr in 1:nrow(Q_month_mat)){
    for (month in 1:frq){
      t_index <- (frq * (yr - 1)) + month
      Storage <- S[t_index]
      if(H > 1){
        Q_forecast <- forecast[t_index,]
      } else {
        Q_forecast <- forecast[t_index]
      }
      evap_forecast <- evap_seas[t_index:(t_index + H - 1)]
      c2g <- x$Bellman[,matrix(rep(1:frq, frq + 1)[-seq(1, frq ^ 2, frq + 1)], ncol = frq)[H, month]]
      R[t_index] <- reservoir::dp_supply(Q = Q_forecast, capacity = capacity,
                                         target = target, surface_area = surface_area,
                                         max_depth = max_depth, evap = evap_forecast,
                                         S_disc = S_disc, R_disc = R_disc, loss_exp = loss_exp,
                                         S_initial = Storage / capacity, c2g = c2g,
                                         plot = FALSE, rep_rrv = FALSE)$releases[1]
      E[t_index] <- GetEvap(s = S[t_index], q = Qfc[t_index], r = R, ev = evap_fc[t_index])
      y[t_index] <- GetLevel(c, S[t_index] * 10 ^ 6)
      if((S[t_index] - R[t_index] + Qfc[t_index] - E[t_index]) > capacity){
        S[t_index + 1] <- capacity
        Sp[t_index] <- S[t_index] - R[t_index] + Qfc[t_index] - E[t_index] - capacity
      }else{
        if((S[t_index] - R[t_index] + Qfc[t_index] - E[t_index]) < 0){
          S[t_index + 1] <- 0
          R[t_index] <- S[t_index] + Qfc[t_index] - E[t_index]
        }else{
          S[t_index + 1] <- S[t_index] - R[t_index] + Qfc[t_index] - E[t_index]
        }
      }
    }
  }
  S <- ts(S[1:(length(S)-1)], start = start(Qfc), frequency = frq)
  R <- ts(R, start = start(Qfc), frequency = frq)
  E <- ts(E, start = start(Qfc), frequency = frq)
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Sp <- ts(Sp, start = start(Qfc), frequency = frq)
  total_penalty <- sum( ( (target - R) / target) ^ loss_exp)
  
  ## SIMULATE THE RESERVOIR WITHOUT FORECAST-INFORMED OPERATION
  xx <- simRes(Qfc, target = target, capacity = capacity,
               surface_area = surface_area, max_depth = max_depth,
               evap = evap_fc, plot = FALSE, S_initial = S_initial,
               policy = x)
  
  total_penalty_sdp <- sum( ( (target - xx$releases) / target) ^ loss_exp)
  
  cost_saving <- 100 * (1 - (total_penalty / total_penalty_sdp))
  
  if(plot) {
    plot(R, ylab = "Controlled release", ylim = c(0, target))
    lines(xx$releases, col = "grey", lty = 2)
    plot(S, ylab = "Storage", ylim = c(0, capacity))
    lines(xx$storage, col = "grey", lty = 2)
    plot(Sp, ylab = "Uncontrolled spill")
    lines(xx$spill, col = "grey", lty = 2)
  }
  
  results <- list(S, R, E, y, Sp, total_penalty, cost_saving)
  names(results) <- c("storage", "releases", "evap_loss", "water_level",
                      "spill", "total_penalty", "percent_cost_saving")
  return(results)
}