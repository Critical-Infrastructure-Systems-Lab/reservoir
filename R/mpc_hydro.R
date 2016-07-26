#' @title Simulate water supply reservoir operations informed by seasonal forecasts.
#' @description For simulating a water supply reservoir operated with rolling horizon, adaptive control (Model Predictive Control).
#' @param Q             time series object. Observed reservoir inflow totals. Recommended units: Mm^3 (Million cubic meters).
#' @param forecast      matrix: N * H, where N is the number of forecast-issue periods and H is the forecast horizon (i.e., number of periods) .
#' @param start_yr      the start year of the forecast. If the 'Q' and 'forecast' parameters have the same start year then leave blank.
#' @param capacity      numerical. The reservoir storage capacity. Recommended units: Mm^3 (Million cubic meters).
#' @param capacity_live numerical. The volume of usable water in the reservoir ("live capacity" or "active storage"). capacity_live <= capacity. Default capacity_live = capacity. Must be in Mm^3.
#' @param surface_area  numerical. The reservoir water surface area at maximum capacity. Recommended units: km^2 (square kilometers).
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap          time series object of equal length to Q, vector of length frequency(Q), or numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param installed_cap numerical. The hydropower plant electric capacity (MW).
#' @param efficiency    numerical. The hydropower plant efficiency. Default is 0.9, but, unless user specifies an efficiency, it will be automatically re-estimated if head and qmax are supplied.
#' @param head          numerical. The maximum hydraulic head of the hydropower plant (m). Can be omitted if qmax is supplied.
#' @param qmax          numerical. The maximum flow into the hydropower plant. Can be omitted and estimated if head is supplied. Must be in volumetric units of Mm^3.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @return Returns a list of reservoir variables as time series for the forecast period. Also returns penalty cost during operating period and cost savings relative to operations without forecasts.
#' @examples #
#' Q <- resX$Q_Mm3
#' forecastQ <- bootcast(Q, start_yr = 1980, H = 3, plot = FALSE)
#' @import stats
#' @export
simcast_hydro <- function(Q, forecast, start_yr, capacity, capacity_live = capacity,
                          surface_area, max_depth, evap, installed_cap, head, qmax,
                          efficiency = 0.9, S_disc = 1000, R_disc = 10,
                          Q_disc = c(0, 0.2375, 0.4750, 0.7125, 0.95, 1),
                          S_initial = 1, plot = TRUE){
  
  frq <- frequency(Q)
  if (is.ts(Q)==FALSE) stop("Q must be seasonal time series object with frequency of 12 or 4")
  if (frq != 12 && frq != 4) stop("Q must have frequency of 4 or 12")
  
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
  
  if ((missing(head) || is.na(head)) && (missing(qmax) || is.na(qmax))) {
    stop("You must enter a value for either head or qmax")
  }
  if (!missing(head) && !missing(qmax) && missing(efficiency) && !is.na(head) && !is.na(qmax)) {
    efficiency <- installed_cap / (9.81 * 1000 * head * (qmax / ((365.25/frq) * 24 * 60 * 60)))
    if (efficiency > 1) {
      warning("Check head, qmax and installed_cap: calculated efficiency exceeds 100 %")
    }
  }
  if (missing(head) || is.na(head)) {
    head <- installed_cap / (efficiency * 9.81 * 1000 * (qmax / ((365.25/frq) * 24 * 60 * 60)))
  }
  if (missing(qmax) || is.na(qmax)){
    qmax <- (installed_cap / (efficiency * 9.81 * 1000 * head)) * ((365.25/frq) * 24 * 60 * 60)
  }

  
  # SET UP (storage-depth-area relationships)----------------------------------------------------- 
  if (missing(max_depth) || is.na(max_depth)){
    c <- sqrt(2) / 3 * (surface_area * 10 ^ 6) ^ (3/2) / (capacity * 10 ^ 6)
    GetLevel <- function(c, V){
      y <- (6 * V / (c ^ 2)) ^ (1 / 3)
      return(y)
    }
    GetArea <- function(c, V){
      Ay <- (((3 * c * V) / (sqrt(2))) ^ (2 / 3))
      return(Ay)
    }
    yconst <- head - GetLevel(c, capacity * 10 ^ 6)
    if (yconst <0){
      capacity_live <- min(capacity_live, capacity - (-yconst) ^ 3 * c ^ 2 / 6 / 10 ^ 6)
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
    yconst <- head - max_depth
    if (yconst <0){
      capacity_live <- min(capacity_live, capacity - (-yconst / max_depth) ^ (2 / c) * capacity )
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
  x <- reservoir::sdp_hydro(Q = Qtr, capacity = capacity, capacity_live = capacity_live,
                            surface_area = surface_area, max_depth = max_depth, evap = evap_tr,
                            installed_cap = installed_cap, head = head, qmax = qmax,
                            efficiency = efficiency, S_disc = S_disc, Q_disc = Q_disc,
                            S_initial = 1, plot = FALSE, tol = 0.999, Markov = FALSE)
  
  # SIMULATE FORECAST-INFORMED MODEL
  message("Beginning simuation...")
  S <- vector("numeric", length(Qfc)); S[1] <- S_initial * capacity
  R <- vector("numeric", length(Qfc))
  E <- vector("numeric", length(Qfc))
  y <- vector("numeric", length(Qfc))
  Sp <- vector("numeric", length(Qfc))
  Power <- vector("numeric", length(Qfc))
  
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
      
      r2g <- x$Bellman[,matrix(rep(1:frq, frq + 1)[-seq(1, frq ^ 2, frq + 1)], ncol = frq)[H, month]]
      R[t_index] <- reservoir::dp_hydro(Q = ts(Q_forecast, frequency = frq), capacity = capacity,
                                        capacity_live = capacity_live, surface_area = surface_area,
                                        evap = evap_forecast, installed_cap = installed_cap,
                                        head = head, qmax = qmax, max_depth = max_depth,
                                        efficiency = efficiency, S_disc = S_disc, R_disc = R_disc,
                                        S_initial = Storage / capacity, r2g = r2g, plot = FALSE)$Release_Mm3[1]

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
      Power[t_index] <- max(efficiency * 1000 * 9.81 * (GetLevel(c,mean(S[t_index:(t_index + 1)]) * (10 ^ 6)) + yconst) * 
                              R[t_index] / (365.25 / frq * 24 * 60 * 60), 0)
    }
  }
  S <- ts(S[1:(length(S) - 1)],start = start(Qfc),frequency = frq)
  R <- ts(R, start = start(Qfc), frequency = frq)
  E <- ts(E, start = start(Q), frequency = frq)
  y <- ts(y, start = start(Q), frequency = frq)


  
  S <- ts(S[1:(length(S)-1)], start = start(Qfc), frequency = frq)
  R <- ts(R, start = start(Qfc), frequency = frq)
  E <- ts(E, start = start(Qfc), frequency = frq)
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Sp <- ts(Sp, start = start(Qfc), frequency = frq)
  Power <- ts(Power, start = start(Q), frequency = frq)
  Energy_MWh <- sum(Power * (365.25 / frq) * 24)

  
  ## SIMULATE THE RESERVOIR WITHOUT FORECAST-INFORMED OPERATION
  #xx <- simRes(Qfc, target = target, capacity = capacity,
  #             surface_area = surface_area, max_depth = max_depth,
  #             evap = evap_fc, plot = FALSE, S_initial = S_initial,
  #             policy = x)
  
  #total_penalty_sdp <- sum( ( (target - xx$releases) / target) ^ loss_exp)
  #
  #cost_saving <- 100 * (1 - (total_penalty / total_penalty_sdp))
  
  if(plot) {
    plot(R, ylab = "Turbined release [Mm3]", ylim = c(0, qmax), main = paste0("Total output = ", round(Energy_MWh/1000000, 3), " TWh"))
    plot(Power, ylab = "Power [MW]", ylim = c(0, installed_cap))
    plot(S, ylab = "Storage [Mm3]", ylim = c(0, capacity))
    plot(Sp, ylab = "Uncontrolled spill [Mm3]")
  }
  
  results <- list(S, R, E, y, Sp, Power, Energy_MWh)
  names(results) <- c("storage", "releases", "evap_loss", "water_level",
                      "spill", "power", "Energy_MWh")
  return(results)
}