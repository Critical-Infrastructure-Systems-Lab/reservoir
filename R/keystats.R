#' @title Generate a range of key water supply reservoir variables
#' @description For quuickly analyzing a range of stats relating to inflows, outflows, storage dynamics and performance.
#' @param Q                  time series object representing the release time series of a reservoir.
#' @param target             numerical. The constant target water delivery.
#' @param capacity           numerical. The active capacity of the reservoir.
#' @return Returns a wide range of statistics relating to the dynamics and performance of the reservoir.
#' @examples keystats(resX$Q_Mm3, target = 50, capacity = resX$cap_Mm3)
#' @import stats
#' @importFrom moments skewness
#' @export
keystats <- function(Q, target, capacity) {
  
  Q_annual <- aggregate(Q)
  
  if(length(Q_annual) < 30) warning("More than 30 years' inflow data recommended for reasonable results.")
  
  sim <- simRes(Q, target = target, capacity = capacity, plot = FALSE)
  QminusR <- Q - target
  
  ## INFLOW-OUTFLOW STATS
  mu <- mean(Q_annual)
  sigma <- sd(Q_annual)
  Cs <- round(skewness(Q_annual), 3) 
  Cv <- round(sd(Q_annual)/mean(Q_annual), 3)
  ACF1 <- round(acf(Q_annual, plot=FALSE)[[1]][2], 3)
  DR <- round(target / mean(Q), 3)
  m <- round((1 - DR) / Cv, 3)
  Hrst <- round(Hurst(Q_annual), 3)
  
  inflow_stats <- c(mu, sigma, Cs, Cv, ACF1, Hrst, DR, m)
  names(inflow_stats) <- c("mAF", "sdAF", "CsAF", "CvAF", "ac1AF", "Hurst", "draft_ratio", "drift")
  
  ## STORAGE
  SR <- round(capacity / mean(Q_annual), 3)
  SS <- round(capacity / sd(Q_annual), 3)
  Peclet <- round((mean(QminusR) * capacity) / (2 * var(QminusR)), 3)
  CP_long <- round(critPeriod(sim$storage, report = "longest"), 3)  
  CP_short <- round(critPeriod(sim$storage, report = "shortest"), 3)
  CP_ave <- round(critPeriod(sim$storage, report = "average"), 3)
  spill_ratio <- round(mean(sim$spill) / mean(Q), 3)
  spill_freq_annual <- round(1 - sum(aggregate(sim$spill)==0)
                             / length(aggregate(sim$spill)), 3)
  spill_freq <- round(1 - sum(sim$spill==0)
                      / length(sim$spill), 3)
  
  storage_stats <- c(SR, SS, Peclet, CP_long, CP_short, CP_ave,
                     spill_ratio, spill_freq_annual, spill_freq)
  names(storage_stats) <- c("storage_ratio", "standardized_storage", "Peclet_number",
                            "CritPeriod_longest", "CritPeriod_shortest", "CritPeriod_average",
                            "spill_ratio", "spill_freq_annual", "spill_freq")
  
  ## PERFORMANCE STATS
  no_fail_yield <- round(yield(Q, capacity = capacity, reliability = 1, plot = FALSE)$Yield, 3)
  yield_99 <- round(yield(Q, capacity = capacity, reliability = 0.99, plot = FALSE)$Yield, 3)
  yield_95 <- round(yield(Q, capacity = capacity, reliability = 0.95, plot = FALSE)$Yield, 3)
  yield_90 <- round(yield(Q, capacity = capacity, reliability = 0.90, plot = FALSE)$Yield, 3)
  rrv_ <- rrv(sim$releases, target)
  
  performance_stats <- c(no_fail_yield, yield_99, yield_95, yield_90,
                         round(rrv_[1], 3), round(rrv_[2], 3),
                         round(rrv_[3], 3), round(rrv_[4], 3), round(rrv_[5], 3))
  
  results <- list(inflow_stats, storage_stats, performance_stats)
  names(results) <- c("Inflow_Outflow_stats", "Storage_stats", "Performance_stats")
  
  return(results)
}
