#' reservoir: Tools for Analysis, Design, and Operation of Water Supply Storages
#'
#' Measure single reservoir performance using resilience, reliability, and vulnerability metrics; compute storage-yield-reliability relationships; determine no-fail Rippl storage with sequent peak analysis; optimize release decisions using determinisitc and stochastic dynamic programming; evaluate inflow characteristics.
#' 
#' @section Analysis and design functions:
#' The \code{\link{Rippl}} function executes the sequent peak algorithm to determine the no-fail storage for given inflow and release time series.
#' The \code{\link{storage}} function gives the design storage for a specified time-based reliability and yield. Similarly, the \code{\link{yield}} function computes yield given the storage capacity.
#' The \code{\link{simRes}} function simulates a reservoir under standard operating policy, or using an optimised policy produced by \code{\link{sdp_supply}}.
#' The \code{\link{rrv}} function returns three reliability measures, relilience, and dimensionless vulnerability for given storage, inflow time series, and target release. Users can assume Standard Operating Policy, or can apply the output of \code{\link{sdp_supply}} to determine the RRV metrics under different operating objectives.
#' The \code{\link{Hurst}} function estimates the Hurst coefficient for an annualized inflow time series.   
#' @section Optimization functions:
#' Dynamic Programming functions find the optimal sequence of releases for a given reservoir. Stochastic Dynamic Programming  functions find the optimal release policy for a given reservoir, based on storage, within-year time period and, optionally, current-period inflow. 
#' For single-objective water supply reservoirs, users may specify a loss exponent parameter for supply deficits and then optimize reservoir release decisions using Dynamic Programming (\code{\link{dp_supply}}) or Stochastic Dynamic Programming (\code{\link{sdp_supply}}). There is also an option to simulate the output of \code{\link{sdp_supply}} using the \code{\link{rrv}} function to validate the policy under alternative inflows or analyze reservoir performance under different operating objectives.
#' The optimal operating policy for hydropower operations can be found using Dynamic Programming \code{\link{dp_hydro}} or Stochastic Dynamic Programming \code{\link{sdp_hydro}}. The target is to maximise total energy output over the length of the input time series of inflows.
#' @section Storage-depth-area relationships:
#' All reservoir analysis and optimization functions, with the exception of \code{\link{Rippl}}, allow the user to incorporate evaporation losses from the reservoir surface. There are two relationships. The simplest is based on the reservoir \code{capacity} and \code{surface_area} variables. A more nuanced relationship is implemeted if the user also provides the \code{max_depth} input variable.
#' Users must use the recommended units when implementing evaporation losses.
#' @docType package
#' @name reservoir
#' @examples \donttest{# 1. Express the distribution of Rippl storage for a known inflow process...
#' 
#' # a) Assume the inflow process follows a lognormal distribution
#' # (meanlog = 0, sdlog = 1):
#' x <- rlnorm(1200)
#' 
#' # b) Convert to a 100-year, monthly time series object beginning Jan 1900
#' x <- ts(x, start = c(1900, 1), frequency = 12)
#' 
#' # c) Begin reservoir analysis... e.g., compute the Rippl storage
#' x_Rippl <- Rippl(x, target = mean(x) * 0.9)
#' no_fail_storage <- x_Rippl$Rippl_storage
#' 
#' # d) Resample x and loop the procedure multiple times to get the
#' # distribution of no-failure storage for the inflow process assuming 
#' # constant release (R) equal to 90 percent of the mean inflow.
#' no_fail_storage <- vector("numeric", 500)
#' for (i in 1:length(no_fail_storage)){
#'   x <- ts(rlnorm(1200), start = c(1900, 1), frequency = 12)
#'   no_fail_storage[i] <- Rippl(x, target = mean(x) * 0.9 ,plot = FALSE)$Rippl_storage
#' }
#' hist(no_fail_storage)
#' 
#' 
#' # 2. Trade off between annual reliability and vulnerability for a given system...
#' 
#' # a) Define the system: inflow time series, storage, and target release.
#' inflow_ts <- resX$Q_Mm3
#' storage_cap <- resX$cap_Mm3
#' demand <- 0.3 * mean(resX$Q_Mm3)
#' 
#' # b) define range of loss exponents to control preference of high reliability
#' # (low loss exponent) or low vulnerability (high loss exponent).
#' loss_exponents <- c(1.0, 1.2, 1.4, 1.6, 1.8, 2)
#' 
#' # c) set up results table
#' pareto_results <- data.frame(matrix(ncol = 2, nrow = length(loss_exponents)))
#' names(pareto_results) <- c("reliability", "vulnerability")
#' row.names(pareto_results) <- loss_exponents
#' 
#' # d) loop the sdp function through all loss exponents and plot results
#' for (loss_f in loss_exponents){
#'  sdp_temp <- sdp_supply(inflow_ts, capacity = storage_cap, target = demand, rep_rrv = TRUE,
#'  S_disc = 200, R_disc = 20, plot = FALSE, loss_exp = loss_f, Markov = TRUE)
#'  pareto_results$reliability[which(row.names(pareto_results)==loss_f)] <- sdp_temp$annual_reliability
#'  pareto_results$vulnerability[which(row.names(pareto_results)==loss_f)] <- sdp_temp$vulnerability
#'  }
#' plot (pareto_results$reliability,pareto_results$vulnerability, type = "b", lty = 3)
#' }

NULL