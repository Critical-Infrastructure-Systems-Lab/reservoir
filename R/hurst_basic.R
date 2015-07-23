#' @title Hurst coefficient estimation
#' @description Hurst coefficient estimation.
#' @param Q   vector or annualized time series object. Net inflows or streamflow totals.
#' @return Returns an estimate of the Hurst coefficient, H (0.5 < H < 1).
#' @references H.E.Hurst (1951) Long-term storage capacity of reservoirs, Transactions of the American Society of Civil Engineers 116, 770-808.
#' @references Pfaff, B. (2008) Analysis of integrated and cointegrated time series with R, Springer, New York. [p.68]
#' @examples Q_annual <- aggregate(ResX_inflow.ts) #convert monthly to annual data
#' Hurst(Q_annual)
#' @import stats
#' @export
Hurst <- function(Q) {
    if (frequency(Q) != 1)
        warning("Q is not an annualized series")
    if (is.ts(Q) == FALSE && is.vector(Q) == FALSE)
        stop("Q must be time series or vector object")
    Q <- Q - mean(Q)
    max.Q <- max(cumsum(Q))
    min.Q <- min(cumsum(Q))
    RS <- (max.Q - min.Q) / sd(Q)
    H <- log(RS) / log(length(Q))
    return(H)
}
