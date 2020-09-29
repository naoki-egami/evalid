#' Partial Conjunction Test
#' @param pvalue a vector of one-sided p-values
#' @param threshold A threshold for the partial conjunction test. If NULL, we compute partial conjunction p-values for all thresholds. Default is NULL.
#' @return \code{pct} returns the following values.
#'  \itemize{
#'    \item \code{pc}: A table of partial conjunction p-values. \code{threshold} indicates each threshold. \code{p_value} shows its corresponding partial conjunction p-value. \code{h_num} indicates the corresponding hypothesis, where each hypothesis is indexed by the original order in which they appear in the argument \code{pvalue}.
#'    \item \code{tpate}: Estimates of the T-PATE
#'  }
#' @description \code{pct} implements the sign-generalization.
#' @references Egami and Hartman. (2020+). Elements of External Validity: Framework, Design, and Analysis
#' @export

pct <- function(pvalue, threshold = NULL){

  n <- length(pvalue)
  if(is.null(threshold) == TRUE){
    u  <- seq(from = 1, to = n)
  }else{
    u <- threshold
  }

  # sorting
  h_num <- seq(1:length(pvalue))[order(pvalue)]
  p <- sort(pvalue)

  # Holm
  pu <- (n - u + 1)*p[u]

  if(length(pu) == 1){
    pc <- pu
  }else{
    # monotonicity
    pc0 <- c()
    pc0[1] <- pu[1]
    for(i in 2:n){
      pc0[i] <- max(pc0[(i-1)], pu[i])
    }
    pc0[pc0 > 1] <- 1
    pc <- as.data.frame(cbind(u, pc0, h_num))
    colnames(pc) <- c("threshold", "p_value", "h_num")
  }
  return(pc)
}
