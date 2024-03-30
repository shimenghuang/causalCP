# ---- Standard binary segmentation ----

#' Standard binary segmentation for causal change points.
#'
#' @details In the case that the estimation function (`esti_fun`) does not
#'   achieve global min at one of the true CCPs when multiple CCPs exist, the
#'   returned vector may contain time points that are not CCPs. One may use
#'   the `prune_candi` function to remove the NCCPs.
#'
#' @return A vector of estimated causal change points.
#' 
#' @export
stdBS_ccp <- function(X, Y, st, en, min_gap = 100,
                      test_fun = NULL,
                      test_param = list(alpha = 0.05),
                      esti_fun = NULL,
                      esti_param = list(intercept = TRUE),
                      verbose = TRUE) {
  message("st: ", st, ", en: ", en)
  if (en - st < min_gap) return(numeric())
  pval <- do.call(test_fun, c(list(X = X, Y = Y, I = st:en), test_param))
  if (pval > test_param$alpha) return(numeric())
  one_est <- do.call(esti_fun, c(list(X = X, Y = Y, I = st:en), esti_param))
  message("one_est: ", one_est)
  est_left <- stdBS_ccp(X, Y, st, one_est-1, min_gap, test_fun,
                        test_param, esti_fun, esti_param, verbose)
  est_right <- stdBS_ccp(X, Y, one_est, en, min_gap, test_fun,
                         test_param, esti_fun, esti_param, verbose)
  return(c(est_left, one_est, est_right))
}

# ---- Seeded binary segmentation ----

#' Calculate the boundaries of seeded intervals.
#'
#' @param n Total number of observations.
#' @param decay A real number in `(1, 2]` which is the inverse of the decay
#'   parameter in `[1/2, 1)`. Default value is `sqrt(2)`.
#' @param min_len Minimum length of the intervals. If NULL, all length are
#'   reserved.
#' @param return_unique A boolean, whether to return only the unique intervals
#'   or not. Default to be `TRUE`.
#'
#' @return A matrix with 2 columns the start and end indices of the seeded
#'   intervals.
#'
#' @references Solt Kovács, Housen Li, Peter Buhlmann, and Axel Munk (2022),
#'   Seeded binary segmentation: a general methodology for fast and optimal
#'   changepoint detection, Biometrika, Volume 110, Issue 1, Pages 249–256,
#'   https://doi.org/10.1093/biomet/asac052
#'   
#' @export
seeded_intervals <- function(n, decay = sqrt(2), min_len = NULL, return_unique = TRUE) {
  n	<- as.integer(n)
  depth	<- ceiling(log(n, base = decay)) # k
  boundary_df <- data.frame(st = 1, en = n, depth = 1)
  for (ii in 2:depth){
    len_interval <- n * (1/decay)^(ii-1) # l_k
    if (!is.null(len_interval) & len_interval < min_len) {
      break
    }
    num_interval <- 2 * ceiling(decay^(ii-1)) - 1 # n_k
    new_intervals <- data.frame(st = floor(seq(1, n-len_interval+1, length.out = num_interval)),
                                en = ceiling(seq(len_interval, n, length.out = num_interval)))
    new_intervals$depth <- ii
    boundary_df	<- rbind(boundary_df, new_intervals)
  }
  if (return_unique) {
    return(unique(boundary_df))
  }
  boundary_df
}

#' Remove intervals containing the given indices.
#'
#' @details This function removes the intervals (rows) in `boundary_df` which
#'   contain the indices in `est_ccp`.
#'
#' @param boundary_df A data.frame returned by `seeded_intervals`.
#' @param est_ccp A vector of indices, e.g., estimated causal change points.
#'
#' @return A data.frame with possibly reduced number of rows.
eliminate_intervals <- function(boundary_df, est_ccp) {
  if (nrow(boundary_df) == 0 | length(est_ccp) == 0) return(boundary_df)
  remove_flag <- sapply(est_ccp, \(ii) {
    sapply(1:nrow(boundary_df), \(rr) {
      (ii > boundary_df[rr,]$st) & (ii <=  boundary_df[rr,]$en)
    })
  })
  # Note: in case of bounadry_df having only one row, the returned remove_flag
  #   will be a vector of length of est_ccp, thus need to enforce the matrix
  remove_flag <- matrix(remove_flag, ncol = length(est_ccp))
  remove_flag <- rowSums(remove_flag) > 0
  boundary_df[!remove_flag,]
}

#' Seeded binary segmentation for causal change points.
#'
#' @details This procedure relies on a test function for whether an interval
#'   contains a change point or not and a estimation function for estimating
#'   the location of a change point.
#'
#' @param X The design matrix of dimension `n_tot x d`.
#' @param Y The response vector of length `n_tot`.
#' @param test_fun A function for testing the existence of a causal change point,
#'   taking values `X`, `Y`, an interval `I`, and possibly other parameters
#'   specified in `test_param`, returns a p-value between 0 and 1.
#' @param test_param A list containing additional parameters used by `test_fun`.
#' @param est_fun A loss function for estimating the location of a change point,
#'   taking values `X`, `Y`, an interval `I`, and possibly other parameters
#'   specified in `est_param`, returns a change point estimate.
#' @param est_param A list containing additional parameters used by `loss_fun`.
#'
#' @return A vector of estimated change points.
#' 
#' @export
seedBS_ccp <- function(X, Y, boundary_df, est_ccp = integer(),
                       test_fun = NULL,
                       test_param = list(alpha = 0.05),
                       esti_fun = NULL,
                       esti_param = list(intercept = TRUE),
                       verbose = TRUE) {
  if (verbose) message("nrow(boundary_df) = ", nrow(boundary_df))
  if (nrow(boundary_df) < 1) {
    return(est_ccp)
  }
  else {
    max_depth <- max(unique(boundary_df$depth))
    if (verbose) message("max_depth = ", max_depth)
    boundary_tmp <- boundary_df[boundary_df$depth == max_depth,]
    pvals <- rep(NA, nrow(boundary_tmp))
    for (jj in 1:nrow(boundary_tmp)) {
      st <- boundary_tmp$st[jj]
      en <- boundary_tmp$en[jj]
      pvals[jj] <- do.call(test_fun, c(list(X = X, Y = Y, I = st:en), test_param))
    }
    ord_pval <- order(pvals) # smallest to largest p-values
    boundary_tmp <- boundary_tmp[ord_pval,]
    boundary_tmp <- boundary_tmp[sort(pvals) < test_param$alpha,]
    while(nrow(boundary_tmp) > 0) {
      st <- boundary_tmp$st[1]
      en <- boundary_tmp$en[1]
      if (verbose) message("estimating in ", st, " to ", en)
      one_est <- do.call(esti_fun, c(list(X = X, Y = Y, I = st:en), esti_param))
      est_ccp <- c(est_ccp, one_est)
      boundary_tmp <- eliminate_intervals(boundary_tmp, one_est)
    }
    boundary_df <- eliminate_intervals(boundary_df, est_ccp)
    boundary_df <- boundary_df[!boundary_df$depth == max_depth,]
    seedBS_ccp(X, Y, boundary_df, est_ccp, test_fun, test_param, esti_fun, esti_param)
  }
}

# ---- Pruning function ----

#' Prune candidates and remove the ones that are not causal change points.
#'
#' @return A vector of estimated causal change points
#' 
#' @export
prune_candi <- function(X, Y, candi, PS, alpha, test_name, ...) {
  if (length(candi) == 0) return(candi)
  n_tot <- nrow(X)
  segs <- get_segs(candi, n_tot)
  max_pvals <- rep(NA, length(candi))
  for (ii in 1:(length(segs)-1)) {
    max_pvals[ii] <- max(two_env_test(X, Y, segs[[ii]], segs[[ii+1]], PS, test_name, ...)$pvals)
  }
  return(candi[max_pvals < alpha])
}
