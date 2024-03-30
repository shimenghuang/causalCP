# This script contains some utility functions based on a fixed graph,
#.  a fairly minimal Markov blanket with multivariate Gaussian distribution.
# Specifically, the graph contains 5 nodes, X1 is a parent of X2, X1 and X2
#.  are parents of X5 (Y), and X3 is a child of X5, which has a co-parent X4.

# ---- theoretical quantities for the Gaussian MB graph G ----

#' Compute the theoretical mean vector based of graph G.
comp_joint_mean <- function(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5) {
  mu1 <- m1
  mu2 <- a12 * m1 + m2
  mu3 <- (b15 * a53 + a12 * b25 * a53) * m1 +
    b25 * a53 * m2 + a53 * m5 + a43 * m4 + m3
  mu4 <- m4
  mu5 <- (b15 + a12 * b25) * m1 + b25 * m2 + m5
  mu <- c(mu1, mu2, mu3, mu4, mu5)
  names(mu) <- c('X1', 'X2', 'X3', 'X4', 'X5')
  return(mu)
}

#' Compute the theoretical covariance matrix based of graph G.
comp_joint_cov <- function(a12, b15, b25, a53, a43, s1, s2, s3, s4, s5) {

  # row X1
  sig11 <- s1^2
  sig12 <- a12 * s1^2
  sig13 <- (b15 * a53 + a12 * b25 * a53) * s1^2
  sig14 <- 0
  sig15 <- (b15 + a12 * b25) * s1^2

  # row X2
  sig22 <- a12^2 * s1^2 + s2^2
  sig23 <- a12 * (b15 * a53 + a12 * b25 * a53) * s1^2 + b25 * a53 * s2^2
  sig24 <- 0
  sig25 <- a12 * (b15 + a12 * b25) * s1^2 + b25 * s2^2

  # row X3
  sig33 <- (b15 * a53 + a12 * b25 * a53)^2 * s1^2 +
    (b25 * a53)^2 * s2^2 + a53^2 * s5^2 + a43^2 * s4^2 + s3^2
  sig34 <- a43 * s4^2
  sig35 <- (b15 * a53 + a12 * b25 * a53) * (b15 + a12 * b25) * s1^2 +
    b25 * a53 * b25 * s2^2 +  a53 * s5^2

  # row X4
  sig44 <- s4^2
  sig45 <- 0

  # row X5 (Y)
  sig55 <- (b15 + a12 * b25)^2 * s1^2 + b25^2 * s2^2 + s5^2

  Sig <- matrix(c(sig11, sig12, sig13, sig14, sig15,
                      0, sig22, sig23, sig24, sig25,
                      0,     0, sig33, sig34, sig35,
                      0,     0,     0, sig44, sig45,
                      0,     0,     0,     0, sig55),
                byrow = TRUE, 5, 5)
  Sig[lower.tri(Sig)] <- t(Sig)[lower.tri(Sig)]
  rownames(Sig) <- c('X1', 'X2', 'X3', 'X4', 'X5')
  colnames(Sig) <- c('X1', 'X2', 'X3', 'X4', 'X5')
  return(Sig)
}

#' Compute the conditional distribution (mean and covariance) based on graph G.
#'
#' @param mu A vector of the joint mean.
#' @param Sig A matrix of the joint covariance.
#' @param dep_idx Indices of the dependent variables.
#' @param cond_idx Indices of the conditioning variables.
#' @param cond_val Values of the conditioning variables. It should have the
#'   same length as `cond_idx`.
#'
#' @return A list of two elements: the conditional mean `mu` and the conditional
#'   covariance matrix `Sig`.
comp_cond_dist <- function(mu, Sig, dep_idx, cond_idx, cond_val) {
  Sig11 <- Sig[dep_idx, dep_idx, drop=FALSE]
  Sig12 <- Sig[dep_idx, cond_idx, drop=FALSE]
  Sig22 <- Sig[cond_idx, cond_idx, drop=FALSE]
  cond_mu <- mu[dep_idx] + Sig12 %*% solve(Sig22, (cond_val-mu[cond_idx]))
  cond_Sig <- Sig11 - Sig12 %*% solve(Sig22, t(Sig12))
  return(list(mu = cond_mu,
              Sig = cond_Sig))
}

#' Compute the population OLS coefficient (with an intercept)
comp_pop_ols <- function(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
                         s1, s2, s3, s4, s5, S) {

  mu <- comp_joint_mean(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5)
  mu_S <- c(mu[S], 1) # add the intercept
  Sig <- comp_joint_cov(a12, b15, b25, a53, a43, s1, s2, s3, s4, s5)
  if (length(S) == 0) {
    Sig_S <- matrix(0, nrow = 1, ncol = 1)
  } else {
    Sig_S <- Sig[S,S]
    Sig_S <- cbind(rbind(Sig_S, 0), 0) # add the intercept
  }
  # Sig_S[length(S) + 1, length(S) + 1] <- 1
  Exx_S <- Sig_S + mu_S %*% t(mu_S)
  Exy_S <- c(Sig[S, 'X5'], 0) + mu_S * mu['X5']
  solve(Exx_S, Exy_S)
}

# ---- simulate data ----

#' Simulate data from graph G with multivariate Gaussian variables.
#'
#' @details Graph G is a dag where Y (X5) has two parent X1 and X2
#'   and one children X3. X2 is also a child of X1, and X3 is also
#'   a child of X4. This can also be done by using
#'   `mvtnorm::rmvnorm(n, mu, Sig)`.
sim_data <- function(n,
                     m1, m2, m3, m4, m5,
                     s1, s2, s3, s4, s5,
                     a12, b15, b25, a53, a43) {
  eps1 <- rnorm(n, m1, s1)
  eps2 <- rnorm(n, m2, s2)
  eps3 <- rnorm(n, m3, s3)
  eps4 <- rnorm(n, m4, s4)
  eps5 <- rnorm(n, m5, s5)

  X1 <- eps1
  X4 <- eps4
  X2 <- a12 * X1 + eps2
  X5 <- b15 * X1 + b25 * X2 + eps5 # X5 is Y
  X3 <- a53 * X5 + a43 * X4 + eps3

  X <- cbind(X1, X2, X3, X4, X5)
  colnames(X) <- c('X1', 'X2', 'X3', 'X4', 'X5')
  return(X)
}

#' Simulate data from graph G with mixed multivariate Gaussian variables.
#'
#' @details Slow implementation at the moment.
sim_data_mix <- function(m1_all, m2_all, m3_all, m4_all, m5_all,
                         s1_all, s2_all, s3_all, s4_all, s5_all,
                         a12_all, b15_all, b25_all, a53_all, a43_all,
                         seed_num) {
  set.seed(seed_num)
  tXY <-  mapply(\(m1, m2, m3, m4, m5,
                   s1, s2, s3, s4, s5,
                   a12, b15, b25, a53, a43) {
    sim_data(1,
             m1, m2, m3, m4, m5,
             s1, s2, s3, s4, s5,
             a12, b15, b25, a53, a43)
  },
  m1_all, m2_all, m3_all, m4_all, m5_all,
  s1_all, s2_all, s3_all, s4_all, s5_all,
  a12_all, b15_all, b25_all, a53_all, a43_all)
  return(t(tXY))
}

#' Generate parameters for the 5 node MB graph.
#'
#' @details A set of initial values for the causal parameters and non-causal
#'   parameters is first generated based on a given seed number. Then, at each
#'   change point based on its type, a (random sub)set of parameters that
#'   correspond to causal or non-causal change will randomly increase or
#'   decrease a chosen percentage. Note that a `true_ccp` refers to time points
#'   where the causal parameters (m5, s5, b15, and b25) change and a `true_nccp`
#'   refers to a time point where the non-causal parameters (m1 to m4, s1 to s4,
#'   a12, a53, and a43) change. There can be overlaps between the two and the
#'   overlaps are considered as CCPs by the definition in the paper.
#'
#' @param n A scalar of total number of time indices.
#' @param true_ccp A vector of integers (time indices) that at least one causal
#'   parameters should change.
#' @param true_nccp A vector of integers (time indices) that at least one
#'   non-causal parameters should change.
#' @param seed_num An integer as seed number for the random change direction.
#'   Ignored if `fixed_sign` is not 0.
#' @param nccp_mean_nc Number of non-causal mean params to be changed (up to 4).
#' @param nccp_mean_sc The scale of the change applied to each of the params
#'   that should have a change.
#' @param fixed_sign A scalar indicating the change direction, can be one of
#'   -1 (only decrease), 1 (only increase), or 0 (random sign). Default to be 0.
#'
#' @return A list of two data.frames. One contains the initial parameters and
#'   one row per change point. Another contains one row per time index
#'   corresponding to the parameters at the time index.
#'
rand_mb_params_step <- function(n, true_ccp, true_nccp, seed_num,
                                nccp_mean_nc = 4, nccp_mean_sc = 0.5,
                                nccp_sig_nc = 4, nccp_sig_sc = 0.5,
                                nccp_coef_nc = 3, nccp_coef_sc = 0.5,
                                ccp_mean_nc = 1, ccp_mean_sc = 0.5,
                                ccp_sig_nc = 1, ccp_sig_sc = 0.5,
                                ccp_coef_nc = 2, ccp_coef_sc = 0.5,
                                fixed_sign = 0, params_init = NULL) {
  if (!(fixed_sign %in% c(-1, 0, 1))) {
    stop("fixed_sign must be either -1, 1, or 0!")
  }
  set.seed(seed_num)
  get_signs <- function(len, fixed_sign) {
    if (fixed_sign) {
      signs <- rep(fixed_sign, len)
    } else {
      signs <- sample(c(-1, 1), len, replace = TRUE, prob = c(0.5, 0.5))
    }
    signs
  }
  cp_lst <- get_cp_types(true_ccp, true_nccp)
  cps <- sort(as.integer(names(cp_lst)))
  seg_lens <- get_seg_lens(cps, n)
  nccp_params_mean <- c("m1", "m2", "m3", "m4")
  nccp_params_sig <- c("s1", "s2", "s3", "s4")
  nccp_params_coef <- c("a12", "a53", "a43")
  ccp_params_mean <- c("m5")
  ccp_params_sig <- c("s5")
  ccp_params_coef <- c("b15", "b25")
  # initial values of the params
  if (is.null(params_init)) {
    params_all <- data.frame(m1 = NA, m2 = NA, m3 = NA, m4 = NA, m5 = NA,
                             s1 = NA, s2 = NA, s3 = NA, s4 = NA, s5 = NA,
                             a12 = NA, a53 = NA, a43 = NA, b15 = NA, b25 = NA)
    params_all[1,] <- 1
  } else {
    params_all <- params_init
  }
  # go through each cps
  for (ii in seq_along(cp_lst)) {
    new_row <- params_all[ii,]
    if ("NCCP" %in% cp_lst[[ii]]) {
      nccp_mean_change <- sample(nccp_params_mean, nccp_mean_nc, replace = FALSE)
      nccp_sig_change <- sample(nccp_params_sig, nccp_sig_nc, replace = FALSE)
      nccp_coef_change <- sample(nccp_params_coef, nccp_coef_nc, replace = FALSE)
      new_row[, nccp_mean_change] <- new_row[, nccp_mean_change] *
        (1 + get_signs(length(nccp_mean_change), fixed_sign) * nccp_mean_sc)
      new_row[, nccp_sig_change] <- new_row[, nccp_sig_change] *
        sqrt(1 + get_signs(length(nccp_sig_change), fixed_sign) * nccp_sig_sc)
      new_row[, nccp_coef_change] <- new_row[, nccp_coef_change] *
        (1 + get_signs(length(nccp_coef_change), fixed_sign) * nccp_coef_sc)
    }
    if ("CCP" %in% cp_lst[[ii]]) {
      ccp_mean_change <- sample(ccp_params_mean, ccp_mean_nc, replace = FALSE)
      ccp_sig_change <- sample(ccp_params_sig, ccp_sig_nc, replace = FALSE)
      ccp_coef_change <- sample(ccp_params_coef, ccp_coef_nc, replace = FALSE)
      new_row[, ccp_mean_change] <- new_row[, ccp_mean_change] *
        (1 + get_signs(length(ccp_mean_change), fixed_sign) * ccp_mean_sc)
      new_row[, ccp_sig_change] <- new_row[, ccp_sig_change] *
        sqrt(1 + get_signs(length(ccp_sig_change), fixed_sign) * ccp_sig_sc)
      new_row[, ccp_coef_change] <- new_row[, ccp_coef_change] *
        (1 + get_signs(length(ccp_coef_change), fixed_sign) * ccp_coef_sc)
    }
    params_all <- rbind(params_all, new_row)
  }
  rownames(params_all) <- c("1", as.character(cps))
  # repeat each row the length of the segment times
  params_df <- lapply(1:nrow(params_all), \(ii) {
    params_all[rep(ii, seg_lens[ii]), ]
  }) %>% dplyr::bind_rows()
  rownames(params_df) <- NULL
  list(params_all = params_all,
       params_df = params_df)
}

#' Initialize the parameter data.frame.
#'
#' @return A data.frame containing all initial values of the MB model and three
#'   additional columns: `idx`, `is_ccp`, and `is_nccp`.
init_mb_params <- function(init_val = 1) {
  params_step <- data.frame(m1 = NA, m2 = NA, m3 = NA, m4 = NA, m5 = NA,
                            s1 = NA, s2 = NA, s3 = NA, s4 = NA, s5 = NA,
                            a12 = NA, a53 = NA, a43 = NA, b15 = NA, b25 = NA)
  params_step[1,] <- init_val
  params_step$idx <- 1
  params_step$is_ccp <- FALSE
  params_step$is_nccp <- FALSE
  rownames(params_step) <- NULL
  return(params_step)
}

#' Obtain parameters at each time point.
#' 
#' @param n A scalar of total number of time indices.
#' @param params_step A data.frame containing the set of parameters, index
#' (`idx`), and change point types (`is_ccp` and `is_nccp`).
#' 
#' @export
get_full_mb_params <- function(n, params_step) {
  cps <- sort(params_step$idx)[-1] # remove the 1 which is not a change point
  seg_lens <- get_seg_lens(cps, n)
  # repeat each row the length of the segment times
  params_df <- lapply(1:nrow(params_step), \(ii) {
    params_step[rep(ii, seg_lens[ii]), ]
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-idx, -is_ccp, -is_nccp)
  rownames(params_df) <- NULL
  return(params_df)
}

