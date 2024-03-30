## Figure 8 ##

# ---- setup on computerome ----
.libPaths(c("~/causalcp_proj/lib", .libPaths()))

# set up dependencies
if (!require("renv")) {
  install.packages("renv", "~/causalcp_proj/lib", repos = "http://cran.us.r-project.org")
}
renv::restore()
if (!require("dHSIC")) {
  install.packages("dHSIC", "~/causalcp_proj/lib", repos = "http://cran.us.r-project.org")
}

# ---- libraries ----
pkg_dir <- getwd()
pkg_dir <- "~/causalcp_proj/causalCPD/causalCPD/"
devtools::install_local(pkg_dir) # Note: only need to install once
devtools::load_all(pkg_dir)

# set up parallel
library(future.apply)
plan(multisession)

# ---- helper functions ----
test_wrap <- function(X, Y, I, PS = get_powerset(colnames(X)),
                      rel_grid = c(1/2), test_name = "chow1",
                      intercept = TRUE, rev_one_rest = FALSE, two_dir = FALSE,
                      save_resids = FALSE, alpha = 0.05) {
  res <- one_vs_rest_test(X[I,,drop=FALSE], Y[I], PS, rel_grid, test_name,
                          intercept, rev_one_rest, two_dir, save_resids)
  min(sapply(res$pvals_all, max)) * length(rel_grid)
}

esti_fun <- function(X, Y, I,
                     eval_gap = 1,
                     end_buff = max(ceiling(length(I) * 0.1), 10),
                     cond_sets = get_powerset(colnames(X))) {
  # Note: for cc_loss t is relative to the interval
  eval_grid <- seq(end_buff, length(I)-end_buff, by = eval_gap)
  vals <- sapply(eval_grid, \(t) cc_loss(X, Y, t, min(I), max(I),
                                         min_seg = 0.1 * length(I),
                                         cond_sets = cond_sets))
  return(min(I) + eval_grid[which.min(vals)])
}

# ---- main script ----
arg <- commandArgs(trailingOnly = TRUE) # not used here
seed_num <- 1000
reps <- 1:200
ns <- c(1000, 2000, 4000)
rel_true_ccp <- c(0.2, 0.8)
rel_true_nccp <- c(0.5)
rel_min_dis <- 0.2
PS <- get_powerset(as.character(1:4))
test_param <- list(test_name = "chow1", PS = PS, alpha = 0.05) # for actual test
esti_param <- list(eval_gap = 1, cond_sets = PS)

future_lapply(reps, \(kk) {
  devtools::load_all(pkg_dir)
  seedbs_all <- vector("list", length(ns))
  seedbs_prune_all <- vector("list", length(ns))
  stdbs_all <- vector("list", length(ns))
  stdbs_prune_all <- vector("list", length(ns))
  bp1_all <- vector("list", length(ns))
  bp2_all <- vector("list", length(ns))
  bp3_all <- vector("list", length(ns))
  seedbs_time <- vector("numeric", length(ns))
  seedbs_prune_time <- vector("numeric", length(ns))
  stdbs_time <- vector("numeric", length(ns))
  stdbs_prune_time <- vector("numeric", length(ns))
  bp1_time <- vector("numeric", length(ns))
  bp2_time <- vector("numeric", length(ns))
  bp3_time <- vector("numeric", length(ns))
  for (ll in seq_along(ns)) {
    n <- ns[ll]
    min_dis <- rel_min_dis * n
    true_ccp <- n * rel_true_ccp + 1
    true_nccp <- n * rel_true_nccp + 1
    true_rcp <- sort(c(true_ccp, true_nccp))
    
    params_init <- init_mb_params()
    params_step <- params_init[rep(1, 1 + length(true_ccp) + length(true_nccp)),]
    rownames(params_step) <- NULL
    params_step$idx <- c(1, sort(unique(c(true_ccp, true_nccp))))
    for (ii in params_step$idx) {
      if (ii %in% true_ccp) params_step[params_step$idx == ii, 'is_ccp'] <- TRUE
      if (ii %in% true_nccp) params_step[params_step$idx == ii, 'is_nccp'] <- TRUE
    }
    params_step[2:4, 'b15'] <- 1
    params_step[2:4, 'b25'] <- 0
    params_step[3:4, 'a53'] <- 2
    params_step[3:4, 'm3'] <- 2
    params_step[3:4, 'm4'] <- 2
    params_step[3:4, 's3'] <- 1
    params_step[4, 'b15'] <- 0
    params_step[4, 'b25'] <- 1
    params_step[4, 's5'] <- 2
    
    # # verify that the CCP are indeed CCPs (calculate the population OLS coef for all S in PS)
    # ps_tmp <- lapply(PS, \(S) {
    #   if (length(S) != 0) {
    #     paste0("X", S)
    #   }
    #   else {
    #     S
    #   }
    # })
    # beta_diff <- rep(NA, length(ps_tmp))
    # for (ss in seq_along(ps_tmp)) {
    #   S <- ps_tmp[[ss]]
    #   message("**", S, "**")
    #   beta1 <- with(params_step['1',], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
    #                                                 s1, s2, s3, s4, s5, S))
    #   message(beta1)
    #   beta2 <- with(params_step['2',], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
    #                                                 s1, s2, s3, s4, s5, S))
    #   message(beta2)
    #   beta_diff[ss] <- any(beta1 != beta2)
    # }
    # beta_diff <- rep(NA, length(ps_tmp))
    # for (ss in seq_along(ps_tmp)) {
    #   S <- ps_tmp[[ss]]
    #   message("**", S, "**")
    #   beta1 <- with(params_step['3',], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
    #                                                 s1, s2, s3, s4, s5, S))
    #   message(beta1)
    #   beta2 <- with(params_step['4',], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
    #                                                 s1, s2, s3, s4, s5, S))
    #   message(beta2)
    #   beta_diff[ss] <- any(beta1 != beta2)
    # }
    
    params_df <- get_full_mb_params(n, params_step)
    dat <- with(params_df,
                sim_data_mix(m1_all = m1, m2_all = m2, m3_all = m3, m4_all = m4, m5_all = m5,
                             s1_all = s1, s2_all = s2, s3_all = s3, s4_all = s4, s5_all = s5,
                             a12_all = a12, a53_all = a53, a43_all = a43,
                             b15_all = b15, b25_all = b25,
                             seed_num = kk))
    colnames(dat) <- as.character(1:5)
    X <- dat[,1:4]
    Y <- dat[,5]
    
    # seededBS
    cur_time <- proc.time()[[1]]
    boundary_df <- seeded_intervals(n, min_len = min_dis)
    seed_est <- seedBS_ccp(X, Y, boundary_df, est_ccp = c(),
                           test_fun = test_wrap, test_param = test_param,
                           esti_fun = esti_fun, esti_param = esti_param,
                           verbose = TRUE)
    seedbs_all[[ll]] <- sort(seed_est)
    seedbs_time[ll] <- proc.time()[[1]] - cur_time
    
    # seededBS with pruning
    seedbs_prune_all[[ll]] <- prune_candi(X, Y, sort(seed_est), PS, alpha = 0.05/length(seed_est),
                                          test_name = "chow1", intercept = TRUE)
    seedbs_prune_time[ll] <- proc.time()[[1]] - cur_time
    
    # stdBS
    cur_time <- proc.time()[[1]]
    std_est <- stdBS_ccp(X, Y, 1, n, min_gap = min_dis,
                         test_fun = test_wrap,
                         test_param = test_param,
                         esti_fun = esti_fun,
                         esti_param = esti_param,
                         verbose = TRUE)
    stdbs_all[[ll]] <- std_est
    stdbs_time[ll] <- proc.time()[[1]] - cur_time
    
    # stdBS with pruning
    stdbs_prune_all[[ll]] <- prune_candi(X, Y, std_est, PS, alpha = 0.05/length(std_est),
                                         test_name = "chow1", intercept = TRUE)
    stdbs_prune_time[ll] <- proc.time()[[1]] - cur_time
    
    # known number of CCPs
    cur_time <- proc.time()[[1]]
    bp1_all[[ll]] <- strucchange::breakpoints(Y ~ X, breaks = 2, h = rel_min_dis)$breakpoints + 1
    bp1_time[ll] <- proc.time()[[1]] - cur_time
    
    # unknown number of CPs
    cur_time <- proc.time()[[1]]
    bp2_all[[ll]] <- strucchange::breakpoints(Y ~ X, h = rel_min_dis)$breakpoints + 1
    bp2_time[ll] <- proc.time()[[1]] - cur_time
    
    # unknown number of CCPs with prune
    bp3_all[[ll]] <- prune_candi(X, Y, bp2_all[[ll]], PS, alpha = 0.05/length(bp2_all[[ll]]) ,
                                 test_name = "chow1", intercept = TRUE)
    bp3_time[ll] <- proc.time()[[1]] - cur_time
    
  }
  res_lst <- list(
    seedbs_all = seedbs_all,
    seedbs_prune_all = seedbs_prune_all,
    stdbs_all = stdbs_all,
    stdbs_prune_all = stdbs_prune_all,
    bp1_all = bp1_all,
    bp2_all = bp2_all,
    bp3_all = bp3_all,
    seedbs_time = seedbs_time,
    seedbs_prune_time = seedbs_prune_time,
    stdbs_time = stdbs_time,
    stdbs_prune_time = stdbs_prune_time,
    bp1_time = bp1_time,
    bp2_time = bp2_time,
    bp3_time = bp3_time
  )
  saveRDS(res_lst,
          file = paste0(pkg_dir, "/inst/output/ccpl_exp3",
                        "_allmethods",
                        "_2ccp_1nccp",
                        "_manualparams",
                        "_seed", seed_num,
                        "_rep", kk,
                        ".rds"))
  
}, future.seed = NULL, future.packages = c("dplyr"))
