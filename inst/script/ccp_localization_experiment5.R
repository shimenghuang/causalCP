## Figure C2 ##

# Experiments with candidate approach: 1 CCPs candidate that is inaccurate
pkg_dir <- "/Users/hrt620/Documents/ku_projects/causalCPD_proj/causalCPD"
devtools::load_all(pkg_dir)

# set up parallel
library(future.apply)
plan(multisession)

reps <- 1:200
seed_num <- NULL
fixed_sign <- 1
if (is.null(seed_num) & fixed_sign == 0) stop("Must give a seed number if fixed_sign is 0!")
rel_errs <- seq(from = -0.2, to = 0.2, by = 0.05)
cp_scale <- 0.5
n <- 4000
PS <- get_powerset(as.character(1:4))
rel_true_nccp <- c()
rel_true_ccp <- 0.5
rel_min_dis <- 0.05
min_dis <- rel_min_dis * n
n_seg <- 3

split_half <- function(x) {
  len1 <- floor(length(x)/2)
  list(x[1:len1],
       x[(len1+1):length(x)])
}

future_lapply(reps, \(kk) {
  devtools::load_all(pkg_dir)
  true_ccp <- n * rel_true_ccp + 1
  true_nccp <- n * rel_true_nccp + 1
  true_rcp <- sort(c(true_ccp, true_nccp))
  params_lst <- rand_mb_params_step(n, true_ccp, true_nccp, seed_num = seed_num,
                                    nccp_mean_nc = 4, nccp_mean_sc = cp_scale,
                                    nccp_sig_nc = 4, nccp_sig_sc = cp_scale,
                                    nccp_coef_nc = 3, nccp_coef_sc = cp_scale,
                                    ccp_mean_nc = 1, ccp_mean_sc = cp_scale,
                                    ccp_sig_nc = 1, ccp_sig_sc = cp_scale,
                                    ccp_coef_nc = 2, ccp_coef_sc = cp_scale,
                                    fixed_sign = fixed_sign)
  
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
  #   beta1 <- with(params_lst$params_all['1',], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
  #                                                           s1, s2, s3, s4, s5, S))
  #   message(beta1)
  #   beta2 <- with(params_lst$params_all[as.character(true_ccp[1]),], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
  #                                                                                 s1, s2, s3, s4, s5, S))
  #   message(beta2)
  #   beta_diff[ss] <- any(beta1 != beta2)
  # }

  dat <- with(params_lst$params_df,
              sim_data_mix(m1_all = m1, m2_all = m2, m3_all = m3, m4_all = m4, m5_all = m5,
                           s1_all = s1, s2_all = s2, s3_all = s3, s4_all = s4, s5_all = s5,
                           a12_all = a12, a53_all = a53, a43_all = a43,
                           b15_all = b15, b25_all = b25,
                           seed_num = kk))
  colnames(dat) <- as.character(1:5)

  for (rel_err in rel_errs) {
    # candidate approach with one true CCP as the candidate
    cand_set <- true_ccp + n * rel_err
    segs <- get_segs(cand_set, n)
    pvals <- matrix(NA, nrow = length(cand_set), ncol = length(PS))
    rownames(pvals) <- as.character(cand_set)
    colnames(pvals) <- sapply(PS, \(l) {paste(l, collapse = "_")})
    for (ii in seq_along(cand_set)) {
      for (jj in seq_along(PS)) {
        S <- PS[[jj]]
        pvals[ii, jj] <- test_chow1(dat[,1:4], dat[,5], segs[[ii]], segs[[ii+1]], S)$p_value
      }
    }
    max_pvals <- data.frame(cp = rownames(pvals),
                            max_pval = apply(pvals, 1, max),
                            max_pval_set = sapply(PS[apply(pvals, 1, which.max)], \(l) {paste(l, collapse = "_")}))
    max_pvals <- max_pvals %>%
      dplyr::mutate(is_ccp = cp %in% true_ccp,
                    cp_scale = cp_scale,
                    n = n,
                    rep = kk)
    rownames(max_pvals) <- NULL

    # test the left and right intervals of the candidate
    idx_left <- 1:(cand_set[1]-1)
    idx_right <- cand_set[1]:n
    segs_left <- split_half(idx_left)
    segs_right <- split_half(idx_right)
    pvals_left <- two_env_test(dat[,1:4], dat[,5], segs_left[[1]], segs_left[[2]], ps = PS)$pvals
    pvals_right <- two_env_test(dat[,1:4], dat[,5], segs_right[[1]], segs_right[[2]], ps = PS)$pvals

    saveRDS(list(max_pvals = max_pvals,
                 pvals_left = pvals_left,
                 pvals_right = pvals_right),
            file = paste0(pkg_dir, "/inst/output/ccpl_exp5",
                          "_fixedsign", fixed_sign,
                          "_cpscale", cp_scale,
                          "_n", n,
                          "_relerr", rel_err,
                          "_seed", seed_num,
                          "_rep", kk,
                          ".rds"))
  }
}, future.seed = NULL)

