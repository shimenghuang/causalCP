## Figure 7 ##

pkg_dir <- "/Users/hrt620/Documents/ku_projects/causalCPD_proj/causalCPD"
devtools::load_all(pkg_dir)

# set up parallel
library(future.apply)
plan(multisession)

# ---- main script ----
arg <- commandArgs(trailingOnly = TRUE)
reps <- 1:200
seed_num <- NULL
fixed_sign <- 1
if (is.null(seed_num) & fixed_sign == 0) stop("Must give a seed number if fixed_sign is 0!")
ns <- c(2000)
cp_scale <- 0.5
rel_min_dis <- 0.05
PS <- get_powerset(as.character(1:4))
rel_true_nccp <- c(0.25, 0.75)
rel_true_ccps <- seq(from = 0.1, to = 0.9, by = 0.1)
rel_min_seg <- 0.1

future_lapply(reps, \(kk) {
  devtools::load_all(pkg_dir)
  for (n in ns) {
    min_dis <- rel_min_dis * n
    min_seg <- rel_min_seg * n
    score_lst <- vector("list", length(rel_true_ccps))
    for (ll in seq_along(rel_true_ccps)) {
      true_ccp <- n * rel_true_ccps[ll] + 1
      true_nccp <- n * rel_true_nccp + 1
      true_rcp <- sort(c(true_ccp, true_nccp))
      grid_pts <- seq(from = min_dis + 1, to = n - min_dis + 1, by = min_dis)
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
      #   beta2 <- with(params_lst$params_all[as.character(true_ccp),], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
      #                                                                              s1, s2, s3, s4, s5, S))
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
      
      # niklas1 score
      score_lst[[ll]] <- sapply(grid_pts, \(t) cc_loss(dat[,1:4], dat[,5], t, 1, n,
                                                       min_seg = min_seg,
                                                       cond_sets = PS,
                                                       intercept = TRUE,
                                                       verbose = FALSE))
    }
    saveRDS(do.call(rbind, score_lst),
            file = paste0(pkg_dir, "/inst/output/ccpl_exp2",
                          "_cclossniklas1",
                          "_relminseg", rel_min_seg,
                          "_fixedsign", fixed_sign,
                          "_mindis", min_dis,
                          "_cpscale", cp_scale,
                          "_n", n,
                          "_seed", seed_num,
                          "_rep", kk,
                          ".rds"))
  }
}, future.seed = TRUE, future.packages = c("dplyr"))


