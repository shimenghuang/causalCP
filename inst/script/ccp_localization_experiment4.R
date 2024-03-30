## Figure C1 ##

# Experiments with candidate approach: 2 CCPs and only one is a candidate
pkg_dir <- "/Users/hrt620/Documents/ku_projects/causalCPD_proj/causalCPD"
devtools::load_all(pkg_dir)

# set up parallel
library(future.apply)
plan(multisession)

reps <- 1:200
seed_num <- NULL
fixed_sign <- 1
if (is.null(seed_num) & fixed_sign == 0) stop("Must give a seed number if fixed_sign is 0!")
ns <- c(250, 1000, 4000)
cp_scale <- 0.5
PS <- get_powerset(as.character(1:4))
rel_true_nccp <- c()
rel_true_ccp <- c(0.25, 0.75)

future_lapply(reps, \(kk) {
  devtools::load_all(pkg_dir)
  for (n in ns) {
    max_pvals_lst <- vector("list", 3)
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
    #                                                                              s1, s2, s3, s4, s5, S))
    #   message(beta2)
    #   beta_diff[ss] <- any(beta1 != beta2)
    # }
    # beta_diff <- rep(NA, length(ps_tmp))
    # for (ss in seq_along(ps_tmp)) {
    #   S <- ps_tmp[[ss]]
    #   message("**", S, "**")
    #   beta1 <- with(params_lst$params_all[as.character(true_ccp[1]),], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
    #                                                           s1, s2, s3, s4, s5, S))
    #   message(beta1)
    #   beta2 <- with(params_lst$params_all[as.character(true_ccp[2]),], comp_pop_ols(a12, b15, b25, a53, a43, m1, m2, m3, m4, m5,
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

    # candidate approach with one true CCP as the candidate
    for (ll in 1:2) {
      cand_set <- true_ccp[ll]
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
      max_pvals_lst[[ll]] <- max_pvals
    }

    # candidate approach with both true CCP as the candidates
    cand_set <- true_ccp
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
    max_pvals_lst[[3]] <- max_pvals

    saveRDS(max_pvals_lst,
            file = paste0(pkg_dir, "/inst/output/ccpl_exp4",
                          "_fixedsign", fixed_sign,
                          "_cpscale", cp_scale,
                          "_n", n,
                          "_seed", seed_num,
                          "_rep", kk,
                          ".rds"))
  }
}, future.seed = NULL)

