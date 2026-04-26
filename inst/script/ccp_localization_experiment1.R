## Figure 6 -- with DPDU (Xu, Wang, Zhao, Yu 2022, arXiv:2207.12453) competitor ##

# ---- setup ----
if (!require("renv")) {
  install.packages("renv", repos = "http://cran.us.r-project.org")
}
renv::restore()
if (!require("dHSIC")) {
  install.packages("dHSIC", repos = "http://cran.us.r-project.org")
}
if (!require("changepoints")) {
  install.packages("changepoints", repos = "http://cran.us.r-project.org")
}

# ---- libraries ----
pkg_dir <- getwd()
devtools::install_local(pkg_dir) # Note: only need to install once
devtools::load_all(pkg_dir)

# set up parallel
library(future.apply)
plan(multisession)

# ---- helper functions ----
# DPDU pipeline: CV.search -> DPDU2 -> local.refine (Xu et al. 2022)
dpdu_fit <- function(X, Y,
                     lambda_set = c(0.3, 0.5, 1, 2),
                     zeta_set = c(10, 15, 20),
                     w = 0.9) {
  temp <- changepoints::CV.search.DPDU.regression(
    y = Y, X = X, lambda_set = lambda_set, zeta_set = zeta_set
  )
  min_idx <- as.vector(arrayInd(which.min(temp$test_error), dim(temp$test_error)))
  cpt_init <- unlist(temp$cpt_hat[min_idx[1], min_idx[2]])
  if (length(cpt_init) == 0) return(integer(0))
  beta_hat <- matrix(unlist(temp$beta_hat[min_idx[1], min_idx[2]]),
                     ncol = length(cpt_init) + 1)
  changepoints::local.refine.DPDU.regression(cpt_init, beta_hat, Y, X, w = w)
}

# ---- main script ----
arg <- commandArgs(trailingOnly = TRUE)
n <- as.integer(arg[1])
message("n: ", n)
seed_num <- as.integer(arg[2])
message("seed_num: ", seed_num)
fixed_sign <- 0
if (is.null(seed_num) & fixed_sign == 0) stop("Must give a seed number if fixed_sign is 0!")
reps <- 1:200
PS <- get_powerset(as.character(1:4))
set.seed(seed_num)
(rel_true_rcp <- c(0.25, 0.5, 0.75))
(rel_true_ccp <- c(0.5)) # Note: more than one does not work directly with the score approach here
true_ccp <- ceiling(n * rel_true_ccp + 1)
true_rcp <- ceiling(n * rel_true_rcp + 1)
true_nccp <- setdiff(true_rcp, true_ccp)
cand_set <- unique(sort(c(true_rcp, true_ccp)))
segs <- get_segs(cand_set, n)
rel_min_dis <- 0.05
min_dis <- rel_min_dis * n # grid distance for evaluating the losses
grid_pts <- ceiling(seq(from = min_dis + 1, to = n - min_dis + 1, by = min_dis))
nccp_scale <- 0.5
ccp_scales <- 0.5 # seq(from = 0.3, to = 0.7, by = 0.2)
rel_min_seg <- 0.1
min_seg <- rel_min_seg * n # 20 # min_seg for the three new losses
n_seg <- 3 # the number of segments m for cc_loss

future_lapply(reps, \(kk) {
  devtools::load_all(pkg_dir)
  max_pvals_lst <- vector("list", length(ccp_scales))
  max_pvals_lst2 <- vector("list", length(ccp_scales))
  score_lst <- vector("list", length(ccp_scales))
  struc_est <- rep(NA, length(ccp_scales))
  dpdu_lst <- vector("list", length(ccp_scales))
  candi1_time <- vector("numeric", length(ccp_scales))
  candi2_time <- vector("numeric", length(ccp_scales))
  ccloss_time <- vector("numeric", length(ccp_scales))
  bp_time <- vector("numeric", length(ccp_scales))
  dpdu_time <- vector("numeric", length(ccp_scales))

  for (ll in seq_along(ccp_scales)) {
    params_lst <- rand_mb_params_step(n, true_ccp, true_nccp, seed_num,
      nccp_mean_nc = 4, nccp_mean_sc = nccp_scale,
      nccp_sig_nc = 4, nccp_sig_sc = nccp_scale,
      nccp_coef_nc = 3, nccp_coef_sc = nccp_scale,
      ccp_mean_nc = 1, ccp_mean_sc = ccp_scales[ll],
      ccp_sig_nc = 1, ccp_sig_sc = ccp_scales[ll],
      ccp_coef_nc = 2, ccp_coef_sc = ccp_scales[ll],
      fixed_sign = fixed_sign
    )

    dat <- with(
      params_lst$params_df,
      sim_data_mix(
        m1_all = m1, m2_all = m2, m3_all = m3, m4_all = m4, m5_all = m5,
        s1_all = s1, s2_all = s2, s3_all = s3, s4_all = s4, s5_all = s5,
        a12_all = a12, a53_all = a53, a43_all = a43,
        b15_all = b15, b25_all = b25,
        seed_num = kk
      )
    )
    colnames(dat) <- as.character(1:5)

    # candidate approach (known RCPs)
    pvals <- matrix(NA, nrow = length(cand_set), ncol = length(PS))
    rownames(pvals) <- as.character(cand_set)
    colnames(pvals) <- sapply(PS, \(l) {
      paste(l, collapse = "_")
    })
    cur_time <- proc.time()[[1]]
    for (ii in seq_along(cand_set)) {
      for (jj in seq_along(PS)) {
        S <- PS[[jj]]
        pvals[ii, jj] <- test_chow1(dat[, 1:4], dat[, 5], segs[[ii]], segs[[ii + 1]], S)$p_value
        # pvals[ii, jj] <- test_chow2(dat[,1:4], dat[,5], segs[[ii]], segs[[ii+1]], S)$p_value
      }
    }
    candi1_time[ll] <- proc.time()[[1]] - cur_time
    max_pvals <- data.frame(
      cp = rownames(pvals),
      max_pval = apply(pvals, 1, max),
      max_pval_set = sapply(PS[apply(pvals, 1, which.max)], \(l) {
        paste(l, collapse = "_")
      })
    )
    max_pvals <- max_pvals %>%
      dplyr::mutate(
        is_ccp = cp %in% true_ccp,
        cp_scale = ccp_scales[ll],
        n = n,
        rep = kk
      )
    rownames(max_pvals) <- NULL
    max_pvals_lst[[ll]] <- max_pvals

    # candidate approach (estimated number of RCPs)
    cur_time <- proc.time()[[1]]
    cand_set2 <- strucchange::breakpoints(dat[, 5] ~ dat[, 1:4], h = 100, breaks = 5)$breakpoints
    if (!anyNA(cand_set2)) {
      segs2 <- get_segs(cand_set2, n)
      pvals2 <- matrix(NA, nrow = length(cand_set2), ncol = length(PS))
      rownames(pvals2) <- as.character(cand_set2)
      colnames(pvals2) <- sapply(PS, \(l) {
        paste(l, collapse = "_")
      })
      for (ii in seq_along(cand_set2)) {
        for (jj in seq_along(PS)) {
          S <- PS[[jj]]
          pvals2[ii, jj] <- test_chow1(dat[, 1:4], dat[, 5], segs2[[ii]], segs2[[ii + 1]], S)$p_value
          # pvals2[ii, jj] <- test_chow2(dat[,1:4], dat[,5], segs2[[ii]], segs2[[ii+1]], S)$p_value
        }
      }
      candi2_time[ll] <- proc.time()[[1]] - cur_time
      max_pvals2 <- data.frame(
        cp = rownames(pvals2),
        max_pval = apply(pvals2, 1, max),
        max_pval_set = sapply(PS[apply(pvals2, 1, which.max)], \(l) {
          paste(l, collapse = "_")
        })
      )
      max_pvals2 <- max_pvals2 %>%
        dplyr::mutate(
          is_ccp = cp %in% true_ccp,
          cp_scale = ccp_scales[ll],
          n = n,
          rep = kk
        )
      rownames(max_pvals2) <- NULL
    } else {
      max_pvals2 <- data.frame(
        cp = n + 1,
        max_pval = 1,
        max_pval_set = "",
        is_ccp = FALSE,
        cp_scale = ccp_scales[ll],
        n = n,
        rep = kk
      )
      candi2_time[ll] <- proc.time()[[1]] - cur_time
    }
    max_pvals_lst2[[ll]] <- max_pvals2

    # niklas1 score
    cur_time <- proc.time()[[1]]
    score_lst[[ll]] <- sapply(grid_pts, \(t) cc_loss(dat[, 1:4], dat[, 5], t, 1, n,
      min_seg = min_seg,
      cond_sets = PS,
      intercept = TRUE,
      verbose = FALSE
    ))
    ccloss_time[ll] <- proc.time()[[1]] - cur_time


    # structural change
    cur_time <- proc.time()[[1]]
    struc_est[ll] <- strucchange::breakpoints(dat[, 5] ~ dat[, 1:4], h = min_seg, breaks = 1)$breakpoints
    bp_time[ll] <- proc.time()[[1]] - cur_time

    # DPDU regression (Xu, Wang, Zhao, Yu 2022, arXiv:2207.12453)
    cur_time <- proc.time()[[1]]
    dpdu_lst[[ll]] <- dpdu_fit(dat[, 1:4], dat[, 5]) + 1
    dpdu_time[ll] <- proc.time()[[1]] - cur_time
  }

  res_lst <- list(
    max_pvals_all = do.call(rbind, max_pvals_lst),
    max_pvals_all2 = do.call(rbind, max_pvals_lst2),
    score_all = do.call(rbind, score_lst),
    struc_all = struc_est,
    dpdu_all = dpdu_lst,
    candi1_time = candi1_time,
    candi2_time = candi2_time,
    ccloss_time = ccloss_time,
    bp_time = bp_time,
    dpdu_time = dpdu_time
  )

  saveRDS(res_lst,
    file = paste0(
      pkg_dir, "/inst/output/ccpl_exp1-1",
      "_allmethods",
      "_relminseg", rel_min_seg,
      "_nseg", n_seg,
      "_fixedsign", fixed_sign,
      "_relmindis", rel_min_dis,
      "_nccpscale", nccp_scale,
      "_n", n,
      "_seed", seed_num,
      "_rep", kk,
      ".rds"
    )
  )
}, future.seed = TRUE, future.packages = c("dplyr"))
