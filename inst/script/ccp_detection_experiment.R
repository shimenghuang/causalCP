## Figure 5 ##

library(dplyr)
devtools::load_all()

# ---- helper functions ----

test_wrap <- function(X, Y, I, PS, rel_grid, test_name = "chow1",
                      intercept = TRUE, rev_one_rest = FALSE, two_dir = FALSE,
                      save_resids = TRUE, alpha = 0.05) {
  res <- one_vs_rest_test(X[I,], Y[I], PS, rel_grid, test_name,
                          intercept, rev_one_rest, two_dir, save_resids)
  min(sapply(res$pvals_all, max))
}

binom_test <- function(x, n) {
  tst <- binom.test(x, n)
  tibble::tibble(est = tst$estimate,
                 lwr = tst$conf.int[1],
                 upr = tst$conf.int[2])
}

# ---- level + power of one split when the CCP is at different locations ----
library(future.apply)
plan(multisession)

pkg_dir <- "/Users/hrt620/Documents/ku_projects/causalCP"
test_name <- "chow1"
fixed_sign <- 1
seed_num <- NULL
n_rep <- 200
rel_true_ccps <- c(NA, seq(from = 0.1, to = 0.9, length.out = 9))
rel_true_nccp <- c()
ns <- c(250, 1000, 4000)
future_lapply(rel_true_ccps, \(rel_true_ccp) {
  devtools::load_all(pkg_dir)
  for (n in ns) {
    message("n = ", n)
    if (is.na(rel_true_ccp)) {
      true_ccp <- numeric()
    } else {
      true_ccp <- n * rel_true_ccp + 1
    }
    (true_nccp <- n * rel_true_nccp + 1)
    (true_rcp <- sort(c(true_ccp, true_nccp)))
    pvals1 <- rep(NA, n_rep)
    pvals2 <- rep(NA, n_rep)
    cp_scale <- ifelse(is.na(rel_true_ccp), 0, 1)
    for (kk in 1:n_rep) {
      params_lst <- rand_mb_params_step(n, true_ccp, true_nccp, seed_num = seed_num,
                                        nccp_mean_nc = 0, nccp_mean_sc = NA,
                                        nccp_sig_nc = 0, nccp_sig_sc = NA,
                                        nccp_coef_nc = 0, nccp_coef_sc = NA,
                                        ccp_mean_nc = 0, ccp_mean_sc = NA,
                                        ccp_sig_nc = 0, ccp_sig_sc = NA,
                                        ccp_coef_nc = 2, ccp_coef_sc = cp_scale,
                                        fixed_sign = fixed_sign)
      dat <- with(params_lst$params_df,
                  sim_data_mix(m1_all = m1, m2_all = m2, m3_all = m3, m4_all = m4, m5_all = m5,
                               s1_all = s1, s2_all = s2, s3_all = s3, s4_all = s4, s5_all = s5,
                               a12_all = a12, a53_all = a53, a43_all = a43,
                               b15_all = b15, b25_all = b25,
                               seed_num = kk))
      colnames(dat) <- as.character(1:5)
      X <- dat[,1:4]
      Y <- dat[,5]
      PS <- get_powerset(as.character(1:4))
      pvals1[kk] <- test_wrap(X, Y, 1:n, PS, rel_grid = c(1/2), test_name = test_name)

      params_lst <- rand_mb_params_step(n, true_ccp, true_nccp, seed_num = seed_num,
                                        nccp_mean_nc = 0, nccp_mean_sc = NA,
                                        nccp_sig_nc = 0, nccp_sig_sc = NA,
                                        nccp_coef_nc = 0, nccp_coef_sc = NA,
                                        ccp_mean_nc = 0, ccp_mean_sc = NA,
                                        ccp_sig_nc = 1, ccp_sig_sc = cp_scale,
                                        ccp_coef_nc = 0, ccp_coef_sc = NA,
                                        fixed_sign = fixed_sign)
      dat <- with(params_lst$params_df,
                  sim_data_mix(m1_all = m1, m2_all = m2, m3_all = m3, m4_all = m4, m5_all = m5,
                               s1_all = s1, s2_all = s2, s3_all = s3, s4_all = s4, s5_all = s5,
                               a12_all = a12, a53_all = a53, a43_all = a43,
                               b15_all = b15, b25_all = b25,
                               seed_num = kk))
      colnames(dat) <- as.character(1:5)
      X <- dat[,1:4]
      Y <- dat[,5]
      PS <- get_powerset(as.character(1:4))
      pvals2[kk] <- test_wrap(X, Y, 1:n, PS, rel_grid = c(1/2), test_name = test_name)
    }
    saveRDS(list(pvals1 = pvals1,
                 pvals2 = pvals2),
            file = paste0(pkg_dir, "/inst/output/ccpd_exp",
                          "_testname", test_name,
                          "_reltrueccp", rel_true_ccp,
                          "_nonccp",
                          "_n", n,
                          "_fixedsign", fixed_sign,
                          "_cpscale", cp_scale,
                          "_seed", seed_num,
                          ".rds"))
  }
}, future.seed = NULL, future.packages = c("dplyr"))

# ---- check results ----
library(dplyr)
library(ggplot2)

seed_num <- ""
fixed_sign <- 1
test_name <- "chow1"
n_rep <- 200
rel_true_ccps <- c(NA, seq(from = 0.1, to = 0.9, length.out = 9))
rel_true_nccp <- c()
ns <- c(250, 1000, 4000)
output_prefix <- "./inst/output/ccpd_exp"
res_all <- lapply(ns, \(n) {
  lapply(seq_along(rel_true_ccps), \(ll) {
    cp_scale <- ifelse(is.na(rel_true_ccps[ll]), 0, 1)
    pvals <- readRDS(paste0(output_prefix,
                            "_testname", test_name,
                            "_reltrueccp", rel_true_ccps[ll],
                            "_nonccp",
                            "_n", n,
                            "_fixedsign", fixed_sign,
                            "_cpscale", cp_scale,
                            "_seed", seed_num,
                            ".rds"))
    df <- data.frame(pval = pvals$pvals1)
    df$rep <- 1:n_rep
    df$rel_true_ccp <- rel_true_ccps[ll]
    df$n <- n
    df
  }) %>%
    bind_rows()
})  %>%
  bind_rows()

# plot all when num_split with error bar (version 3)
res_all_sum <- res_all %>%
  dplyr::group_by(n, rel_true_ccp) %>%
  dplyr::mutate(rej = pval < 0.05,
                rel_true_ccp = dplyr::if_else(is.na(rel_true_ccp), "no CCP", as.character(rel_true_ccp)),
                rel_true_ccp = factor(rel_true_ccp, levels = c("no CCP", as.character(seq(from = 0.1, to = 0.9, by = 0.1)))),
                n = factor(n, levels = rev(ns))) %>%
  dplyr::summarize(power = binom_test(sum(rej), length(rej)), .groups = "drop") %>%
  dplyr::mutate(alpha = 0.05)

res_all_sum %>%
  ggplot2::ggplot(aes(x = rel_true_ccp, y = power$est, shape = n, group = n)) +
  ggplot2::geom_point(position = position_dodge(0.55)) +
  ggplot2::geom_line(data = res_all_sum %>%
                       dplyr::filter(rel_true_ccp != "no CCP"),
                     aes(x = rel_true_ccp, y = power$est),
                     alpha = 0.2, linewidth = 0.2,
                     position = position_dodge(0.55)) +
  facet_grid(.~ rel_true_ccp != "no CCP", scales = "free_x", space = "free") +
  ggplot2::geom_errorbar(aes(ymin = power$lwr, ymax = power$upr),
                         position = position_dodge(0.55),
                         width = 0.2,
                         alpha = 0.5) +
  ggplot2::geom_hline(data = res_all_sum %>%
                        dplyr::filter(rel_true_ccp == "no CCP"),
                      aes(yintercept = alpha), color = 'red', lty = 2, 
                      linewidth = 0.2) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::ylab("fraction of rejection") +
  ggplot2::xlab("relative location of the true CCP") +
  ggplot2::labs(shape = "n") +
  theme_bw() +
  theme(text = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 10, family = "serif"),
        axis.text = element_text(size = 10, family = "serif"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave(paste0("inst/output/plots/ccpd_exp",
              "_testname", test_name,
              "_nonccp",
              "_fixedsign", fixed_sign,
              "_cpscale", 1,
              "_seed", seed_num,
              ".pdf"),
       width = 4.9, height = 2.5)

