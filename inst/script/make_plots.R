devtools::load_all()
library(dplyr)
library(ggplot2)

# ---- helper functions ----

#' Load experiment results based on one seed.
get_all_est_one_seed_tmp <- function(output_prefix,
                                     rel_min_dis,
                                     nccp_scale,
                                     seed_num, reps, n) {
  max_pvals_df_seed <- vector("list", length(reps))
  max_pvals2_df_seed <- vector("list", length(reps))
  score_df_seed <- vector("list", length(reps))
  score2_df_seed <- vector("list", length(reps))
  struc_df_seed <- vector("list", length(reps))
  for (kk in reps) {
    res_all_lst <- readRDS(paste0(output_prefix,
                                  "_relmindis", rel_min_dis,
                                  "_nccpscale", nccp_scale,
                                  "_n", n,
                                  "_seed", seed_num,
                                  "_rep", kk,
                                  ".rds"))
    max_pvals_df_seed[[kk]] <- res_all_lst$max_pvals_all %>%
      dplyr::mutate(est_ccp = max_pval < 0.05,
                    rep = kk,
                    seed_num = seed_num) %>%
      dplyr::slice_min(max_pval) %>%
      dplyr::mutate(est_value = dplyr::if_else(est_ccp,
                                               as.double(as.character(cp)),
                                               n + 1)) %>%
      dplyr::select(seed_num, cp_scale, n, rep, est_value) %>%
      dplyr::mutate(time = res_all_lst$candi1_time)

    max_pvals2_df_seed[[kk]] <- res_all_lst$max_pvals_all2 %>%
      dplyr::mutate(est_ccp = max_pval < 0.05,
                    rep = kk,
                    seed_num = seed_num) %>%
      dplyr::slice_min(max_pval) %>%
      dplyr::mutate(est_value = dplyr::if_else(est_ccp,
                                               as.double(as.character(cp)),
                                               n + 1)) %>%
      dplyr::select(seed_num, cp_scale, n, rep, est_value) %>%
      dplyr::mutate(time = res_all_lst$candi2_time)

    min_dis <- rel_min_dis * n
    grid_pts <- ceiling(seq(from = min_dis + 1, to = n - min_dis + 1, by = min_dis))
    score_df_seed[[kk]] <- data.frame(est_value = grid_pts[which.min(res_all_lst$score_all)],
                                      cp_scale = nccp_scale,
                                      n = n,
                                      rep = kk,
                                      seed_num = seed_num) %>%
      dplyr::select(seed_num, cp_scale, n, rep, est_value) %>%
      dplyr::mutate(time = res_all_lst$ccloss_time)

    struc_df_seed[[kk]] <- data.frame(est_value = dplyr::if_else(is.na(res_all_lst$struc_all),
                                                                 n + 1,
                                                                 res_all_lst$struc_all + 1),
                                      cp_scale = nccp_scale,
                                      n = n,
                                      rep = kk,
                                      seed_num = seed_num) %>%
      dplyr::select(seed_num, cp_scale, n, rep, est_value) %>%
      dplyr::mutate(time = res_all_lst$bp_time)
  }
  max_pvals_df_seed <- do.call(rbind, max_pvals_df_seed)
  max_pvals2_df_seed <- do.call(rbind, max_pvals2_df_seed)
  score_df_seed <- do.call(rbind, score_df_seed)
  struc_df_seed <- do.call(rbind, struc_df_seed)
  # add method names and combine them in one data.frame
  max_pvals_df_seed$method <- "OracleCandidates"
  max_pvals2_df_seed$method <- "BPCandidates"
  score_df_seed$method <- "CSL"
  struc_df_seed$method <- "Breakpoints"
  rbind(max_pvals_df_seed,
        max_pvals2_df_seed,
        score_df_seed,
        struc_df_seed)
}

get_est <- function(output_prefix, fixed_sign, cp_scale, seed_num, reps, ns,
                    methods) {
  est_all <- vector("list", length(reps))
  time_all <- vector("list", length(reps))
  for (kk in reps) {
    res_lst <- readRDS(paste0(output_prefix,
                              "_2ccp_1nccp",
                              "_manualparams",
                              "_seed", seed_num,
                              "_rep", kk,
                              ".rds"))
    res_df_lst <- lapply(1:7, \(rr) {
      res <- res_lst[[rr]]
      df <- do.call(rbind, lapply(seq_along(res), \(ll) {
        cbind(res[[ll]], rep(ns[ll], length(res[[ll]])))
      })) |> as.data.frame()
      colnames(df) <- c("est", "n")
      df %>%
        dplyr::mutate(method = methods[rr])
    })
    res_df <- do.call(rbind, res_df_lst)
    res_df[['rep']] <- kk
    est_all[[kk]] <- res_df
    
    time_df <- lapply(8:length(res_lst), \(rr) {
      res <- res_lst[[rr]]
      df <- data.frame(n = ns, time = res, method = methods[rr-7])
    }) %>%
      dplyr::bind_rows()
    time_all[[kk]] <- time_df
  }
  est_df <- do.call(rbind, est_all)
  time_df <- do.call(rbind, time_all)
  return(list(est_df = est_df,
              time_df = time_df))
}

# ---- Figure 6 ----
fixed_sign <- 0
nccp_scale <- 0.5
ns <- c(250, 1000, 4000)
(rel_true_rcp <- c(0.25, 0.5, 0.75))
(rel_true_ccp <- c(0.5))
rel_true_nccp <- setdiff(rel_true_rcp, rel_true_ccp)
ccp_scales <- 0.5
seed_num <- 1000
reps <- 1:200
n_seg <- 3
rel_min_seg <- 0.1
rel_min_dis <- 0.05
res_df <- lapply(ns, \(n) {
  min_dis <- n * rel_min_dis
  min_seg <- rel_min_seg * n
  grid_pts <- ceiling(seq(from = min_dis + 1, to = n - min_dis + 1, by = min_dis))
  output_prefix <- paste0("./inst/output/ccpl_exp1",
                          "_allmethods",
                          "_relminseg", rel_min_seg,
                          "_nseg", n_seg,
                          "_fixedsign", fixed_sign)
  get_all_est_one_seed_tmp(output_prefix,
                           # min_dis,
                           rel_min_dis,
                           nccp_scale,
                           seed_num, reps, n)
  }) %>%
    dplyr::bind_rows()

# plot the estimates
methods <- c("CSL", "OracleCandidates", "BPCandidates", "Breakpoints")
p <- res_df %>%
  dplyr::mutate(method = factor(method, levels = rev(methods)),
                rel_est_value = dplyr::if_else(est_value == n + 1,
                                               1, (est_value-1) / n),
                n = ordered(n, levels = rev(ns))) %>%
  ggplot(aes(x = rel_est_value, y = method,
             color = factor(method))) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, "no CCP \ndetected"),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0, 1)) +
  geom_vline(xintercept = rel_true_ccp, linetype = 2, color = 'red', linewidth = 0.5) +
  geom_vline(xintercept = rel_true_nccp, lty = 2, color = 'lightgrey') +
  ggbeeswarm::geom_quasirandom(dodge.width = 1, alpha = 0.5, size = 1) +
  facet_grid(n ~ ., switch = "y") +
  ylab("total number of time points") +
  xlab("relative location") +
  scale_color_brewer(name = "",
                     labels = c("CSL" = "LossCS",
                                "OracleCandidates" = "Prune-Oracle",
                                "BPCandidates" = "Prune-BP",
                                "Breakpoints" = "BP-1"),
                     type = "qual", palette = 6, direction = -1,
                     guide = guide_legend(reverse = TRUE))

extra_df <- res_df %>%
  dplyr::filter(n %in% ns) %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(no_est = sum(est_value == n + 1), .groups = "drop") %>%
  dplyr::filter(no_est > 0) %>%
  dplyr::mutate(rel_est_value = 1 - 0.05,
                n = ordered(n, levels = rev(ns))) # for location of the labels
p +
  geom_text(data = extra_df,
            size = 3.5,
            family = "serif",
            mapping = aes(x = rel_est_value, y = method, label = no_est),
            show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,-5,0),
        legend.text = element_text(size = 9, family = "mono"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 10, family = "serif"),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 10, family = "serif"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 5,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))

ggsave(glue::glue("./inst/output/plots/ccpl_exp1",
                  "_allmethods",
                  "_relminseg", rel_min_seg,
                  "_nseg", n_seg,
                  "_fixedsign", fixed_sign,
                  "_relmindis", rel_min_dis,
                  "_seed", seed_num,
                  ".pdf"),
       width = 6.16, height = 2.9)

res_df %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(median_time = median(time))

# ---- Figure 7 ----
seed_num <- ""
reps <- 1:200
cp_scale <- 0.5
rel_min_dis <- 0.05
fixed_sign <- 1
n <- 2000
rel_min_seg <- 0.1 # 0.05, 0.1
n_seg <- 3
rel_true_ccps <- seq(from = 0.1, to = 0.9, by = 0.1)
est_df_lst <- vector("list")
min_dis <- rel_min_dis * n
grid_pts <- seq(from = min_dis + 1, to = n - min_dis + 1, by = min_dis)
methods <- c("cc_loss", "cc_loss_noscale",
             "cc_loss_niklas1", "cc_loss_niklas2",
             "cc_loss_jonas1", "cc_loss_jonas2", "cc_loss_jonas3")
output_prefix <- glue::glue("./inst/output/ccpl_exp2_cclossniklas1_relminseg{rel_min_seg}_fixedsign{fixed_sign}_mindis{min_dis}_cpscale{cp_scale}_n{n}")

est_df <- lapply(reps, \(kk) {
  res <- readRDS(paste0(output_prefix,
                        "_seed", seed_num,
                        "_rep", kk, ".rds"))
  est_df <- data.frame(est_value = apply(res, 1, \(vals) {grid_pts[which.min(vals)]}),
                       rel_true_ccp = rel_true_ccps,
                       n = n,
                       rep = kk,
                       method = methods[3]) %>%
    dplyr::select(n, rep, est_value, rel_true_ccp, method)
}) %>%
  dplyr::bind_rows()

est_df <- est_df %>%
  dplyr::mutate(rel_est_ccp = (est_value - 1)/n,
                rel_rand_ccp = sample.int(n, dplyr::n(), replace = TRUE)/n) %>%
  dplyr::mutate(rel_est_ccp = dplyr::if_else(rel_est_ccp == 1,
                                             rel_rand_ccp,
                                             rel_est_ccp),
                method = factor(method, levels = methods)) %>%
  dplyr::select(-rel_rand_ccp) %>%
  identity()

est_df %>%
  ggplot(aes(x = rel_true_ccp,
             y = rel_est_ccp)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, color = "#984EA3", size = 1) +
  geom_abline(intercept = 0, slope = 1, color = 'grey') +
  geom_vline(xintercept = 0.25, lty = 2, color = 'lightgrey') +
  geom_vline(xintercept = 0.75, lty = 2, color = 'lightgrey') +
  geom_hline(yintercept = 0.25, lty = 2, color = 'lightgrey') +
  geom_hline(yintercept = 0.75, lty = 2, color = 'lightgrey') +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme(panel.grid.minor = element_blank()) +
  xlab("true relative location") +
  ylab("estimated relative location") +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,-5,0),
        strip.background = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 10, family = "serif"),
        axis.text = element_text(size = 9, family = "serif"),
        text = element_text(size = 10, family = "serif"),
        plot.margin = margin(t = 1,  # Top margin
                             r = 1,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))

ggsave(glue::glue("./inst/output/plots/ccpl_exp2",
                  "_relminseg", rel_min_seg,
                  "_fixedsign", fixed_sign,
                  "_relmindis", rel_min_dis,
                  "_cpscale", cp_scale,
                  ".pdf"),
       width = 3.2, height = 3)

# ---- Figure 8 ----
output_prefix <- "./inst/output/ccpl_exp3_allmethods"
seed_num <- 1000 # seed not actually used in the manual parameter experiments
fixed_sign <- 1
if (is.null(seed_num) & fixed_sign == 0) stop("Must give a seed number if fixed_sign is 0!")
ns <- c(1000, 2000, 4000)
rel_true_ccp <- c(0.2, 0.8)
rel_true_nccp <- c(0.5)
cp_scale <- 0.5
reps <- 1:200
methods <- c("LossCS-SeedBS", "LossCS-SeedBS\n-Prune",
             "LossCS-StdBS", "LossCS-StdBS\n-Prune",
             "BP-2", "BP", "Prune-BP")
res_lst <- get_est(output_prefix, fixed_sign, cp_scale, seed_num, reps, ns, methods)
est_df <- res_lst$est_df %>%
  dplyr::mutate(method = factor(method,
                                levels = c("BP-2", "BP", "Prune-BP",
                                           "LossCS-SeedBS", "LossCS-SeedBS\n-Prune",
                                           "LossCS-StdBS", "LossCS-StdBS\n-Prune")))

# plot the relative location of the estimates
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)
est_df %>%
  dplyr::mutate(rel_est_ccp = est / n,
                n = ordered(n, levels = rev(ns))) %>%
  ggplot(aes(x = rel_est_ccp)) +
  geom_histogram(position = "identity", binwidth = 0.05, alpha = 0.85) +
  geom_vline(xintercept = rel_true_ccp, color = 'red', lty = 2, linewidth = 0.3) +
  geom_vline(xintercept = rel_true_nccp, color = 'grey', lty = 2, linewidth = 0.3) +
  facet_grid(n ~ method) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                     guide = guide_axis(angle = 45)) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw() +
  xlab("relative location") +
  theme(panel.spacing = unit(0.1, "lines"),
        legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,-5,0),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text.x = element_text(size = 7, family = "mono"),
        strip.text.y = element_text(size = 7, family = "serif"),
        panel.grid.major.x = element_blank(), # remove vertical grid
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.1),
        panel.grid.minor.y = element_line(linewidth = 0.1),
        axis.text = element_text(size = 7, family = "serif"),
        text = element_text(size = 7, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))

ggsave(glue::glue("./inst/output/plots/ccpl_exp3_allmethods",
                  "_2ccp_1nccp",
                  "_fixedsign", fixed_sign,
                  "_cpscale", cp_scale,
                  "_seed", seed_num,
                  ".pdf"),
       width = 6.13, height = 3)

tbl_time <- res_lst$time_df %>%
  dplyr::group_by(method, n) %>%
  dplyr::summarise(median_time = median(time)) %>%
  tidyr::pivot_wider(id_cols = "n", names_from = "method", values_from = "median_time")

print(xtable::xtable(tbl_time), booktabs = TRUE, include.rownames = FALSE)

# ---- Figure C1 ----
output_prefix <- "./inst/output/ccpl_exp4"
fixed_sign <- 1
cp_scale <- 0.5
seed_num <- ""
reps <- 1:200
ns <- c(250, 1000, 4000) # c(500, 1000, 2000, 4000)
candset <- c("k1", "k2", "k1, k2") # c("1st", "2nd", "both")

res_all <- lapply(reps, \(kk) {
  lapply(ns, \(n) {
    res_lst <- readRDS(paste0(output_prefix,
                              "_fixedsign", fixed_sign,
                              "_cpscale", cp_scale,
                              "_n", n,
                              "_seed", seed_num,
                              "_rep", kk,
                              ".rds"))
    lapply(seq_along(res_lst), \(ll) {
      res <- res_lst[[ll]]
      res$cand <- candset[ll]
      if (candset[ll] == "k1") {
        res$est <- res$max_pval < 0.05
        res$rel_loc <- 0.25
      } else if (candset[ll] == "k2") {
        res$est <- res$max_pval < 0.05
        res$rel_loc <- 0.75
      }
      else {
        res$est <- res$max_pval < 0.05
        res$rel_loc <- c(0.25, 0.75)
      }
      res
    }) %>%
      dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows()
}) %>%
  dplyr::bind_rows()

binom_test <- function(x, n) {
  tst <- binom.test(x, n)
  tibble::tibble(est = tst$estimate,
                 lwr = tst$conf.int[1],
                 upr = tst$conf.int[2])
}

res_all %>%
  dplyr::mutate(cand = factor(cand, levels = candset)) %>%
  dplyr::group_by(n, cand, rel_loc) %>%
  dplyr::summarise(power = binom_test(sum(est), length(est))) %>%
  dplyr::mutate(n = factor(n)) %>%
  ggplot(aes(x = n, y = power$est, color = cand, shape = factor(rel_loc))) +
  geom_pointrange(aes(ymin = power$lwr, ymax = power$upr),
                  position = position_dodge(0.2)) +
  theme_bw() +
  labs(color = "candidate set", shape = "relative location",
       y = "power", x = "sample size") +
  theme(legend.position = "right",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,-5,0),
        legend.text = element_text(hjust = 0),
        text = element_text(size = 11, family = "serif"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 2)) +
  scale_color_discrete(labels = c("k1" = parse(text = "k[1]"),
                                  "k2" = parse(text = "k[2]"),
                                  "k1, k2" = parse(text = "k[1]*','~k[2]")))

ggsave(glue::glue("./inst/output/plots/ccpl_exp4",
                  "_fixedsign", fixed_sign,
                  "_cpscale", cp_scale,
                  "_seed", seed_num,
                  ".pdf"),
       width = 3.68, height = 2.5)

# ---- Figure C2 ----
output_prefix <- "./inst/output/ccpl_exp5"
fixed_sign <- 1
cp_scale <- 0.5
seed_num <- ""
n <- 4000
reps <- 1:200
rel_errs <- seq(from = -0.2, to = 0.2, by = 0.05)
rel_min_dis <- 0.05
min_dis <- rel_min_dis * n
n_seg <- 3
grid_pts <- seq(from = min_dis + 1, to = n - min_dis + 1, by = min_dis)

res_all <- lapply(reps, \(kk) {
  lapply(rel_errs, \(rel_err) {
    res_lst <- readRDS(paste0(output_prefix,
                              "_fixedsign", fixed_sign,
                              "_cpscale", cp_scale,
                              "_n", n,
                              "_relerr", rel_err,
                              "_seed", seed_num,
                              "_rep", kk,
                              ".rds"))
    est <- res_lst$max_pvals$max_pval < 0.05
    data.frame(est = est,
               n = n, rep = kk,
               rel_err = rel_err,
               max_pval_left = max(res_lst$pvals_left),
               max_pval_right = max(res_lst$pvals_right))
  }) %>%
    dplyr::bind_rows()
}) %>%
  dplyr::bind_rows()
    
# plot percentage of left or right reject there is an invariant set
res_sum1 <- res_all %>%
  dplyr::mutate(rel_err = ordered(rel_err),
                rej_left = max_pval_left < 0.05,
                rej_right = max_pval_right < 0.05) %>%
  dplyr::group_by(n, rel_err) %>%
  dplyr::summarise(pct_ccp = sum(est)/length(reps),
                   pct_rej_left = sum(rej_left)/length(reps),
                   pct_rej_right = sum(rej_right)/length(reps))

res_sum2 <- res_all %>%
  dplyr::mutate(rel_err = ordered(rel_err),
                rej_left = max_pval_left < 0.05,
                rej_right = max_pval_right < 0.05,
                rej = rej_left | rej_right) %>%
  dplyr::filter(est == TRUE) %>%
  dplyr::group_by(n, rel_err) %>%
  dplyr::summarise(pct_rej = sum(rej)/n(),
                   pct_acc = 1 - pct_rej) %>%
  tidyr::pivot_longer(cols = starts_with("pct"),
                      names_to = "category",
                      values_to = "percentage")

res_sum2 %>%
  ggplot() +
  geom_bar(aes(x = rel_err, y = percentage, fill = category),
           data = res_sum2,
           stat = "identity",
           inherit.aes = FALSE,
           alpha = 0.8,
           width = 0.5) +
  geom_point(aes(x = rel_err, y = pct_ccp, shape = "reject"),
             data = res_sum1,
             size = 3) +
  theme_bw() +
  labs(x = "relative error",
       y = "percentage",
       shape = "",
       fill = "exist inv. sets (left or right)") +
  scale_fill_hue(labels = c("not rejected" , "rejected")) +
  scale_shape_manual(labels = '\"not a CCP\" rejected', values = "x") +
  theme(text = element_text(size = 11, family = "serif"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 2)) +
  guides(fill = guide_legend(reverse=TRUE))
    

ggsave(glue::glue("./inst/output/plots/ccpl_exp5",
                  "_shiftedcandidate",
                  "_fixedsign", fixed_sign,
                  "_cpscale", cp_scale,
                  "_seed", seed_num,
                  ".pdf"),
       width = 4.9, height = 2.5)

