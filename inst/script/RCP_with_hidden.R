library(strucchange)
library(dplyr)
library(ggplot2)

# ---- helper function ----

make_plot <- function(
    res_df, n_tot, levels, n_sim = 50, bin_w = 50,
    red_vdash = NULL, gray_vdash = NULL) {
  first_bin_center <- bin_w / 2
  na_center <- first_bin_center - bin_w
  df_na <- res_df %>%
    filter(is.na(est)) %>%
    mutate(est = na_center)
  df_na$d <- factor(df_na$d, levels = levels, labels = levels, ordered = TRUE)

  df_num <- res_df %>% filter(!is.na(est))
  df_num$d <- factor(df_num$d, levels = levels, labels = levels, ordered = TRUE)

  # Custom labeller function
  custom_labeller <- function(value) paste0("beta[2] - beta[1] == ", value)

  # Convert it to a labeller for ggplot
  my_labeller <- as_labeller(custom_labeller, default = label_parsed)

  ggplot() +
    # NA bar aligned with binwidth
    geom_bar(
      data = df_na,
      aes(x = est),
      width = bin_w,
      position = "dodge",
      color = "gray",
      alpha = 0.75
    ) +
    # Histogram for numeric
    geom_histogram(
      data = df_num,
      aes(x = est),
      binwidth = bin_w,
      boundary = 0, # start bins at 0
      position = "dodge",
      color = "gray",
      alpha = 0.75
    ) +
    geom_vline(xintercept = red_vdash, linetype = 2, color = "red") +
    geom_vline(xintercept = gray_vdash, linetype = 2, color = "gray") +
    coord_cartesian(xlim = c(na_center - bin_w / 2, n_tot)) +
    scale_x_continuous(
      breaks = c(na_center, seq(0, n_tot, 100)),
      labels = c("No RCP", "", seq(0, n_tot, 100)[-1])
    ) +
    scale_y_continuous(
      breaks = seq(0, n_sim, 25),
      labels = seq(0, n_sim, 25),
      limits = c(0, 50)
    ) +
    facet_wrap(d ~ ., nrow = 4, labeller = my_labeller) +
    labs(x = "Estimate") +
    # scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.margin = margin(b = -5),
      panel.grid.minor.x = element_blank(),
      # strip.text = element_blank(),
      # strip.background = element_blank(),
      # strip.placement = "outside",
      strip.switch.pad.grid = unit(0, "pt"),
      strip.switch.pad.wrap = unit(0, "pt")
    )
}

# ---- Localizing one RCP given a stationary hidden confounder ----

set.seed(42)

num_sim <- 50
n1 <- 400
n2 <- 500
dd <- c(0, 0.5, 1, 1.5)
res <- lapply(1:num_sim, \(ii) {
  lst <- lapply(dd, \(d) {
    H <- rnorm(n1 + n2)
    X <- H + rnorm(n1 + n2)
    b1 <- 0.5
    b2 <- b1 + d
    Y1 <- H[1:n1] + X[1:n1] * b1 + rnorm(n1)
    Y2 <- H[(n1 + 1):(n1 + n2)] + X[(n1 + 1):(n1 + n2)] * b2 + rnorm(n2)
    Y <- c(Y1, Y2)
    est <- breakpoints(Y ~ X, h = 0.2, breaks = 1)$breakpoints
    if (anyNA(est)) {
      data.frame(d = d, est = NA, est_num = NA)
    } else {
      data.frame(d = d, est = sort(est), est_num = seq_len(length(est)))
    }
  })
  tmp <- do.call(rbind, lst)
  tmp$sim <- ii
  tmp
})
res_df <- do.call(rbind, res)
write.csv(res_df,
  "../output/rcp_with_stationary_hidden_res.csv",
  row.names = FALSE
)

make_plot(res_df,
  n_tot = n1 + n2,
  levels = as.character(dd), bin_w = 20, red_vdash = n1
)
ggsave("../output/plots/rcp_with_stationary_hidden_hist.pdf",
  width = 4, height = 4
)

# ---- Localizing one RCP given a non-stationary hidden ----

set.seed(42)

num_sim <- 50
n1 <- 300 # before RCP
n2 <- 300 # before hidden has a mean shift
n3 <- 300
dd <- c(0, 0.5, 1, 1.5)
res <- lapply(1:num_sim, \(ii) {
  lst <- lapply(dd, \(d) {
    H1 <- rnorm(n1 + n2, mean = 0)
    H2 <- rnorm(n3, mean = 1)
    H <- c(H1, H2)
    X <- H + rnorm(n1 + n2 + n3)
    b1 <- 0.5
    b2 <- b1 + d
    Y1 <- H[1:n1] + X[1:n1] * b1 + rnorm(n1)
    Y2 <- H[(n1 + 1):(n1 + n2 + n3)] +
      X[(n1 + 1):(n1 + n2 + n3)] * b2 + rnorm(n2)
    Y <- c(Y1, Y2)
    est <- breakpoints(Y ~ X, h = 0.2, breaks = 2)$breakpoints
    if (anyNA(est)) {
      data.frame(d = d, est = NA, est_num = NA)
    } else {
      data.frame(d = d, est = sort(est), est_num = seq_len(length(est)))
    }
  })
  tmp <- do.call(rbind, lst)
  tmp$sim <- ii
  tmp
})
res_df <- do.call(rbind, res)
write.csv(res_df,
  "../output/rcp_with_nonstationary_hidden_res.csv",
  row.names = FALSE
)

make_plot(res_df,
  n_tot = n1 + n2 + n3,
  levels = dd,
  bin_w = 50,
  red_vdash = n1, gray_vdash = n1 + n2
)
ggsave("../output/plots/rcp_with_nonstationary_hidden_hist.pdf",
  width = 4, height = 4
)

# ---- make plot ----

# # Choose binwidth for numeric part
# bin_w <- 50
# first_bin_center <- bin_w / 2

# # Make NA bin center one binwidth to the left
# (na_center <- first_bin_center - bin_w)

# df_na <- res_df %>%
#   filter(is.na(est)) %>%
#   mutate(est = na_center)
# df_na$d <- factor(df_na$d, levels = dd)

# df_num <- res_df %>% filter(!is.na(est))
# df_num$d <- factor(df_num$d, levels = dd)

# ggplot() +
#   # Histogram for numeric
#   geom_histogram(
#     data = df_num,
#     aes(x = est, fill = d),
#     binwidth = bin_w,
#     boundary = 0, # start bins at 0
#     position = "dodge",
#     color = "gray"
#   ) +
#   # NA bar aligned with binwidth
#   geom_bar(
#     data = df_na,
#     aes(x = est, fill = d),
#     width = bin_w,
#     position = "dodge",
#     color = "gray"
#   ) +
#   geom_vline(xintercept = 400, linetype = 2, color = "red") +
#   coord_cartesian(xlim = c(na_center - bin_w / 2, 1000)) +
#   scale_x_continuous(
#     breaks = c(na_center, seq(0, 1000, 200)),
#     labels = c("No RCP", "", seq(0, 1000, 200)[-1])
#   ) +
#   labs(x = "Estimate", fill = bquote(beta[2] - beta[1])) +
#   scale_fill_brewer(palette = "Set2") +
#   theme_minimal() +
#   theme(legend.position = "top")
