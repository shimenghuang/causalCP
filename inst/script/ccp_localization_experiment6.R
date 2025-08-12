library(future.apply)
plan(multisession)

test_wrap_chow <- function(X, Y, I, PS = get_powerset(colnames(X)),
                           rel_grid = c(1 / 3, 2 / 3), test_name = "chow2",
                           intercept = TRUE, rev_one_rest = FALSE, two_dir = FALSE,
                           save_resids = FALSE, alpha = 0.05) {
  res <- one_vs_rest_test(
    X[I, , drop = FALSE], Y[I], PS, rel_grid, test_name,
    intercept, rev_one_rest, two_dir, save_resids
  )
  min(sapply(res$pvals_all, max)) # * length(rel_grid)
}

esti_fun <- function(X, Y, I,
                     eval_gap = ceiling(length(I) * 0.01),
                     end_buff = max(ceiling(length(I) * 0.1), 10),
                     min_seg = ncol(X) + 1,
                     cond_sets = get_powerset(colnames(X)),
                     return_score = FALSE) {
  # eval_grid <- seq(end_buff, length(I)-end_buff, by = eval_gap)
  eval_grid <- seq(min(I) + end_buff, max(I) - end_buff, by = eval_gap)
  eval_grid_rel <- eval_grid - min(I) + 1
  vals <- sapply(eval_grid_rel, \(t) {
    cc_loss(X, Y, t,
      min(I), max(I),
      min_seg = min_seg,
      cond_sets = cond_sets
    )
  })
  if (return_score) {
    return(list(
      scores = vals,
      eval_grid = eval_grid,
      est = eval_grid[which.min(vals)]
    ))
  } else {
    return(eval_grid[which.min(vals)])
  }
}

#' Higher dimensional example with one causal change point and one non-causal change point
#'
dgp_hd <- function(n1, n2, n3, d, idx_y = NULL,
                   num_parents = NULL,
                   num_children = NULL,
                   num_connect = NULL,
                   multi_factor = c(2, 2),
                   seed_num = NULL) {
  if (!is.null(seed_num)) set.seed(seed_num)

  n <- n1 + n2 + n3
  if (is.null(idx_y)) idx_y <- ceiling(median(1:d))

  # connectivity
  B <- matrix(0, d, d)
  tot_connect <- sum(upper.tri(B[-idx_y, -idx_y]))
  if (is.null(num_connect)) num_connect <- floor(tot_connect * 0.5)
  message("number of additional edges: ", num_connect)
  B[-idx_y, -idx_y][upper.tri(B)[-idx_y, -idx_y]] <- sample(c(
    rep(0, tot_connect - num_connect),
    rep(1, num_connect)
  ))
  if (!is.null(num_parents)) {
    parents <- sample(1:idx_y, num_parents, replace = FALSE)
    B[parents, idx_y] <- 1
    B[-parents, idx_y] <- 0
  }
  if (!is.null(num_children)) {
    children <- sample(1:idx_y, num_children, replace = FALSE)
    B[idx_y, children] <- 1
    B[idx_y, -children] <- 0
  }

  Id <- diag(d)
  eps <- matrix(rnorm(n * d, 0, 0.5), n, d)
  # eps[(n1+n2+1):n,-idx_y] <- eps[(n1+n2+1):n,-idx_y] * 0.5
  B1 <- B * 0.5
  B1[, idx_y] <- B1[, idx_y] * 1
  B2 <- B1
  B2[, idx_y] <- B2[, idx_y] * multi_factor[1]
  B3 <- B2
  B3[, -idx_y] <- B3[, -idx_y] * multi_factor[2]
  X1 <- eps[1:n1, ] %*% solve(Id - B1)
  X2 <- eps[(n1 + 1):(n1 + n2), ] %*% solve(Id - B2)
  X3 <- eps[(n1 + n2 + 1):n, ] %*% solve(Id - B3)
  Y <- c(X1[, idx_y], X2[, idx_y], X3[, idx_y])
  X <- rbind(X1[, -idx_y], X2[, -idx_y], X3[, -idx_y])
  colnames(X) <- seq_len(ncol(X))

  beta_star <- cbind(B1[, idx_y], B2[, idx_y], B3[, idx_y])

  return(list(
    X = X,
    Y = Y,
    B = B,
    parents = which(B[-idx_y, idx_y] != 0),
    beta_star = beta_star
  ))
}

# # ---- one experiment ----

# dat <- dgp_hd(300, 300, 300,
#   d = 1001,
#   num_parents = 2, num_children = 1, num_connect = 5,
#   multi_factor = c(4, 4),
#   seed_num = NULL
# )
# # dat$B
# which(dat$B != 0, arr.ind = TRUE)
# dat$parents
# dat$beta_star[dat$parents, ]
# dim(dat$X)
# plot(dat$Y)

# fit <- glmnet::glmnet(dat$X, dat$Y, alpha = 1)
# cvfit <- glmnet::cv.glmnet(dat$X, dat$Y, nfolds = 10)
# (selected_vars <- which(as.numeric(coef(cvfit, s = "lambda.1se"))[-1] != 0))
# selected_vars <- as.character(selected_vars)
# dat$parents

# # PS <- get_powerset(colnames(dat$X))
# PS <- get_powerset(selected_vars)
# # test_wrap_chow(dat$X, dat$Y, 1:nrow(dat$X), PS, rel_grid = c(1 / 2))
# test_wrap_chow(dat$X, dat$Y, 1:nrow(dat$X), PS, rel_grid = c(1 / 3, 2 / 3))
# test_wrap_chow(dat$X, dat$Y, 1:nrow(dat$X), PS, rel_grid = c(1 / 4, 2 / 4, 3 / 4))

# ss <- c(100)
# idx <- 1:nrow(dat$X)
# eval_grid <- seq(
#   from = max(ss) * 2 + 1,
#   to = nrow(dat$X[idx, ]) - max(ss) * 2 - 1,
#   by = 10
# )
# ccloss_vals <- sapply(eval_grid, \(t) cc_loss(as.matrix(dat$X[idx, ]),
#   dat$Y[idx],
#   t, 1, nrow(dat$X[idx, ]),
#   min_seg = ss[1],
#   cond_sets = PS
# ))
# plot(eval_grid, ccloss_vals, type = "l")
# eval_grid[which.min(ccloss_vals)]

# test_param <- list(PS = PS, alpha = 0.05) # for actual test
# esti_param <- list(
#   eval_gap = 10,
#   min_seg = 49,
#   end_buff = 100,
#   cond_sets = PS,
#   return_score = TRUE
# )
# boundary_df <- seeded_intervals(nrow(dat$X[idx, ]), min_len = 400)
# res <- seedBS_ccp(as.matrix(dat$X[idx, ]),
#   dat$Y[idx],
#   boundary_df,
#   test_fun = test_wrap_chow, test_param = test_param,
#   esti_fun = esti_fun, esti_param = esti_param,
#   return_score = TRUE,
#   verbose = TRUE
# )
# res$est
# sort(res$est)

# --- repetitions ----

# strucchange::breakpoints(dat$Y ~ dat$X)
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
# multi_num <- 4
multi_num <- as.numeric(args[1])
message("multi_num: ", multi_num)
n123 <- c(300, 300, 300)
num_connect <- 5
reps <- 1:200
# seeds <- future_lapply(seq_along(reps),
#                        FUN = function(x) .Random.seed,
#                        future.chunk.size = Inf, future.seed = 42L)
res_all <- future_lapply(reps, \(rr) {
  devtools::load_all()

  # seed_num <- .Random.seed[-1][rr]
  # message("seed", seed_num)

  max_regen <- 10
  num_regen <- 0
  while (num_regen < max_regen) {
    num_regen <- num_regen + 1
    dat <- tryCatch(
      {
        dgp_hd(n123[1], n123[2], n123[3],
          d = 1001,
          num_parents = 2, num_children = 1, num_connect = num_connect,
          multi_factor = rep(multi_num, 2),
          seed_num = NULL
        )
      },
      error = function(e) {
        message("Error: ", e$message)
        NULL
      }
    )
    if (!is.null(dat)) {
      message("dgp run ", num_regen, " time(s)")
      break
    }
  }

  message("dat$X: ", nrow(dat$X), " x ", ncol(dat$X))
  message("dat$Y: ", length(dat$Y))

  cvfit <- glmnet::cv.glmnet(dat$X, dat$Y, alpha = 1, nfolds = 10, dfmax = 10)
  (selected_vars <- which(as.numeric(coef(cvfit, s = "lambda.1se"))[-1] != 0))
  selected_vars <- as.character(selected_vars)
  message("number of selected vars: ", length(selected_vars))
  message("selected vars: ", paste(selected_vars, collapse = " "))

  # # if there are more than 10 varaibles selected by Lasso, then move on
  # if (length(selected_vars) > 10) {
  #   tibble::tibble(rep = rr, pval = NA, est1 = NA, est2 = NA)
  #   next
  # }

  # get powerset of Lasso selected variables
  PS <- get_powerset(selected_vars)

  # test
  pval <- test_wrap_chow(dat$X, dat$Y, 1:nrow(dat$X), PS,
    rel_grid = c(1 / 4, 2 / 4, 3 / 4)
  )

  # minimizer of loss
  ss <- c(100)
  idx <- 1:nrow(dat$X)
  eval_grid <- seq(
    from = max(ss) * 2 + 1,
    to = nrow(dat$X[idx, ]) - max(ss) * 2 - 1,
    by = 10
  )
  ccloss_vals <- sapply(eval_grid, \(t) cc_loss(
    as.matrix(dat$X[idx, ]),
    dat$Y[idx],
    t, 1, nrow(dat$X[idx, ]),
    min_seg = ss[1],
    cond_sets = PS
  ))
  est1 <- eval_grid[which.min(ccloss_vals)]

  # seeded binseg
  test_param <- list(PS = PS, alpha = 0.05) # for actual test
  esti_param <- list(
    eval_gap = 10,
    min_seg = 49,
    end_buff = 100,
    cond_sets = PS,
    return_score = TRUE
  )
  boundary_df <- seeded_intervals(nrow(dat$X[idx, ]), min_len = 400)
  res <- seedBS_ccp(as.matrix(dat$X[idx, ]),
    dat$Y[idx],
    boundary_df,
    test_fun = test_wrap_chow, test_param = test_param,
    esti_fun = esti_fun, esti_param = esti_param,
    return_score = TRUE,
    verbose = TRUE
  )
  est2 <- res$est

  # seeded binseg + prune
  tibble::tibble(rep = rr, pval = pval, est1 = est1, est2 = est2)
}, future.seed = TRUE)

res_df <- do.call(rbind, res_all)
saveRDS(res_df, paste0(
  "../output/",
  "hd_ccp_nccp_strength", multi_num,
  "_num_connect", num_connect,
  "_ns", paste(n123, collapse = "_"), ".rds"
))

# ---- plot results ----

hd_ccp_nccp_strength3 <- readRDS("../output/hd_ccp_nccp_strength3_num_connect5_ns300_300_300.rds")
hd_ccp_nccp_strength3.5 <- readRDS("../output/hd_ccp_nccp_strength3.5_num_connect5_ns300_300_300.rds")
hd_ccp_nccp_strength4 <- readRDS("../output/hd_ccp_nccp_strength4_num_connect5_ns300_300_300.rds")
hd_ccp_nccp_strength2.5$strength <- 2.5
hd_ccp_nccp_strength3$strength <- 3
hd_ccp_nccp_strength3.5$strength <- 3.5
hd_ccp_nccp_strength4$strength <- 4
res_comb <- rbind(
  hd_ccp_nccp_strength3,
  hd_ccp_nccp_strength3.5,
  hd_ccp_nccp_strength4
)

res_comb_long <- res_comb %>%
  dplyr::mutate(strength = factor(as.integer(factor(strength)))) %>%
  tidyr::pivot_longer(3:4, names_to = "method", values_to = "est") %>%
  dplyr::group_by(rep, strength, method) %>%
  dplyr::distinct() %>%
  dplyr::group_by(rep, strength, method) %>%
  dplyr::arrange(est, .by_group = TRUE) %>%
  dplyr::mutate(num_est = factor(1:n()))

res_comb_long %>%
  dplyr::filter(method == "est1") %>%
  ggplot(aes(x = strength, y = pval, group = strength)) +
  ggbeeswarm::geom_quasirandom(width = 0.2, size = 0.75, alpha = 0.5) +
  geom_boxplot(alpha = 0.1, width = 0.45) +
  geom_hline(yintercept = 0.05, linetype = 2, color = "red") +
  labs(y = "p-value") +
  theme_bw() +
  theme(text = element_text(family = "serif", size = 11))

ggsave("../output/plots/highdimensional_pvals.pdf",
  width = 4, height = 2.5
)

res_comb_long %>%
  ggplot(aes(x = est, y = strength, group = method, color = method, shape = num_est)) +
  ggbeeswarm::geom_quasirandom(size = 2, width = 0.1, dodge.width = 0.6, alpha = 0.75) +
  geom_vline(xintercept = 301, linetype = 2, color = "red") +
  scale_color_discrete(name = "Method", labels = c("LossCS", "LossCS-SeedBS")) +
  scale_shape_manual(name = "CCP estimate", values = c(20, 4), labels = c("1st est.", "2nd est.")) +
  labs(x = "estimated CCP") +
  theme_bw() +
  theme(
    text = element_text(family = "serif", size = 11),
    legend.position = "top",
    legend.box.margin = margin(-3, 0, -3, 0),
    legend.margin = margin(-3, 0, -3, 0),
    legend.box = "vertical",
    legend.spacing.y = unit(0.1, "cm")
  ) +
  guides(
    colour = guide_legend(keywidth = 0.1, nrow = 1),
    shape = guide_legend(keywidth = 0.1, nrow = 1),
  )

ggsave("../output/plots/highdimensional_estimates.pdf",
  width = 4, height = 2.5
)
