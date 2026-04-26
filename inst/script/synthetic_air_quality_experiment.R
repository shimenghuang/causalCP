#' Specifying lags in formulas
#'
#' @param x Numeric vector or time series
#' @param k Number of lags
#'
#' @importFrom data.table shift
#' @export
lags <- function(x, k = 1) {
  ret <- data.table::shift(x, k, fill = NA)
  if (is.list(ret))
    ret <- do.call("cbind", ret)
  ret
}

#' @param Y should be a matrix of dimension n x 1
lag_data <- function(X, Y, ks) {
  orig_names <- c(colnames(X), names(Y))
  X_lag_lst <- lapply(ks, \(k) {
    if (k == 0) {
      X_lagk <- cbind(apply(X, 2, lags, k = k))
    } else {
      X_lagk <- cbind(apply(X, 2, lags, k = k), lags(Y, k = k))
    }
  })
  X_lag <- do.call(cbind, X_lag_lst)[-(1:max(ks)),]
  Y_lag <- Y[-(1:max(ks)),,drop=FALSE]
  if (length(setdiff(ks, 0)) > 0) {
    colnames(X_lag) <- c(colnames(X),
                         paste0(rep(orig_names, sum(ks != 0)),
                                rep(setdiff(ks, 0), each = length(orig_names))))
  } else {
    colnames(X_lag) <- colnames(X)
  }
  names(Y_lag) <- names(Y)

  return(list(X_lag = X_lag,
              Y_lag = Y_lag))
}

test_wrap_seqicp <- function(X, Y, I, PS, test, link = sum, grid = NULL, alpha = 0.05) {
  if (is.null(grid)) {
    grid <- c(0, round(nrow(X)/3), round(2 * nrow(X)/3), nrow(X))
  }
  pvals <- sapply(PS, \(S) {
    seqICP::seqICP.s(X, Y, S, test = test,
                     par.test = list(grid = grid,
                                     complements = TRUE,
                                     link = link, # max
                                     alpha = alpha,
                                     B = 999,
                                     permutation = FALSE))$p.value
  })
  return(max(pvals))
}

test_wrap_chow <- function(X, Y, I, PS = get_powerset(colnames(X)),
                           rel_grid = c(1/3, 2/3), test_name = "chow2",
                           intercept = TRUE, rev_one_rest = FALSE, two_dir = FALSE,
                           save_resids = FALSE, alpha = 0.05) {
  res <- one_vs_rest_test(X[I,,drop=FALSE], Y[I], PS, rel_grid, test_name,
                          intercept, rev_one_rest, two_dir, save_resids)
  min(sapply(res$pvals_all, max)) # * length(rel_grid)
}

get_season <- function(dates) {
  spring <- 3L
  summer <- 6L
  autumn <- 9L
  winter <- 12L
  seasons <- sapply(dates, \(date) {
    mon <- as.numeric(format(date, "%m"))
    if (is.na(mon)) stop("month should be between 1 and 12")
    if (mon < spring | mon >= winter) season <- "winter"
    else if (mon >= spring & mon < summer) season <- "spring"
    else if (mon >= summer & mon < autumn) season <- "summer"
    else if (mon >= autumn & mon < winter) season <- "autumn"
    season
  })
  return(factor(seasons, levels = c("spring", "summer", "autumn", "winter")))
}

# ----

set.seed(624)
dates <- seq.Date(from = as.Date("2011-03-15"),
                  to = as.Date("2018-03-14"),
                  by = "day")
years <- lubridate::year(dates)
months <- lubridate::month(dates)
days <- lubridate::day(dates)
seasons <- get_season(dates)
df_all <- data.frame(date = dates, year = years, month = months, day = days, season = seasons)
df_all <- df_all %>%
  dplyr::mutate(heating = 0,
                heating = dplyr::if_else((month == 3) & (day < 15), 1, heating),
                heating = dplyr::if_else((month == 11) & (day >= 15), 1, heating),
                heating = dplyr::if_else((month >= 1) & (month < 3), 1, heating),
                heating = dplyr::if_else((month > 11) & (month <= 12), 1, heating),
                # fuel = 1,
                fuel = dplyr::if_else(year >= 2015, 0.5, 1),
                # dewp = 0.5 + arima.sim(n(), model = list(ar = c(0.5, 0.1), sd = sqrt(0.1))),
                dewp = arima.sim(n = n(),
                                 model = list(ar = c(0.8)), # , ma = c(-0.2, 0.25)
                                 rand.gen = function(n, ...) 0.2 * rt(n, df = 10)),
                wspm = 1 + arima.sim(n(),
                                     # model = list(ar = c(0.3), sd = sqrt(0.1)),
                                     model = list(ar = c(0.3)),
                                     rand.gen = function(n, ...) 0.5 * rt(n, df = 10)),
                lpm2p5 = 2 * fuel * heating +
                  fuel * arima.sim(n(),
                                   # model = list(ar = c(0.5), sd = sqrt(0.1)),
                                   model = list(ar = c(0.5)),
                                   rand.gen = function(n, ...) rt(n, df = 10)) +
                  dewp - 0.5 * dplyr::lag(as.numeric(dewp)) - wspm - 0.3 * dplyr::lag(as.numeric(wspm)),
                lpm10 = dplyr::if_else((year %in% c(2011, 2012, 2016, 2017)) & (month %in% c(9, 10, 11)),
                                       rnorm(n(), 2, 0.5), rnorm(n(), 0.5, sqrt(0.5))) + lpm2p5)
df_all <- df_all[-1,]
head(df_all)
head(as_tibble(df_all[,c(6, 8, 9, 11)]))
plot(df_all$date, df_all$lpm2p5)
plot(df_all$date, exp(df_all$lpm2p5))
plot(df_all$date, df_all$lpm10)

ks <- 0:1
dat_lag_lst <- lag_data(as_tibble(df_all[,c(6, 8, 9, 11)]),
                        data.frame(lpm2p5 = df_all$lpm2p5), ks = ks)
X_lag <- dat_lag_lst$X_lag
Y_lag <- dat_lag_lst$Y_lag
dat_lag <- cbind(X_lag, Y_lag)
plot(resid(lm(Y_lag$lpm2p5 ~ X_lag)))

# test_hansen2000(X_lag, Y_lag$lpm2p5, B = 10)

# ---- detection ----

# Marginal: CUSUM test
CPAT::CUSUM.test(Y_lag$lpm2p5)$p.value
# [1] 8.717863e-05 # fuel \neq 1
# [1] 5.28067e-07 # fuel = 1

# Regression: F test (chow test)
fm <- reformulate(colnames(X_lag), colnames(Y_lag))
strucchange::sctest(strucchange::Fstats(fm, data = dat_lag))$p.value
test_wrap_chow(X_lag, Y_lag$lpm2p5, I = 1:nrow(X_lag),
               PS = colnames(X_lag), rel_grid = c(1/3, 2/3))
# [1] 0 # fuel \neq 1
# [1] 2.74679e-06 # fuel = 1

# Oracle regression: F test (chow test)
fm1 <- reformulate(colnames(X_lag)[!grepl("lpm10", colnames(X_lag))], colnames(Y_lag))
strucchange::sctest(strucchange::Fstats(fm1, data = dat_lag))$p.value
test_wrap_chow(X_lag, Y_lag$lpm2p5, I = 1:nrow(X_lag),
               PS = colnames(X_lag)[!grepl("lpm10", colnames(X_lag))], rel_grid = c(1/3, 2/3))
# [1] 5.259429e-11
# 0.5405055 # fuel = 1

# Causal: our test
PS_orig <- get_powerset(c("dewp",
                          "wspm",
                          "heating",
                          "lpm10")) # [,-c(1,4,5)]
PS <- lapply(seq_along(PS_orig), \(ii) {
  c(PS_orig[[ii]], colnames(X_lag)[-c(1:4)])
})
# PS <- get_powerset(colnames(X_lag))
test_wrap_chow(X = as.matrix(X_lag), Y = Y_lag$lpm2p5,
               rel_grid = c(1/3, 2/3),
               1:nrow(X_lag), PS = PS, test_name = "chow2")
# [1] 1.132638e-05
# [1] 0.6593862 # fuel = 1

# ---- localization ----

# Marginal: sbs
sbs_obj <- wbs::sbs(Y_lag$lpm2p5)
plot(sbs_obj)
sbs_res <- wbs::changepoints(sbs_obj)
sort(df_all$date[-1][sbs_res$cpt.th[[1]]])
# [1] "2011-05-03" "2011-11-15" "2012-03-15" "2012-04-07" "2012-06-05" "2012-07-10" "2012-11-11"
# [8] "2013-03-07" "2013-04-26" "2013-11-12" "2014-03-19" "2014-08-27" "2014-09-30" "2014-11-11"
# [15] "2015-03-24" "2015-11-20" "2016-03-22" "2016-11-10" "2017-03-28" "2017-11-14"

# [1] "2011-05-03" "2011-11-15" "2012-03-15" "2012-11-11" "2013-03-07" "2013-11-12" "2014-03-10"
# [8] "2014-09-30" "2014-11-11" "2015-03-24" "2015-11-10" "2016-03-22" "2016-04-08" "2016-12-11"
# [15] "2017-03-28" "2017-11-14"

# Marginal: wbs
wbs_obj <- wbs::wbs(Y_lag$lpm2p5)
plot(wbs_obj)
wbs_res <- wbs::changepoints(wbs_obj)
sort(df_all$date[-1][wbs_res$cpt.th[[1]]])
# [1] "2011-05-03" "2011-11-18" "2012-03-15" "2012-03-22" "2012-04-08" "2012-06-05" "2012-07-10"
# [8] "2012-11-14" "2012-12-16" "2013-03-07" "2013-04-26" "2013-09-07" "2013-09-14" "2013-11-14"
# [15] "2014-02-05" "2014-03-19" "2014-05-10" "2014-05-17" "2014-07-11" "2014-08-27" "2014-09-30"
# [22] "2014-11-14" "2015-01-02" "2015-03-24" "2015-11-20" "2016-03-23" "2016-04-08" "2016-12-11"
# [29] "2017-03-19" "2017-07-09" "2017-11-14"

# Regression: bai-perron with dynamic programming
fm <- reformulate(colnames(X_lag), colnames(Y_lag))
(bp_res <- strucchange::breakpoints(fm, data = dat_lag, breaks = 3))
sort(df_all$date[-1][bp_res$breakpoints])
# [1] "2012-11-30" "2014-12-29"
# [1] "2012-11-30" "2016-08-29" # fuel = 1
plot(bp_res)

# Oracle regression: bai-perron with dynamic programming
fm1 <- reformulate(colnames(X_lag)[!grepl("lpm10", colnames(X_lag))], colnames(Y_lag))
(bp_res1 <- strucchange::breakpoints(fm1, data = dat_lag))
sort(df_all$date[-1][bp_res1$breakpoints])
# [1] "2014-12-09"
# [1] NA # fuel = 1
plot(bp_res1)

# Oracle regression: sbs on residual of oracle regression
sbs_obj1 <- wbs::sbs(resid(lm(fm1, data = dat_lag)))
plot(sbs_obj1)
sbs_res1 <- wbs::changepoints(sbs_obj1)
sort(df_all$date[-1][sbs_res1$cpt.th[[1]]])

# Oracle regression: wbs on residual of oracle regression
wbs_obj1 <- wbs::wbs(resid(lm(fm1, data = dat_lag)))
plot(wbs_obj1)
wbs_res1 <- wbs::changepoints(wbs_obj1, penalty="mbic.penalty")
sort(df_all$date[-1][wbs_res1$cpt.ic[[1]]])
sort(df_all$date[-1][wbs_res1$cpt.th[[1]]])

# CausalCP: candidates
bp_res$breakpoints
test_wrap_chow(X_lag, Y_lag$lpm2p5, 1:1384, PS = PS, rel_grid = c(1/3, 2/3))
test_wrap_chow(X_lag, Y_lag$lpm2p5, 626:nrow(X_lag), PS = PS, rel_grid = c(1/3, 2/3))

# CausalCP: LossCS
ss <- c(100)
idx <- 1:nrow(X_lag)
eval_grid <- seq(from = max(ss) * 2 + 1,
                 to = nrow(X_lag[idx,]) - max(ss) * 2 - 1,
                 by = 1)
ccloss_vals <- sapply(eval_grid, \(t) cc_loss_niklas1(as.matrix(X_lag[idx,]),
                                                      Y_lag$lpm2p5[idx],
                                                      t, 1, nrow(X_lag[idx,]),
                                                      min_seg = ss[1],
                                                      cond_sets = PS))
plot(df_all[-(1:max(ks)),][idx,]$date[eval_grid], ccloss_vals, type = 'l')
df_all$date[-1][eval_grid[which.min(ccloss_vals)]]
# [1] "2015-02-07"
# [1] # fuel = 1

# ---- old code ----

(efp_res <- strucchange::efp(fm, data = dat_lag, h = 0.2, type = "OLS-CUSUM"))
plot(efp_res)
width <- unique(strucchange::boundary(efp_res))
efp_df <- data.frame(date = df_all$date[-1],
                     efp = efp_res$process[-1],
                     lwr = -width,
                     upr = width)
efp_df %>%
  ggplot(aes(x = date, y = efp)) +
  geom_line() +
  geom_hline(yintercept = c(-width, width), color = 'red', linetype = 2) +
  theme_bw()

strucchange::sctest(strucchange::Fstats(fm, data = dat_lag))
strucchange::breakpoints(strucchange::Fstats(fm, data = dat_lag))
PS_orig <- get_powerset(c("dewp",
                          "wspm",
                          "heating",
                          "lpm10")) # [,-c(1,4,5)]
PS <- lapply(seq_along(PS_orig), \(ii) {
  c(PS_orig[[ii]], colnames(X_lag)[-c(1:4)])
})
# PS <- get_powerset(colnames(X_lag))
test_wrap_chow(X = as.matrix(X_lag), Y = Y_lag$lpm2p5,
               rel_grid = c(1/3, 2/3),
               1:nrow(X_lag), PS = PS, test_name = "chow2")

ss <- c(120)
idx <- 1:nrow(X_lag)
eval_grid <- seq(from = max(ss) * 2 + 1,
                 to = nrow(X_lag[idx,]) - max(ss) * 2 - 1,
                 by = 1)
# eval_grid <- eval_grid[50:length(eval_grid)]
ccloss_vals <- sapply(eval_grid, \(t) cc_loss_niklas1(as.matrix(X_lag[idx,]),
                                                      Y_lag$lpm2p5[idx],
                                                      t, 1, nrow(X_lag[idx,]),
                                                      min_seg = ss[1],
                                                      cond_sets = PS))
plot(df_all[-(1:max(ks)),][idx,]$date[eval_grid], ccloss_vals, type = 'l')
# plot(df_daily$date[-(1:3)][eval_grid], ccloss_vals, type = 'l')
df_all$date[-1][eval_grid[which.min(ccloss_vals)]]

test_param <- list(PS = PS, alpha = 0.05) # for actual test
esti_param <- list(eval_gap = 1,
                   min_seg = 30,
                   end_buff = 60,
                   cond_sets = PS,
                   return_score = TRUE)
boundary_df <- seeded_intervals(nrow(X_lag), min_len = 366)
res <- seedBS_ccp(X_lag, Y_lag$lpm2p5, boundary_df,
                  test_fun = test_wrap_chow, test_param = test_param,
                  esti_fun = esti_fun, esti_param = esti_param,
                  return_score = TRUE,
                  verbose = TRUE)
res$est
df_all$date[-1][res$est]

idx <- 1401:1924
test_wrap_chow(X = as.matrix(X_lag)[idx,], Y = Y_lag$lpm2p5[idx],
               rel_grid = c(1/3, 2/3),
               1:nrow(X_lag[idx,]), PS = PS, test_name = "chow2")

saveRDS(df_all, file = "synthetic_air_quality_data2.rds")
