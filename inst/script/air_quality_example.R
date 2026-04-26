library(dplyr)
library(ggplot2)
devtools::load_all()

# ---- helpers ----

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

esti_fun <- function(X, Y, I,
                     eval_gap = ceiling(length(I) * 0.01),
                     end_buff = max(ceiling(length(I) * 0.1), 10),
                     min_seg = ncol(X) + 1,
                     cond_sets = get_powerset(colnames(X)),
                     return_score = FALSE) {
  # eval_grid <- seq(end_buff, length(I)-end_buff, by = eval_gap)
  eval_grid <- seq(min(I) + end_buff, max(I)-end_buff, by = eval_gap)
  eval_grid_rel <- eval_grid - min(I) + 1
  vals <- sapply(eval_grid_rel, \(t) {
    cc_loss_niklas1(X, Y, t,
                    min(I), max(I),
                    min_seg = min_seg,
                    cond_sets = cond_sets)
  })
  if (return_score) {
    return(list(scores = vals,
                eval_grid = eval_grid,
                est = eval_grid[which.min(vals)]))
  } else {
    return(eval_grid[which.min(vals)])
  }
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

# ---- load and process data ----

df_orig <- read.csv2("beijing+multi+site+air+quality+data/PRSA_Data_20130301-20170228/PRSA_Data_Wanshouxigong_20130301-20170228.csv", header = TRUE, sep = ",")

df_orig <- df_orig %>%
  tidyr::unite(date, c("year", "month", "day"), sep = "-") %>%
  dplyr::mutate(date = zoo::as.Date(date, format = "%Y-%m-%d")) %>%
  dplyr::mutate(across(c(3:13, 15), as.numeric))
dim(df_orig)
df_orig <- df_orig %>%
  dplyr::filter(date < zoo::as.Date("2017-01-01", format = "%Y-%m-%d") &
                  date > zoo::as.Date("2014-01-01", format = "%Y-%m-%d")) %>%
  dplyr::mutate(season = get_season(date),
                weekday = weekdays(date),
                weekday = factor(weekday, levels = c("Monday", "Tuesday", "Wednesday",
                                                     "Thursday", "Friday", "Saturday", "Sunday")),
                year = lubridate::year(date),
                weeknum = lubridate::week(date))
df_orig$heating <- 0
df_orig$heating[(df_orig$date > zoo::as.Date("2013-11-15", format = "%Y-%m-%d")) &
                  (df_orig$date < zoo::as.Date("2014-03-15", format = "%Y-%m-%d"))] <- 1
df_orig$heating[(df_orig$date > zoo::as.Date("2014-11-15", format = "%Y-%m-%d")) &
                  (df_orig$date < zoo::as.Date("2015-03-15", format = "%Y-%m-%d"))] <- 1
df_orig$heating[(df_orig$date > zoo::as.Date("2015-11-15", format = "%Y-%m-%d")) &
                  (df_orig$date < zoo::as.Date("2016-03-15", format = "%Y-%m-%d"))] <- 1
df_orig$heating[(df_orig$date > zoo::as.Date("2016-11-15", format = "%Y-%m-%d")) &
                  (df_orig$date < zoo::as.Date("2017-03-15", format = "%Y-%m-%d"))] <- 1
plot(df_orig$date, df_orig$heating)

df_daily <- df_orig %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(across(c(3:12, 14), \(x) mean(x, na.rm = TRUE)),
                   season = unique(season),
                   weekday = unique(weekday),
                   year = unique(year),
                   heating = unique(heating),
                   weeknum = unique(weeknum)) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na()
plot(df_daily$date, df_daily$PM2.5)
plot(df_daily$date, df_daily$PM10)
plot(df_daily$date, df_daily$WSPM)

df_weekly <- df_orig %>%
  dplyr::group_by(year, weeknum) %>%
  dplyr::summarise(across(c(PM2.5, PM10, SO2, NO2, CO, O3, TEMP, PRES, DEWP, RAIN), \(x) mean(x, na.rm = TRUE)),
                   season = unique(season)) %>%
  dplyr::ungroup() %>%
  tidyr::drop_na()

df_use <- df_daily # choose which frequency to use

df_use %>%
  ggplot(aes(x = date, y = log(PM2.5), color = factor(heating))) +
  # ggplot(aes(x = date, y = log(PM2.5))) +
  geom_point(size = 1) +
  theme_bw() +
  labs(x = "date",
       y = "log(PM2.5)") +
  scale_color_discrete(name = "Heating") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.box.margin = margin(-3, 0, -3, 0),
        legend.position = "top",
        legend.margin = margin(-3, 0, -3, 0),
        text = element_text(size = 11, family = "serif"))

ggsave("marginal_lpm2p5_nocolor.pdf",
       width = 3.68, height = 2.75)

ggsave("marginal_lpm2p5.pdf",
       width = 3.68, height = 2.75)

X <- df_use %>%
  # dplyr::select(PM10, SO2, NO2, CO, O3, TEMP, PRES, DEWP, RAIN) %>%
  dplyr::select(PM10, SO2, NO2, CO, # O3,
                # TEMP, PRES, RAIN,
                # weekday, season,
                DEWP, WSPM, heating) %>%
  dplyr::mutate(PM10 = log(PM10),
                SO2 = log(SO2),
                NO2 = log(NO2),
                CO = log(CO)) %>%
  dplyr::mutate(across(.cols = 1:6,
                       .fn = \(x) as.numeric(scale(x)))) %>%
  dplyr::select(-PM10)
# weekday <- factor(df_use$weekday, levels = c("Monday", "Tuesday", "Wednesday",
#                                              "Thursday", "Friday", "Saturday", "Sunday"))
# season <- factor(df_use$season, levels = c("spring", "summer", "autumn", "winter"))
# trend_vars <- model.matrix(~ season + weekday)
# X <- cbind(X, trend_vars[,2:ncol(trend_vars)])
Y <- df_use %>%
  dplyr::select(PM2.5) %>%
  as.matrix() %>%
  log() %>%
  scale() %>%
  as.data.frame()

pairs(cbind(Y, X))

fm <- reformulate(colnames(X), colnames(Y))
mod0 <- lm(fm, data = cbind(X, Y))
summary(mod0)
plot(resid(mod0))
plot(pacf(resid(mod0)))

# # Marginal: CUSUM test
# CPAT::CUSUM.test(Y$PM2.5)$p.value
#
# # Regression: F test (chow test)
# fm <- reformulate(colnames(X), colnames(Y))
# strucchange::sctest(strucchange::Fstats(fm, data = cbind(X, Y)))$p.value
#
# # Our test
# PS <- get_powerset(colnames(X))
# test_wrap_chow(X = as.matrix(X),
#                Y = Y$PM2.5,
#                rel_grid = c(1/5, 2/5, 3/5, 4/5),
#                1:nrow(X), PS = PS, test_name = "chow2")

# ---- lag data ----

ks <- 0:1
dat <- lag_data(X, Y, ks = ks)
X_lag <- dat$X_lag
Y_lag <- dat$Y_lag
mod1 <- lm(PM2.5 ~ ., data = cbind(X_lag, Y_lag))
summary(mod1)
plot(resid(mod1))
plot(pacf(resid(mod1)))

# Marginal: CUSUM test
plot(df_use$date[-c(1:max(ks))], Y_lag$PM2.5)
CPAT::CUSUM.test(Y_lag$PM2.5)$p.value

# Regression: F test (chow test)
fm <- reformulate(colnames(X_lag), colnames(Y_lag))
strucchange::sctest(strucchange::Fstats(fm, data = cbind(X_lag, Y_lag)))$p.value

PS_orig <- get_powerset(colnames(X))
PS <- lapply(seq_along(PS_orig), \(ii) {
  c(PS_orig[[ii]], colnames(X_lag)[-(1:ncol(X))])
})
test_wrap_chow(X = as.matrix(X_lag), Y = Y_lag$PM2.5,
               rel_grid = c(1/5, 2/5, 3/5, 4/5),
               1:nrow(X_lag), PS = PS, test_name = "chow2")

ss <- c(120)
idx <- 1:nrow(X_lag)
eval_grid <- seq(from = max(ss) * 2 + 1,
                 to = nrow(X_lag[idx,]) - max(ss) * 2 - 1,
                 by = 5)
ccloss_vals <- sapply(eval_grid, \(t) cc_loss_niklas1(as.matrix(X_lag[idx,]),
                                                      Y_lag$PM2.5[idx],
                                                      t, 1, nrow(X_lag[idx,]),
                                                      min_seg = ss[1],
                                                      cond_sets = PS))
plot(df_use[-(1:max(ks)),][idx,]$date[eval_grid], ccloss_vals, type = 'l')

test_param <- list(PS = PS, alpha = 0.05) # for actual test
esti_param <- list(eval_gap = 1,
                   min_seg = 20,
                   end_buff = 41,
                   cond_sets = PS,
                   return_score = TRUE)
boundary_df <- seeded_intervals(nrow(X_lag), min_len = 120)
res <- seedBS_ccp(as.matrix(X_lag), Y_lag$PM2.5, boundary_df,
                  test_fun = test_wrap_chow, test_param = test_param,
                  esti_fun = esti_fun, esti_param = esti_param,
                  return_score = TRUE,
                  verbose = TRUE)
res$est
sort(df_use$date[res$est])

sapply(PS, \(S) {
  test_chow2(X_lag, Y_lag$PM2.5, 1:806, 807:nrow(X_lag), S)$p_value
})

# ---- without certain variables ----

pairs(cbind(Y, X))
df_use %>%
  ggplot(aes(x = date, y = DEWP, color = season)) +
  geom_point()

X_new <- X %>% dplyr::select(-SO2, -NO2, -CO, -heating)
ks <- 0:1
dat <- lag_data(X_new, Y, ks = ks)
X_lag <- as.data.frame(dat$X_lag)
Y_lag <- dat$Y_lag
mod1 <- lm(PM2.5 ~ ., data = cbind(X_lag, Y_lag))
summary(mod1)
plot(resid(mod1))
plot(pacf(resid(mod1)))

# Marginal: CUSUM test
CPAT::CUSUM.test(Y_lag$PM2.5)$p.value

# Regression: F test (chow test)
fm <- reformulate(colnames(X_lag), colnames(Y_lag))
strucchange::sctest(strucchange::Fstats(fm, data = cbind(X_lag, Y_lag)))$p.value

PS_orig <- get_powerset(colnames(X_new))
PS <- lapply(seq_along(PS_orig), \(ii) {
  c(PS_orig[[ii]], colnames(X_lag)[-(1:ncol(X_new))])
})
test_wrap_chow(X = as.matrix(X_lag), Y = Y_lag$PM2.5,
               rel_grid = c(1/5, 2/5, 3/5, 4/5),
               1:nrow(X_lag), PS = PS, test_name = "chow2")

# CCPs
test_param <- list(PS = PS, alpha = 0.05) # for actual test
esti_param <- list(eval_gap = 1,
                   min_seg = 30,
                   end_buff = 61,
                   cond_sets = PS,
                   return_score = TRUE)
boundary_df <- seeded_intervals(nrow(X_lag), min_len = 120)
res <- seedBS_ccp(as.matrix(X_lag), Y_lag$PM2.5, boundary_df,
                  test_fun = test_wrap_chow, test_param = test_param,
                  esti_fun = esti_fun, esti_param = esti_param,
                  return_score = TRUE,
                  verbose = TRUE)
sort(df_use$date[res$est])
# [1 lag] "2014-11-05" "2015-03-05" "2016-04-26"
# [2 lag] "2014-11-02" "2016-04-22"
# [3 lag] "2014-11-03" "2015-03-18" "2015-11-04" "2016-04-22"
# [4 lag] "2014-11-11" "2015-10-24" "2016-04-21"
# [5 lag] "2015-04-13" "2015-10-20" "2016-04-27"
# [6 lag] "2014-11-11" "2015-04-28" "2015-10-23" "2016-04-26"
# [7 lag] "2014-11-10" "2015-04-27" "2016-04-25"

# RCPs
bp_res <- strucchange::breakpoints(Y_lag$PM2.5 ~ as.matrix(X_lag), h = 0.1)
plot(bp_res)
sort(df_use$date[bp_res$extract.breaks(RSS.table = bp_res$RSS.table, 8)])

# RCPs (our loss)
test_param <- list(PS = colnames(X_lag), alpha = 0.05) # for actual test
esti_param <- list(eval_gap = 1,
                   min_seg = 30,
                   end_buff = 61,
                   cond_sets = colnames(X_lag),
                   return_score = TRUE)
boundary_df <- seeded_intervals(nrow(X_lag), min_len = 120)
res <- seedBS_ccp(as.matrix(X_lag), Y_lag$PM2.5, boundary_df,
                  test_fun = test_wrap_chow, test_param = test_param,
                  esti_fun = esti_fun, esti_param = esti_param,
                  return_score = TRUE,
                  verbose = TRUE)
sort(df_use$date[res$est])

# MCPs
sbs_obj <- wbs::wbs(Y_lag$PM2.5)
plot(sbs_obj)
sbs_res <- wbs::changepoints(sbs_obj)
sort(df_use$date[sbs_res$cpt.th[[1]]])

# ---- run all localizations and save results ----

ccp_all <- vector("list", 7)
rcp_all <- vector("list", 7)
X_new <- X %>% dplyr::select(-SO2, -NO2, -CO, -heating)
for (kmax in 1:7) {
  message("---- kmax: ", kmax, " ----")
  dat <- lag_data(X_new, Y, ks = 0:kmax)
  X_lag <- as.data.frame(dat$X_lag)
  Y_lag <- dat$Y_lag

  PS_orig <- get_powerset(colnames(X_new))
  PS <- lapply(seq_along(PS_orig), \(ii) {
    c(PS_orig[[ii]], colnames(X_lag)[-(1:ncol(X_new))])
  })

  # CCPs
  test_param <- list(PS = PS, alpha = 0.05) # for actual test
  esti_param <- list(eval_gap = 1,
                     min_seg = 30,
                     end_buff = 61,
                     cond_sets = PS,
                     return_score = TRUE)
  boundary_df <- seeded_intervals(nrow(X_lag), min_len = 120)
  res <- seedBS_ccp(as.matrix(X_lag), Y_lag$PM2.5, boundary_df,
                    test_fun = test_wrap_chow, test_param = test_param,
                    esti_fun = esti_fun, esti_param = esti_param,
                    return_score = TRUE,
                    verbose = TRUE)
  ccp_all[[kmax]] <- sort(df_use$date[res$est])

  # RCPs
  bp_res <- strucchange::breakpoints(Y_lag$PM2.5 ~ as.matrix(X_lag), h = 0.1)
  tmp <- summary(bp_res)
  opt_rss_ncp <- as.integer(names(which.min(tmp$RSS[1,])))
  rcp_all[[kmax]] <- sort(df_use$date[bp_res$extract.breaks(RSS.table = bp_res$RSS.table, opt_rss_ncp)])
}

# MCPs
set.seed(42) # wbs uses random intervals
sbs_obj <- wbs::wbs(Y_lag$PM2.5)
plot(sbs_obj)
sbs_res <- wbs::changepoints(sbs_obj)
mcp_all <- sort(df_use$date[sbs_res$cpt.th[[1]]])

saveRDS(list(ccp_all = ccp_all,
             rcp_all = rcp_all,
             mcp_all = mcp_all),
        "meteorological_ccp_rcp_mcp.rds")


# ---- (full covariates) going through different number of lags ----

kmax_all <- c(1,2,3,4,5,6,7)
# mcp_pvals <- rep(NA, length(kmax_all))
# rcp_pvals <- rep(NA, length(kmax_all))
# ccp_pvals <- rep(NA, length(kmax_all))

cp_pvals_lst <- lapply(kmax_all, \(kmax) {
  ks <- 0:kmax
  dat <- lag_data(X, Y, ks = ks)
  X_lag <- dat$X_lag
  Y_lag <- dat$Y_lag
  # mod1 <- lm(PM2.5 ~ ., data = cbind(X_lag, Y_lag))
  # summary(mod1)
  # plot(resid(mod1))
  # plot(pacf(resid(mod1)))

  # Marginal: CUSUM test
  # plot(df_use$date[-c(1:max(ks))], Y_lag$PM2.5)
  mcp_pval <- CPAT::CUSUM.test(Y_lag$PM2.5)$p.value

  # Regression: F test (chow test)
  fm <- reformulate(colnames(X_lag), colnames(Y_lag))
  rcp_pval <- strucchange::sctest(strucchange::Fstats(fm, data = cbind(X_lag, Y_lag)))$p.value

  # CCP
  PS_orig <- get_powerset(colnames(X))
  PS <- lapply(seq_along(PS_orig), \(ii) {
    c(PS_orig[[ii]], colnames(X_lag)[-(1:ncol(X))])
  })
  ccp_pval <- test_wrap_chow(X = as.matrix(X_lag), Y = Y_lag$PM2.5,
                 rel_grid = c(1/5, 2/5, 3/5, 4/5),
                 1:nrow(X_lag), PS = PS, test_name = "chow2")

  data.frame(mcp = mcp_pval, rcp = rcp_pval, ccp = ccp_pval)
})

cp_pvals_df1 <- do.call(rbind, cp_pvals_lst)
print(cp_pvals_df1, digits = 3)

# ---- (less covariates) going through different number of lags ----

kmax_all <- c(1,2,3,4,5,6,7)
X_new <- X %>% dplyr::select(-SO2, -NO2, -CO, -heating)

cp_pvals_lst <- lapply(kmax_all, \(kmax) {
  ks <- 0:kmax
  dat <- lag_data(X_new, Y, ks = ks)
  X_lag <- as.data.frame(dat$X_lag)
  Y_lag <- dat$Y_lag
  # mod1 <- lm(PM2.5 ~ ., data = cbind(X_lag, Y_lag))
  # summary(mod1)
  # plot(resid(mod1))
  # plot(pacf(resid(mod1)))

  # Marginal: CUSUM test
  # plot(df_use$date[-c(1:max(ks))], Y_lag$PM2.5)
  mcp_pval <- CPAT::CUSUM.test(Y_lag$PM2.5)$p.value

  # Regression: F test (chow test)
  fm <- reformulate(colnames(X_lag), colnames(Y_lag))
  rcp_pval <- strucchange::sctest(strucchange::Fstats(fm, data = cbind(X_lag, Y_lag)))$p.value

  # CCP
  PS_orig <- get_powerset(colnames(X_new))
  PS <- lapply(seq_along(PS_orig), \(ii) {
    c(PS_orig[[ii]], colnames(X_lag)[-(1:ncol(X_new))])
  })
  ccp_pval <- test_wrap_chow(X = as.matrix(X_lag), Y = Y_lag$PM2.5,
                             rel_grid = c(1/5, 2/5, 3/5, 4/5),
                             1:nrow(X_lag), PS = PS, test_name = "chow2")

  data.frame(mcp = mcp_pval, rcp = rcp_pval, ccp = ccp_pval)
})

cp_pvals_df2 <- do.call(rbind, cp_pvals_lst)
print(cp_pvals_df2, digits = 3)
cp_pvals_df1$model <- "All covariates"
cp_pvals_df2$model <- "Only meterological"
cp_pvals_df <- rbind(cp_pvals_df1, cp_pvals_df2) %>%
  dplyr::mutate(nlag = rep(1:7, 2),
                nlag = factor(nlag)) %>%
  tidyr::pivot_longer(1:3, names_to = "cp", values_to = "pval") %>%
  dplyr::mutate(cp = toupper(cp),
                cp = factor(cp, levels = c("MCP", "RCP", "CCP")))

cp_pvals_df %>%
  ggplot(aes(x = nlag, y = log(pval), color = cp, group = interaction(cp, model))) +
  geom_point() +
  geom_line() +
  facet_grid(model ~., scales = "free_y") +
  geom_hline(yintercept = log(0.05), color = 'red', linetype = 2) +
  labs(x = "Number of lags",
       y = "log(p-value)") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.box.margin = margin(-3, 0, -3, 0),
        legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(-3, 0, -3, 0),
        text = element_text(size = 11, family = "serif"))

ggsave("realdata_airquality_pvals.pdf",
       width = 3.68, height = 2.75)
