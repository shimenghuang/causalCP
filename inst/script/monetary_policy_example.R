library(ggplot2)

pkg_dir <- "/Users/hrt620/Documents/ku_projects/causalCPD_proj/causalCPD"

lag_vars <- function(X, Y, shifts = c(0L, 1L)) {
  X_new <- lapply(shifts, \(shift) {
    Xk <- dplyr::lag(X, n = shift)
    if (shift > 0) {
      Yk <- dplyr::lag(Y, n = shift)
      Xk <- cbind(Xk, Yk)
      colnames(Xk) <- c(paste0(colnames(X), "_lag", shift), paste0("Y_lag", shift))
    } else {
      colnames(Xk) <- paste0(colnames(X), "_lag", shift)
    }
    Xk
  }) %>%
    dplyr::bind_cols()
  X_new <- X_new[(max(shifts) + 1):nrow(X),]
  col_names <- colnames(X_new)
  X_new <- as.matrix(X_new)
  colnames(X_new) <- col_names
  Y_new <- Y[(max(shifts) + 1):length(Y)]
  return(list(X = X_new, Y = Y_new))
}

drop_vars <- function(X, names) {
  keep_names <- setdiff(colnames(X), names)
  return(X[,keep_names])
}

esti_fun <- function(X, Y, I,
                     eval_gap = 1,
                     end_buff = max(ceiling(length(I) * 0.1), 10),
                     cond_sets = get_powerset(colnames(X)),
                     min_seg = 30) {
  # Note: for cc_loss t is relative to the interval
  eval_grid <- seq(end_buff, length(I)-end_buff, by = eval_gap)
  vals <- sapply(eval_grid, \(t) cc_loss_niklas1(X, Y, t, min(I), max(I),
                                                 min_seg = min_seg,
                                                 cond_sets = cond_sets))
  return(min(I) + eval_grid[which.min(vals)])
}

test_wrap <- function(X, Y, I, PS = get_powerset(colnames(X)),
                      rel_grid = c(1/2), test_name = "chow1",
                      intercept = TRUE, rev_one_rest = FALSE, two_dir = FALSE,
                      save_resids = FALSE, alpha = 0.05) {
  res <- one_vs_rest_test(X[I,,drop=FALSE], Y[I], PS, rel_grid, test_name,
                          intercept, rev_one_rest, two_dir, save_resids)
  min(sapply(res$pvals_all, max)) # * length(rel_grid)
}

test_wrap <- function(X, Y, I, min_seg = 30) {
  length_out <- max(length(I)/min_seg, 3)
  grid_pts <- seq(from = 0, to = length(I), length.out = length_out)
  res <- seqICP(as.matrix(X[I,,drop = FALSE]), Y[I], test = "variance",
                par.test = list(grid = grid_pts,
                                complements = TRUE,
                                link = sum,
                                alpha = 0.05,
                                B = 100,
                                permutation = FALSE),
                par.model = list(pknown = TRUE, p = 3),
                stopIfEmpty = FALSE, silent = TRUE,
                model = "ar")
  return(min(res$test.results$p.value))
}

# ---- load data ----

load(paste0(pkg_dir, "/inst/data/monetarypolicy.Rda"))

###
# Extract variables
###

date <- MPdata[,1]
EXMA <- MPdata[,2]
EXME <- MPdata[,3]
EXMAUS <- MPdata[,4]
EXMEUS <- MPdata[,5]
FCI <- MPdata[,6]
Gold <- MPdata[,7]
RIMF <- MPdata[,8]
CHFsec <- MPdata[,9]
MA <- MPdata[,10]
OA <- MPdata[,11]
Assets <- MPdata[,12]
CPI <- MPdata[,13]
BCI <- MPdata[,14]
CMR <- MPdata[,15]
EuroGDP <- MPdata[,16]
SwissGDP <- MPdata[,17]

##
# Transform variables
##

totOA <- Assets-MA-RIMF-CHFsec
ldtotOAr <- diff(log(totOA/Assets))

# log returns
ldEXME <- diff(log(1/EXME))
# log differences of GDP with currency removed for Swiss GDP
ldGDPeu <- diff(log(EuroGDP))
ldGDPch <- diff(log(SwissGDP/EXMA))
# log-differences of fraction of balance sheet
ldFCIr <- diff(log(FCI/Assets))
ldRIMFr <- diff(log(RIMF/Assets))
ldCHFr <- diff(log(CHFsec/Assets))
ldGoldr <- diff(log(Gold/Assets))
ldOAr <- diff(log(OA/Assets))
ldMAr <- diff(log(MA/Assets))
# compute inflation from CPI
dCPI <- diff(CPI)/(CPI[-1])
# difference of call money rate (log not possible due to negative rates)
dCMR <- diff(CMR)

##
# Collect relevant variables
##

MPdata2 <- cbind(ldEXME,dCMR,ldFCIr,ldRIMFr,ldMAr,ldCHFr,ldtotOAr,ldGDPch,ldGDPeu,dCPI)
MPdata2 <- scale(MPdata2)
dat_all <- as.data.frame(MPdata2)

# ---- compare the full set and a subset of covariates ----
library(seqICP)

X <- dat_all[,-1]
Y <- dat_all$ldEXME
ps <- 2:4

# original set of covariates (test for CCP)
pvals1 <- rep(NA, length(ps))
for (ii in seq_along(ps)) {
  p <- ps[ii]
  message("p = ", p)
  set.seed(1)
  res1 <- seqICP(as.matrix(X), Y, test = "variance",
                 par.test = list(grid = c(0,70,140,216),
                                 complements = TRUE,
                                 link = sum,
                                 alpha = 0.05,
                                 B = 999,
                                 permutation = FALSE),
                 par.model = list(pknown = TRUE, p = p),
                 stopIfEmpty = FALSE, silent = TRUE,
                 model = "ar")
  pvals1[ii] <- max(res1$test.results$p.value)
}

# without X2, X7, or X2 and X7 (test for CCP)
rm_vars <- c(2) # c(7), c(2,7)
pvals2 <- rep(NA, length(ps))
for (ii in seq_along(ps)) {
  p <- ps[ii]
  message("p = ", p)
  set.seed(1)
  res2 <- seqICP(as.matrix(X)[,-rm_vars], Y, test = "variance",
                 par.test = list(grid = c(0,70,140,216),
                                 complements = TRUE,
                                 link = sum,
                                 alpha = 0.05,
                                 B = 999,
                                 permutation = FALSE),
                 par.model = list(pknown = TRUE, p = p),
                 stopIfEmpty = FALSE, silent = TRUE,
                 model = "ar")
  pvals2[ii] <- max(res2$test.results$p.value)
}

# original set of covariates (test for RCP)
set.seed(1)
res1 <- seqICP.s(as.matrix(X), Y, S = 1:ncol(X), test = "variance",
                 par.test = list(grid = c(0,70,140,216),
                                 complements = TRUE,
                                 link = sum,
                                 alpha = 0.05,
                                 B = 999,
                                 permutation = FALSE),
                 par.model = list(pknown = TRUE, p = 4),
                 model = "ar")
res1$p.value

# without X2, X7, or X2 and X7 (test for RCP)
set.seed(1)
rm_vars <- c(2) # c(7), c(2,7)
res2 <- seqICP.s(as.matrix(X)[,-rm_vars], Y, S = 1:ncol(X[,-rm_vars]), test = "variance",
               par.test = list(grid = c(0,70,140,216),
                               complements = TRUE,
                               link = sum,
                               alpha = 0.05,
                               B = 999,
                               permutation = FALSE),
               par.model = list(pknown = TRUE, p = 4),
               model = "ar")
res2$p.value
