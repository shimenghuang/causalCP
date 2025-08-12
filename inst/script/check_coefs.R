library(dplyr)
library(ggplot2)
library(seqICP)

# ---- helper functions ----

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
  X_new <- X_new[(max(shifts) + 1):nrow(X), ]
  col_names <- colnames(X_new)
  X_new <- as.matrix(X_new)
  colnames(X_new) <- col_names
  Y_new <- Y[(max(shifts) + 1):length(Y)]
  return(list(X = X_new, Y = Y_new))
}

drop_vars <- function(X, names) {
  keep_names <- setdiff(colnames(X), names)
  return(X[, keep_names])
}

# ---- load data ----

pkg_dir <- "."
load(paste0(pkg_dir, "/inst/data/monetarypolicy.Rda"))

###
# Extract variables
###

date <- MPdata[, 1]
EXMA <- MPdata[, 2]
EXME <- MPdata[, 3]
EXMAUS <- MPdata[, 4]
EXMEUS <- MPdata[, 5]
FCI <- MPdata[, 6]
Gold <- MPdata[, 7]
RIMF <- MPdata[, 8]
CHFsec <- MPdata[, 9]
MA <- MPdata[, 10]
OA <- MPdata[, 11]
Assets <- MPdata[, 12]
CPI <- MPdata[, 13]
BCI <- MPdata[, 14]
CMR <- MPdata[, 15]
EuroGDP <- MPdata[, 16]
SwissGDP <- MPdata[, 17]

##
# Transform variables
##

totOA <- Assets - MA - RIMF - CHFsec
ldtotOAr <- diff(log(totOA / Assets))

# log returns
ldEXME <- diff(log(1 / EXME))
# log differences of GDP with currency removed for Swiss GDP
ldGDPeu <- diff(log(EuroGDP))
ldGDPch <- diff(log(SwissGDP / EXMA))
# log-differences of fraction of balance sheet
ldFCIr <- diff(log(FCI / Assets))
ldRIMFr <- diff(log(RIMF / Assets))
ldCHFr <- diff(log(CHFsec / Assets))
ldGoldr <- diff(log(Gold / Assets))
ldOAr <- diff(log(OA / Assets))
ldMAr <- diff(log(MA / Assets))
# compute inflation from CPI
dCPI <- diff(CPI) / (CPI[-1])
# difference of call money rate (log not possible due to negative rates)
dCMR <- diff(CMR)

##
# Collect relevant variables
##

MPdata2 <- cbind(ldEXME, dCMR, ldFCIr, ldRIMFr, ldMAr, ldCHFr, ldtotOAr, ldGDPch, ldGDPeu, dCPI)
MPdata2 <- scale(MPdata2)
dat_all <- as.data.frame(MPdata2)

X <- dat_all[, -1]
Y <- dat_all$ldEXME

p <- 2
set.seed(1)
res1 <- seqICP(as.matrix(X), Y,
  test = "variance",
  par.test = list(
    grid = c(0, 70, 140, 216),
    complements = TRUE,
    link = sum,
    alpha = 0.05,
    B = 999,
    permutation = FALSE
  ),
  par.model = list(pknown = TRUE, p = p),
  stopIfEmpty = FALSE, silent = TRUE,
  model = "ar"
)
(S <- res1$parent.set) # 2 7

coef_tab <- data.frame(res1$coefficients)
colnames(coef_tab) <- c("coef", "lwr", "upr")
coef_ci <- tibble::tibble(coef_tab) %>%
  dplyr::mutate(coef_ci = paste0(round(coef, 3), " (", round(lwr, 3), ", ", round(upr, 3), ")")) %>%
  dplyr::pull(coef_ci)
coef_names <- c(
  "Intercept", paste0("X", 1:ncol(X)),
  paste0(paste0(c("Y", paste0("X", 1:ncol(X))), "_lag"), rep(1:p, each = 1 + ncol(X)))
)
coef_ci <- cbind(coef_names, coef_ci)
write.csv(coef_ci, file = "inst/data/coef_ci.csv", row.names = FALSE)

tab <- as.data.frame(coef_ci) |>
  tidyr::separate(coef_names, c("coef", "lag"), "_") |>
  tidyr::pivot_wider(names_from = lag, values_from = coef_ci) |>
  dplyr::rename("Parameter" = coef, "lag0" = `NA`)

knitr::kable(tab, booktabs = TRUE, format = "latex", align = "lrrrrr")
