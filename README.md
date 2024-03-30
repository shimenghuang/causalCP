# Causal change point detection and localization

This repository contains code and data for reproducing the results in [Huang, S., Peters, J., & Pfister, N. (2024). Causal Change Point Detection and Localization. arXiv preprint arXiv:2403.12677](https://arxiv.org/abs/2403.12677).

## Installation

The developmenet R package `causalCP` can be installed by 
```r
# install.packages("remotes")
remotes::install_github("shimenghuang/causalCP")
```

## A simple usage example

```r
dgp <- function(n1 = 500, n2 = 500,
                alpha1 = 1, alpha2 = 1,
                beta1 = 1, beta2 = 2,
                sy1 = 1, sy2 = 1) {
  X1 <- rnorm(n1, 1)
  X2 <- rnorm(n2, 2)
  X <- c(X1, X2)
  y1 <- c(alpha1 + X1 * beta1 + rnorm(n1, 0, sy1))
  y2 <- c(alpha2 + X2 * beta2 + rnorm(n2, 0, sy2))
  y <- c(y1, y2)
  return(list(X = matrix(X, ncol = 1), y = y))
}

dat <- dgp()
cond_sets <- get_powerset(1)

# CCP detection
max(sapply(cond_sets, \(s) {
  test_chow1(dat$X, dat$y, 1:500, 501:1000, s)$p_value
}))

# CCP localization
eval_grid <- seq(10, 990, by = 5)
vals <- sapply(eval_grid, \(t) cc_loss(dat$X, dat$y, t, 1, 1000,
                                       min_seg = 30,
                                       cond_sets = cond_sets))
eval_grid[which.min(vals)]
```


