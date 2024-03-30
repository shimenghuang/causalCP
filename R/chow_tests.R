#' Chow test (eq 29 in Chow 1960)
#'
#' @return A list of the test statistic and the p-value.
#' 
#' @export
test_chow1 <- function(X, Y, I1, I2, S, intercept = TRUE) {
  if (!intercept & length(S) == 0) stop("Cannot fit empty set without intercept.")
  p <- ifelse(intercept, length(S) + 1, length(S))
  df1 <- p
  df2 <- length(c(I1, I2)) - 2*p

  if (intercept) {
    X0 <- cbind(1, X[c(I1, I2),S,drop=FALSE])
    X1 <- cbind(1, X[I1,S,drop=FALSE])
    X2 <- cbind(1, X[I2,S,drop=FALSE])
  } else {
    X0 <- X[c(I1, I2),S,drop=FALSE]
    X1 <- X[I1,S,drop=FALSE]
    X2 <- X[I2,S,drop=FALSE]
  }
  y0 <- Y[c(I1, I2)]
  y1 <- Y[I1]
  y2 <- Y[I2]
  mod0 <- lm(y0 ~ X0 - 1)
  mod1 <- lm(y1 ~ X1 - 1)
  mod2 <- lm(y2 ~ X2 - 1)

  num <- sum( (predict(mod1) - X1 %*% coef(mod0))^2 ) +
    sum( (predict(mod2) - X2 %*% coef(mod0))^2 )
  den <- sum(resid(mod1)^2) + sum(resid(mod2)^2)
  stat <- num / den * df2 / df1
  pval <- pf(stat, df1 = df1, df2 = df2, lower.tail = F)

  return(list(test_stat = stat,
              p_value = pval,
              resids = resid(mod0)))
}

#' Chow test (eq 30 in Chow 1960)
#'
#' @return A list of the test statistic and the p-value.
#' 
#' @export
test_chow2 <- function(X, Y, I1, I2, S, intercept = TRUE) {
  if (!intercept & length(S) == 0) stop("Cannot fit empty set without intercept.")
  p <- ifelse(intercept, length(S) + 1, length(S))
  df1 <- length(I2)
  df2 <- length(I1) - p

  if (intercept) {
    X0 <- cbind(1, X[c(I1, I2),S,drop=FALSE])
    X1 <- cbind(1, X[I1,S,drop=FALSE])
    X2 <- cbind(1, X[I2,S,drop=FALSE])
  }
  else {
    X0 <- X[c(I1, I2),S,drop=FALSE]
    X1 <- X[I1,S,drop=FALSE]
    X2 <- X[I2,S,drop=FALSE]
  }
  y0 <- Y[c(I1, I2)]
  y1 <- Y[I1]
  y2 <- Y[I2]
  mod0 <- lm(y0 ~ X0 - 1)
  mod1 <- lm(y1 ~ X1 - 1)

  num <- sum( (predict(mod1) - c(X1 %*% coef(mod0)))^2 ) +
    sum( (y2 - c(X2 %*% coef(mod0)))^2 )
  den <- sum( (y1 - c(X1 %*% coef(mod1)))^2 )
  stat <- num / den * df2 / df1
  pval <- pf(stat, df1 = df1, df2 = df2, lower.tail = F)

  return(list(test_stat = stat,
              p_value = pval,
              resids = resid(mod0)))
}
