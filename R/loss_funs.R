# ---- splitting functions ----

#' Get the list of indices after splitting an interval evenly into sub-intervals.
#'
#' @param st Starting index (inclusive).
#' @param en Ending index (inclusive).
#' @param n_seg Number of segments.
#'
#' @return A list of sequence (vector) of indices.
#' 
#' @export
get_splits <- function(st, en, n_seg = 3) {
  full_seq <- st:en
  len_seg <- floor((en-st+1)/n_seg)
  idx_lst <- lapply(seq_len(n_seg-1), \(ii) {
    (st + (ii-1) * len_seg):(st + ii * len_seg - 1)
  })
  idx_lst[[n_seg]] <- setdiff(full_seq, c(unlist(idx_lst)))
  return(idx_lst)
}

#' Get the list of indices of one-vs-rest after splitting evenly into sub-intervals.
#'
#' @details In the inner list the shorter interval comes first.
#'
#' @param st The start index (inclusive).
#' @param en The end index (inclusive).
#' @param n_split Number of splits.
#'
#' @return A list of lists where each inner list contains two vectors.
#' 
#' @export
get_one_vs_rest <- function(st, en, n_seg = 3) {
  idx_lst <- get_splits(st, en, n_seg)
  lst_one_vs_rest <- lapply(1:n_seg, \(ii) {
    idx_lst_rm <- idx_lst
    idx_lst_rm[[ii]] <- NULL
    list(idx_lst[[ii]], c(unlist(idx_lst_rm)))
  })
  return(lst_one_vs_rest)
}

#' Get the list of indices of rest-vs-one after splitting evenly into sub-intervals.
#'
#' @details In the inner list the longer interval comes first.
#'
#' @param st Starting index (inclusive).
#' @param en Ending index (inclusive).
#' @param n_split
#'
#' @return A list of lists where each inner list contains two vectors.
#' 
#' @export
get_rest_vs_one <- function(st, en, n_split = 3) {
  idx_lst <- get_splits(st, en, n_split)
  lst_rest_vs_one <- lapply(1:n_split, \(ii) {
    idx_lst_rm <- idx_lst
    idx_lst_rm[[ii]] <- NULL
    list(c(unlist(idx_lst_rm)), idx_lst[[ii]])
  })
  return(lst_rest_vs_one)
}

# ---- calculate loss ----

#' Fit a linear model on a subset of covariates and return the residuals.
#'
#' @details Given two sets of data, the 1st set is used as the training set
#'   and the 2nd is used as the test set. Both the training residual and the
#'   test residuals are returned.
#'
#' @param X_tr A matrix of dimension `n_tr x d`.
#' @param Y_tr A vector of length `n_tr`.
#' @param X_te A matrix of dimension `n_te x d`.
#' @param Y_te A vector of length `n_te`.
#' @param S A vector of either indices or names corresponding to the columns of
#'   `X_tr` and `X_te`.
#' @param intercept Whether to include an intercept in the linear regressions.
#'
#' @return A list of two vectors of length `n_tr` and `n_te` respectively.
#'
get_res <- function(X_tr, Y_tr, X_te, Y_te, S, intercept = TRUE) {
  if (intercept) {
    XS_tr <- cbind(1, X_tr[,S,drop=FALSE])
    XS_te <- cbind(1, X_te[,S,drop=FALSE])
  }
  mod <- lm(Y_tr ~ XS_tr - 1)
  res_in <- resid(mod)
  res_out <- c(Y_te - XS_te %*% coef(mod))
  return(list(res_in = res_in,
              res_out = res_out))
}

#' Calculate the causal change loss on one side of the split.
#'
#' @param split_lst A list of splits. Each element is a list of two vectors of
#'   indices, where the 1st set of indices can be seen as the training set
#'   and the 2nd set can be seen as the test set.
#' @param cond_sets The conditioning sets. Default to the power set of all covariates.
#' @param intercept Whether to include an intercept in the linear regressions.
#'
#' @return A scalar value.
cc_loss_one_side <- function(X, Y, split_lst,
                             cond_sets = get_powerset(colnames(X)),
                             intercept = TRUE) {
  interval <- sort(unique(unlist(split_lst)))
  # vals <- matrix(0, nrow = length(cond_sets), ncol = length(split_lst))
  vals <- rep(0, length(cond_sets))
  for (ii in seq_along(cond_sets)) {
    S <- cond_sets[[ii]]
    for (jj in seq_along(split_lst)) {
      idx_lst <- split_lst[[jj]]
      # res_lst <- get_res(X[interval,,drop=FALSE], Y[interval],
      #                    X[idx_lst,,drop=FALSE], Y[idx_lst],
      #                    S, intercept)
      res_lst <- get_res(X[idx_lst[[1]],,drop=FALSE], Y[idx_lst[[1]]],
                         X[idx_lst[[2]],,drop=FALSE], Y[idx_lst[[2]]],
                         S, intercept)
      # vals[ii, jj] <- abs(mean(res_lst$res_out))
      vals[ii] <- vals[ii] + (mean(res_lst$res_out^2) - mean(res_lst$res_in^2))^2
    }
  }
  return(vals)
}

#' Causal change loss function.
#'
#' @param X The covariate matrix of dimension `n x d`.
#' @param Y The response vector of length `n`.
#' @param i The time index at which to evaluate the loss, relative to `st`.
#' @param st The start index (inclusive).
#' @param en The end index (inclusive).
#' @param split_setting A list of two elements: `type` and `n_split`. `type` can
#'   be one of "one_vs_rest" or "rest_vs_one"; `n_split` is the number of splits
#'   over the interval `st:(st+i-1)` and `(st+i):n` for calculating the loss.
#' @param cond_sets The conditioning sets. Default to the power set of all covariates.
#' @param intercept Whether to include an intercept in the linear regressions.
#'
#' @examples
#' tmp <- sapply(eval_grid, \(t) cc_loss(X, Y, t, 1, n,
#'   split_setting = list(type = "one_vs_rest", n_split  = 3), cond_sets = PS))
#'
#' @return A scalar value.
#' 
#' @export
cc_loss <- function(X, Y, i, st, en,
                    min_seg = 20,
                    cond_sets = get_powerset(colnames(X)),
                    intercept = TRUE,
                    verbose = FALSE) {
  n_seg_left <- floor((i-1)/min_seg)
  n_seg_right <- floor((en-st-i+2)/min_seg)
  if (n_seg_left < 2) {
    vals_left <- 0
  } else {
    split_lst_left <- get_one_vs_rest(st, st+i-1, n_seg_left)
    vals_left <- cc_loss_one_side(X, Y, split_lst_left,
                                  cond_sets = cond_sets,
                                  intercept = intercept)
  }
  if (n_seg_right < 2) {
    vals_right <- 0
  } else {
    split_lst_right <- get_one_vs_rest(st+i, en, n_seg_right)
    vals_right <- cc_loss_one_side(X, Y, split_lst_right,
                                   cond_sets = cond_sets,
                                   intercept = intercept)
  }
  return((min(vals_left) + min(vals_right))/(n_seg_left + n_seg_right))
}

get_segs_one_vs_rest_all <- function(search_grid, n_tot, rev_one_rest = FALSE) {
  segs_individuals <- get_segs(search_grid, n_tot)
  segs <- lapply(segs_individuals, function(m) {
    if(rev_one_rest) {
      list(setdiff(1:n_tot, m), m)
    } else {
      list(m, setdiff(1:n_tot, m))
    }
  })
  return(segs)
}

#' Compare multiple environment with one vs rest.
#'
#' @details Segment the observations according to relative grid point locations and
#'   compare one vs the rest.
#'   
#' @param X_seg A (segment of the) covariate matrix.
#' @param Y_seg A (segment of the) response vector`.
#' @param PS A list of conditioning sets e.g. returned by `get_powerset`.
#' @param rel_grid A vector of scalars in (0, 1) representing the relative grid, e.g., c(1/3, 2/3)
#' @param intercept A boolean whether to include an intercept in the regression.
#' @param test_name A string "chow1" or "chow2".
#' @param two_dir A boolean whether to test both one vs rest and rest vs one.
#' @param save_resids A boolean whether to return residuals.
#' 
#' @return A list containing p-value and residuals (which could be NULL if 
#'   `save_resids` is FALSE).
#' 
#' @export 
one_vs_rest_test <- function(X_seg, Y_seg, PS, rel_grid, test_name = "chow1",
                             intercept = TRUE, rev_one_rest = FALSE,
                             two_dir = FALSE, save_resids = FALSE) {
  n_obs <- nrow(X_seg)
  tmp_grid <- ceiling(n_obs * rel_grid + 1)
  segs <- get_segs_one_vs_rest_all(tmp_grid, n_obs, rev_one_rest)
  pvals_all <- vector("list", length(segs))
  for (ii in seq_along(segs)) {
    seg <- segs[[ii]]
    env1 <- seg[[1]]
    env2 <- seg[[2]]
    if (two_dir) {
      pvals_all[[ii]] <- vector("list", 2)
      res_lst1 <- two_env_test(X_seg, Y_seg, env1, env2, PS, test_name, intercept, save_resids)
      res_lst2 <- two_env_test(X_seg, Y_seg, env2, env1, PS, test_name, intercept, save_resids)
      pvals_all[[ii]][[1]] <- res_lst1$pvals
      pvals_all[[ii]][[2]] <- res_lst2$pvals
      # when save_resids == FALSE, this will just extract NULL
      if (!rev_one_rest & ii == 1) {
        # 1 and 2 have the same mod0 and only should save from the 1st comparison
        #. if not reversed, the rest are from reordered X and Y segments
        resids_all <- res_lst1$resids
      }
      else if (rev_one_rest & ii == (length(rel_grid) + 1)) {
        # 1 and 2 have the same mod0 and only should save from the last comparison
        #. if reversed, the rest are from reordered X and Y segments
        resids_all <- res_lst1$resids
      }
    }
    else {
      res_lst <- two_env_test(X_seg, Y_seg, env1, env2, PS, test_name, intercept, save_resids)
      pvals_all[[ii]] <- res_lst$pvals
      if (!rev_one_rest & ii == 1) {
        resids_all <- res_lst$resids
      }
      else if (rev_one_rest & ii == (length(rel_grid) + 1)) {
        resids_all <- res_lst$resids
      }
    }
  }
  return(list(pvals_all = pvals_all,
              resids_all = resids_all))
}

two_env_test <- function(X_seg, Y_seg, env1, env2, ps, test_name = "chow1",
                         intercept = TRUE, save_resids = FALSE) {
  test_fun <- eval(parse(text = paste0("test_", test_name)))
  # ps <- get_powerset(colnames(X_seg))
  # if (!intercept) {
  #   ps[[1]] <- NULL
  # }
  pvals <- rep(NA, length(ps))
  resids <- vector("list", length(ps))
  for (ii in seq_along(ps)) {
    S <- ps[[ii]]
    res_lst <- test_fun(X_seg, Y_seg, env1, env2, S, intercept = intercept)
    pvals[ii] <- res_lst$p_value
    resids[[ii]] <- res_lst$resids
  }
  if (save_resids) {
    return(list(pvals = pvals,
                resids = resids))
  }
  else {
    return(list(pvals = pvals))
  }
}
