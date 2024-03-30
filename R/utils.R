#' Generate all subsets of a set.
#'
#' @details This function is used to generate all potential invariant sets.
#'
#' @param S A vector with no duplicated values representing a set.
#'
#' @return A list of all subsets of `S`.
#' 
#' @export
get_powerset <- function(S) {
  n <- length(S)
  masks <- 2^(1:n-1)
  sets <- lapply(1:2^n-1, function(u) S[bitwAnd(u, masks) != 0])
  return(sets)
}

#' Construct segments of indicies based on split points.
#' @param cps A vector of integers, the split points (the points themselves are
#'   considered as the starting index of the next segment).
#' @param n_tot An integer, total number of indicies.
#'
#' @return A list of length \code{length(cps) + 1}, each containing the indices
#'   in the segment.
#'   
#' @export
get_segs <- function(cps, n_tot) {
  cps_su <- sort(unique(cps))
  seg_start <- c(1, cps_su)
  seg_end <- c(cps_su-1, n_tot)
  segs <- mapply(\(l,r) {l:r}, seg_start, seg_end, SIMPLIFY = FALSE)
  return(segs)
}

#' Calculate the length of each segment.
#'
#' @param cps_su An integer vector of sorted and unique change points.
#' @param n_obs An integer, number of observations.
#'
#' @return An integer vector of \code{length(cps_su) + 1}.
#' 
#' @export
get_seg_lens <- function(cps_su, n_obs) {
  seg_lens <- unname(diff(c(1, cps_su, n_obs)))
  seg_lens[length(seg_lens)] <- seg_lens[length(seg_lens)] + 1
  return(seg_lens)
}


#' Summarize the types of the change points.
#'
#' @details Given a vector of time points of causal and non-causal changes,
#'   sort them and returns the type(s) of each change point. There can be
#'   overlaps between the two sets, which are considered as CCPs by our
#'   definition.
#'
#' @param loc_ccp A integer vector of time points where causal mechanism changes.
#' @param loc_nccp A integer vector of time points where non-causal mechanism
#'   changes.
#'
#' @return A list containing the type(s) of each change point. The types are
#'   "CCP" or "NCCP". The change points are sorted from earliest to latest.
#'   
#' @export
get_cp_types <- function(loc_ccp, loc_nccp) {
  # get the type(s) of change points (Note:an idx could be both a CCP and NCCP)
  cps <- c(loc_ccp, loc_nccp)
  names(cps) <- c(rep("CCP", length(loc_ccp)), rep("NCCP", length(loc_nccp)))
  cps_su <- sort(unique(cps)) # sorted and unique change points
  cp_types <- vector("list", length = length(cps_su))
  for (ii in seq_along(cps_su)) {
    cp_types[[ii]] <- names(cps[cps == cps_su[ii]])
  }
  names(cp_types) <- cps_su
  return(cp_types)
}

