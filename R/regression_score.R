#' Calculate perturbation effect
#' @description
#' Calculate perturbation effect by ridge regression in single-cell CRISPR screening data.
#'
#'
#' @param count A gene count matrix.
#' @param sgrna A perturbation-cell data frame. Noted that the data frame contain perturbation and cell parameters.
#' @param ctrl A parameter for non-perturbed cells which is recommended as NT.
#' @param nfeatures High-variable gene number. Default: 2000.
#' @param n_permut The number of permutations for p value. Default: NULL.
#' @param paramt The Regularization coefficient. Default: 0.01.
#' @param filt A parameter for filtering low expression genes. Default: 0.01.
#' @param mt.percent A parameter for filtering cells with high expression of mitochondrial genes. Default: 5.
#' @param slot A parameter for normalizing data. Default: 'scale'.
#'
#' @return A list for regression result.
#' @export
#'
#' @examples
#' count <- PerturbScore::count_data
#' sgrna <- PerturbScore::perturbation_data
#' reg_score <- regression_score(count,sgrna,ctr = "NT", nfeatures = 2000,n_permut = 1000,paramt = 0.01,filt = 0.01,mt.percent = 5, slot = 'scale')
regression_score <- function(count, sgrna, ctrl, nfeatures = 2000, n_permut = NULL, paramt = 0.01, filt = 0.01, mt.percent = 5, slot = 'scale'){
  if (!is.null(n_permut)) {
    permutation <- as.integer(n_permut)
  } else {
    permutation <- 1000
  }

  sgr_count <- sgrna
  if (sum(colnames(sgr_count) %in% c("cell", "perturbation")) != 2) {
    stop("cell, perturbation column names not found in sgrna file.")
  }

  tar_count <- count

  mat <- sum(sgr_count[, which(colnames(sgr_count)=="cell")] %in% colnames(tar_count))
  if (mat == 0) {
    stop("cell names in expression matrix and sgrna matrix do not match")
  }

  ind <- cell_info(tar_count, sgr_count)#

  mat_reg = gene_regression(count = tar_count, cell_ident = ind, selected = nfeatures, ctrl = ctrl, gene_frac = filt, slot = slot,mt.percent = mt.percent)

  Xmat = mat_reg[[1]]
  Ymat = mat_reg[[2]]

  perscore = ridge_permutation(Xmat, Ymat, lamb = paramt, n_permut = permutation)

  return(perscore)

}
