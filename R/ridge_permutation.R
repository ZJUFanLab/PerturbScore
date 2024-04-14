#' Gene regression
#' @description
#' The effect of perturbation is calculated by gene ridge regression.
#'
#' @param X A cell identity matrix calculated by gene_regression function.
#' @param Y A gene count matrix calculated by gene_regression function.
#' @param lamb The Regularization coefficient. Default: 0.01.
#' @param n_permut The number of permutations for p value. Default: 1000.
#'
#' @return A list for regression result.
#' @export
#'
#' @examples
#' perscore <- ridge_permutation(X, Y, lamb = 0.01, n_permut = 1000)
ridge_permutation <- function(X, Y, lamb = 0.01, n_permut = 1000){
  mat1 = (t(X) %*% X) + lamb * diag(ncol(X))
  mat2 = solve(mat1) %*% t(X) %*% Y
  mat3 = matrix(rep(0, ncol(mat2) * nrow(mat2)), nrow = nrow(mat2))
  rownames(mat3) = rownames(mat2)
  colnames(mat3) = colnames(mat2)
  for (nm in 1:n_permut) {
    cells = sample(rownames(Y), nrow(Y))
    Xm = X[cells, ]
    Ym = Y
    rownames(Ym) = cells
    mat_tes = (t(Xm) %*% Xm) + lamb * diag(ncol(Xm))
    mat_rad = solve(mat_tes) %*% t(Xm) %*% Ym
    mat3 = mat3 + (abs(mat_rad) > abs(mat2)) * 1
  }
  mat3 = mat3/n_permut
  return(list(mat2, mat3))
}
