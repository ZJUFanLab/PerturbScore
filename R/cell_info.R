#' Construct perturbation-cell matrix
#'
#' @description
#' Construct perturbation-cell matrix. The cell is assigned with TRUE which is infected by corresponding perturbation
#'
#' @param count A gene-cell matrix.
#' @param sgrna A perturbation data frame which contains perturbation and cell parameters.
#'
#' @return A perturbation-cell matrix
#' @export
#'
#' @examples
#' count <- PerturbScore::count_data
#' sgrna <- PerturbScore::perturbation_data
#' ind <- cell_info(count = count, sgrna = count)
cell_info <- function(count, sgrna){
  rnm = unique(sgrna[,which(colnames(sgrna)=="cell")])
  cnm = unique(sgrna[,which(colnames(sgrna)=="perturbation")])

  sc_count<-count

  rnm = rnm[!is.na(rnm)]
  rnm = rnm[rnm %in% colnames(sc_count)]
  if (length(rnm) == 0) {
    stop("Cell names do not match in expression matrix and barcode.")
  }

  cnm = cnm[!is.na(cnm)]

  ind = matrix(rep(FALSE, length(rnm) * length(cnm)), nrow = length(rnm))
  rownames(ind) = rnm
  colnames(ind) = cnm
  row <- sgrna[, 'cell']
  col <- sgrna[, 'perturbation']
  test <- (row %in% rnm) & (col %in% cnm)
  idx <- cbind(row[test], col[test])
  ind[idx]  <- TRUE
  return(ind)
}
