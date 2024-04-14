#' Construct data for gene ridge regression and permutation
#'
#' @param count A gene count matrix.
#' @param cell_ident A perturbation-cell matrix calculated by cell_info function.
#' @param selected High-variable gene number. Default: 2000.
#' @param ctrl A parameter for control cell. Non-perturbed cells are recommended to assign NT.
#' @param gene_frac A parameter for filtering low expression genes. Default: 0.01.
#' @param mt.percent A parameter for filtering cells with high expression of mitochondrial genes. Default: 5.
#' @param slot A parameter for normalizing data. Default: 'scale'.
#'
#' @return A list for gene regression. Xmat represents cell identity matrix that contains cell and perturbation parameters. Ymat represents gene count matrix.
#' @export
#' @importFrom stats quantile
#' @examples
#' mat_reg <- gene_regression(count, ind, selected = 2000,ctrl="NT",
#'                             gene_frac = 0.01, slot = 'scale')
gene_regression<-function(count, cell_ident, selected = NULL,ctrl, gene_frac = 0.01, mt.percent=5, slot = 'scale'){
  threshd = 0.95
  raw_count = count
  if (slot=="count"){
    sc_count<-count
  }
  if (slot=="scale"){
    sc_count <-apply(count, 2, function(x){ log(x/sum(x)*1e4 +1) })
    sc_count <- scale(sc_count, center = T, scale = T)
  }

  sel_genes = rownames(sc_count)
  sel_genes = sel_genes [!is.na(sel_genes)]

  if (is.null(selected) == FALSE) {
    count_filter <- count_filtered(raw_count, frac = gene_frac, nfeatures = selected,  mt.percent = mt.percent)
    sel_genes <- sel_genes[sel_genes %in% rownames(count_filter)]
    if (length(sel_genes) == 0) {
      stop("No genes left after selected")
    }
  }

  sel_genes = sel_genes[sel_genes %in% rownames(raw_count) ]
  raw_count2 = as.matrix(raw_count[sel_genes,])
  sel_genes = rownames(raw_count2)[which(rowSums(raw_count2 != 0) >= ncol(raw_count2) * gene_frac)]

  if (length(sel_genes) == 0) {
    stop("No genes left after filtered")
  }

  sel_cells = rownames(cell_ident)
  sel_cells = sel_cells[sel_cells %in% colnames(sc_count)]

  sel_cells = sel_cells[!is.na(sel_cells)]

  Ymat = as.matrix(sc_count[sel_genes, sel_cells])
  Ymat = t(Ymat)

  tgf = colnames(cell_ident)
  tgf[tgf %in% ctrl] = "NT"
  tgp = as.factor(tgf)

  Xmat = matrix(rep(0, length(sel_cells) * length(unique(tgp))), nrow = length(sel_cells))
  rownames(Xmat) = sel_cells
  colnames(Xmat) = levels(tgp)
  for (cnl in colnames(cell_ident)) {
    cellns = which(cell_ident[, cnl] == TRUE)
    if (cnl %in% ctrl) {
      Xmat[cellns, "NT"] = 1
    } else {
      Xmat[cellns, cnl] = 1
    }
  }
  Xmat[, "NT"] = 1

  Ymat_out = apply(Ymat, 2, function(X) {
    return(quantile(X, probs = threshd))
  })
  out_mat = t(matrix(rep(Ymat_out, nrow(Ymat)), ncol = nrow(Ymat)))
  Ymat_cor = ifelse(Ymat > out_mat, out_mat, Ymat)
  Ymat = Ymat_cor

  return(list(Xmat, Ymat))
}
