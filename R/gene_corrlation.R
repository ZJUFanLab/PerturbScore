#' Compute genes correlation
#' @description
#' Calculate correlation coefficient between genes.
#'
#' @param count A gene count matrix.
#'
#' @return A gene correlation matrix
#' @export
#' @importFrom stats cor
#' @examples
#' gene_corr <- gene_corrlation(count = count_filter)
gene_corrlation <- function(count){
  corr <- cor(t(count), method = "spearman")
  return(corr)
}
