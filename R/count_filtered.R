#' Filter genes
#' @description
#' Filter low expression genes and remove cells with high expression of mitochondiral genes.
#'
#'
#' @param count A gene count matrix. Before input the count data,it is recommended to finish the imputation of count matrix.
#' @param frac A parameter for filtering low expression genes. By default, only genes that have expression in at least that fractions of cells are kept. Default: 0.01.
#' @param nfeatures A parameter for extracting high-variable genes for gene co-expression network construction. By default, the logarithmic curve between mean and variance is fitted based on the loess. Then order the standardized variances and extract nfeatures genes. Default: 2000. NULL: extract all genes.
#' @param mt.percent A parameter for filtering cells with high expression of mitochondrial genes. Default: 5. Filtering out cells that express more than 5% of mitochondrial genes. NULL: No cells are filtered out.
#'
#' @return A gene count matrix
#' @export
#'
#' @importFrom stats loess sd
#' @examples
#' count count <- PerturbScore::count_data
#' count_filter <- count_filtered(count, frac = 0.01, nfeatures = 2000, mt.percent = 5)
#'
count_filtered <- function(count, frac = 0.01, nfeatures = 2000, mt.percent = 5){

  obj <- count

  if(is.null(frac)){
    stop("Please select proper value for frac.")
  }

  if(!is.null(mt.percent)){
    mt_genes <- grep(pattern = "^MT-", x = rownames(obj), value = T)
    percent.featureset <- colSums(obj[mt_genes,])/colSums(obj) * 100
  }
  else {
    percent.featureset <- 0
  }

  sel_gene = rownames(obj)[which(rowSums(obj != 0) >= ncol(obj) * frac)]
  if(length(sel_gene) == 0){
    stop("No genes left after selected. Please check options again.")
  }

  con_sel <- obj[sel_gene,]
  sd_row<-c()
  if(is.null(nfeatures)){
    con_sel <- con_sel[,percent.featureset<5]
    return(con_sel)
  }
  else {
    exp_filtered <- con_sel[, percent.featureset<5]

    genemean <- as.data.frame(rowMeans(exp_filtered))
    genesd <- as.data.frame(apply(exp_filtered, 1, sd))
    genedata <- cbind(genemean, genesd^2)
    genedata <- cbind(rownames(exp_filtered), genedata)
    colnames(genedata) <- c("gene", "mean", "variance")

    fit <- loess(formula = log10(x = variance) ~ log10(x = mean), data = genedata, span = 0.3)

    genedata$variance_adjust <- 10 ^ fit$fitted
    genedata$variance_standard <- genedata$variance / genedata$variance_adjust
    genedata <- genedata[order(genedata$variance_standard,decreasing = T),]
    genes <- genedata$gene[1:nfeatures]
    exp_filtered <- exp_filtered[which(rownames(exp_filtered) %in% genes),]
    return(exp_filtered)
  }
}
