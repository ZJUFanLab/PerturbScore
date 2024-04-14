#' Extract gene correlation coefficient
#' @description
#' Extract gene network parameters included normalized betweenness centrality and gene connectivity.
#'
#' @param corr A correlation matrix calculated by gene_corrlation function.
#' @param do.normalization A parameter whether to  normalize betweenness centrality. Default: T.
#' @param abs A parameter Whether to take the absolute value of the result. Default: T
#'
#' @return A gene network parameters data frame.
#' @export
#' @importFrom igraph graph_from_adjacency_matrix  degree
#' @examples
#' paras <- paras_extraction(corr = gene_corr, do.normalization = T, abs = T)
paras_extraction <- function(corr, do.normalization = T, abs = T){
  test <-abs(corr)
  test[test<0.6]<-0

  g <- igraph::graph_from_adjacency_matrix(test,mode = "undirected",weighted = TRUE,diag = FALSE)

  if(do.normalization) {
    sca_deg<-igraph::degree(g,normalized = T)
  }
  else {
    sca_deg<-igraph::degree(g,normalized = F)
  }
  names(sca_deg)<-rownames(test)

  corsum<-colSums(test)-1
  names(corsum)<-rownames(test)

  fc=sca_deg*corsum
  fc<-data.frame(fc)

  if(abs){
    fc<-abs(fc)
    return(fc)
  }
  else {
    return(fc)
  }
}
