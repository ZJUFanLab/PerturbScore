#' Calculate perturbscore
#' @description
#' Use ridge regression and gene network parameters to calculate the actual perturbation effect in single-cell CRISPR screening data.
#'
#'
#' @param regscore_list A regression result list calculated by regression_score function.
#' @param net_para Gene network parameters calculated by fc_extraction
#'
#' @return A result for perturbscore
#' @export
#'
#' @examples
#' pert_score <- perturbscore(regscore_list = reg_score, net_para = paras)
perturbscore <- function(regscore_list, net_para){
  if(length(regscore_list) != 2){
    stop("Please choose proper regscore_list from regression_score")
  }

  regscore<-regscore_list[[1]]
  regscore_p<-regscore_list[[2]]

  regobj<-regscore[,order(colnames(regscore))]
  netobj<-as.matrix(net_para)
  netobj<-netobj[order(rownames(netobj)),]

  regnam<-colnames(regobj)
  netnam<-names(netobj)

  con<-intersect(regnam, netnam)
  regobj<-regobj[,con]
  netobj<-netobj[con]

  netobj[netobj==0]<-1

  pertscore<-regobj
  for(i in 1:dim(regobj)[2]){
    pertscore[,i]<-regobj[,i]*netobj[i]
  }

  regscore_p<-regscore_p[,which(colnames(regscore_p) %in% colnames(pertscore))]

  return(list(pertscore, regscore_p))
}
