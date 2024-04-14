#' Dot plot
#' @description
#' Dot plot for perturbation effects
#'
#'
#' @param pert_score Result from perturbscore function.
#' @param num A parameter for extracting top phenotype genes by pvalue. Default: NULL.
#' @param features Interested phenotype gene list. Default: NULL.
#' @param perturbations Interested perturbation list. Default: NULL.
#'
#' @return A dot plot
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr gather
#' @examples
#' p <- pert_feature_plot(pert_score,num=2,features=NULL,perturbations=NULL)
#' p
pert_feature_plot<-function(pert_score,num=NULL,features=NULL,perturbations=NULL){
  res_score<-data.frame(pert_score[[1]])
  res_score<-res_score[,order(colnames(res_score))]
  res_p<-data.frame(pert_score[[2]])
  res_p<-res_p[,order(colnames(res_p))]
  res_score<-cbind(rownames(res_score),res_score)
  res_p<-cbind(rownames(res_p),res_p)
  colnames(res_score)[1]<-"perturbation"
  colnames(res_p)[1]<-"perturbation"

  res_score<-tidyr::gather(res_score,key = "gene",value = "perturbscore",-"perturbation")
  res_p<-tidyr::gather(res_p,key = "gene",value = "pvalue",-"perturbation")
  res_score<-cbind(res_score,res_p$pvalue)
  colnames(res_score)[4]<-"pvalue"
  res_score$pvalue[which(res_score$pvalue==0)]<-0.001
  res_score$logpvalue<- -log10(res_score$pvalue)

  p_data<-data.frame(perturbation = res_score$perturbation,
                     gene = res_score$gene,
                     perturbscore=res_score$perturbscore,
                     pvalue = res_score$pvalue,
                     logpvalue = res_score$logpvalue)

  if((is.null(features))&(is.null(perturbations))&(is.null(num))){
    stop("Please choose proper num, features and perturbations")
  }
  if(is.null(num)==FALSE){
    p_data<-p_data %>% dplyr::group_by(perturbation) %>% dplyr::slice_min(pvalue,n=num,with_ties=FALSE)
  }
  if(is.null(features)==FALSE){
    p_data<-p_data[which(p_data$gene %in% features),]
  }
  if(is.null(perturbations)==FALSE){
    p_data<-p_data[which(p_data$perturbation %in% perturbations),]
  }

  ggplot(p_data ,aes(gene,perturbation,size =logpvalue ))+
    geom_point(shape=21,aes(fill= perturbscore),position =position_dodge(0))+
    theme_minimal()+
    xlab(NULL)+ylab(NULL) +
    scale_size_continuous(range=c(1,6),name = "-log(Pvalue)",
                          limits=c(0,4), breaks = c(1,2,3))+
    scale_fill_gradient(high= "#E54924", low = "#498EA4",name="PerturbScore",
                        limit=c(-3,5), breaks=c(-3,1,5))+
    theme_bw()+
    theme(legend.position = "right",legend.box = "vertical", #图例位置
          legend.margin=margin(t= 1, unit='cm'),
          axis.text.x  = element_text(color="black",size=12,angle = 45, vjust = 1,hjust = 1),#x轴
          axis.text.y  = element_text(color="black",size=11),#y轴
          axis.title.y=element_text(vjust=1,size=14),
          panel.grid = element_blank(),
          legend.key = element_blank()
    )+
    labs(y="Perturbations",x = "Phenotypes",title = paste(paste("Top",num,sep = ""),"perturbed phenotype genes"," "))
}
