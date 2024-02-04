
#' Function to generate weights for signatures if multiple signatures found
#' @param exprs_df a data frame of expression data
#' @return vector of weights
#' @export
weight_signature <- function(exprs_df) {
    min_cor <- 0.05
    if (ncol(exprs_df) == 1) {
        return(c(1))
    } else if (ncol(exprs_df) == 2) {
        return(c(0.5,0.5))
    }

    cor_df <- cor(exprs_df, method="spearman")
    mean_cor <- (rowSums(cor_df) - 1) / (length(cor_df) - 1)
    weights <- sapply(mean_cor,function(x) max(x,min_cor))
    weights <- weights / sum(weights)
    return(weights)
}

#' Function to weight signatures if multiple signatures found
#' @param z_scores a data frame of expression data
#' @param weights a vector of weights
#' @return vector of weighted gene expression data
#' @export
stouffer <- function(z_scores, weights) {
    zs <- rep(0, ncol(z_scores))
    sum <- 0
    for (i in 1:ncol(z_scores)) {
        sum <- sum + (z_scores[,i] * weights[i])
    }
    return(sum / sqrt(sum(weights^2)))
}

#' Function to generate data list from metadata objects
#' @param cellDF dataframe subset of cmap cell data
#' @param pertDF dataframe subset of cmap perturbation data
#' @param gctxFile filepath for a gctx file to subset
#' @param phase which phase of data are we using
#' @return a list of dataframes cellsxgenesxperturbations
#' @export
cmapRdataExtract<-function(cellDF,pertDF,sigInfoFile,gctxFile,phase=3){
  require(cmapR)
  #genes
  if(phase==3){
    data(cmapGene)
    cmapGene$gene_id=as.character(cmapGene$gene_id)
    cmapGene=cmapGene[cmapGene$feature_space=="landmark",]
    rids=cmapGene$gene_id
  }else if(phase==1 | phase ==2){
    data(cmapGeneP1)
    cmapGene=cmapGene[cmapGene$pr_is_lm==1,]
    rids=as.character(cmapGene$pr_gene_id)
  }
#sig
  sigInfo=read.csv(sigInfoFile,stringsAsFactors=FALSE,sep="\t")
  #cells
  cellIDs=unique(cellDF$cell_iname)
  #perts
  pertIDs=unique(pertDF$pert_id)
  colDat= read_gctx_meta(gctxFile, dim="col")
  rowDat= read_gctx_meta(gctxFile, dim="row")
  res=list()
  for(i in 1:length(cellIDs)){
    if(phase==3){
      cids=sigInfo[sigInfo$cell_iname==cellIDs[i]&sigInfo$pert_id%in%pertIDs,"sig_id"]
    }else{
      cids=sigInfo[sigInfo$cell_id==cellIDs[i]&sigInfo$pert_id%in%pertIDs,"sig_id"]
    }
    if(length(cids)>0){
      gctx=parse_gctx(gctxFile,cid=cids,rid=rids)
      #res[[cellIDs[i]]]=mat(gctx)
      temp=mat(gctx)
      save(temp, file=paste0(i,".rda"))
    }else{res[[cellIDs[i]]]=NA}
  }
  res
}
