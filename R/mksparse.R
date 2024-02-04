#' Function to create folds of cmap dataset
#' @param cmapDS cmap dataset array cellsxgenesxperts
#' @param nfolds number of folds to create
#' @param oneCellPercent numeric percent of cells that a fold can have
#' @param onePertPercent numeric percent of perturbagens a fold can have
#' @return foldDS array conforming to cmapDS that has fold information
#' @export
mkCmapFold<-function(cmapDS, nfolds, oneCellPercent, onePertPercent){
  ncells=length(cmapDS)
  nperts=nrow(cmapDS[[1]])
  ngenes=ncol(cmapDS[[1]])
  cnames=names(cmapDS)
  gnames=colnames(cmapDS[[1]])
  pnames=rownames(cmapDS[[1]])
  res=matrix(ncol=nperts, nrow=ncells, data=sample(1:nfolds,ncells*nperts, replace=TRUE))
  res
}

#' Function to create folds of cmap dataset
#' @param cmapFold cmap fold dataset from mkCmapFold
#' @return bool does the fold dataset conform to restrictions?
ckCmapFolds<-function(cmapFold, nfolds,oneCellPercent, onePertPercent){
  for(i in 1:nfolds){
    cSums=apply(cmapFold,2,function(x){sum(x==i)})
    rSums=apply(cmapFold,1,function(x){sum(x==i)})
    if(any(100*(cSums/nrow(cmapFold))>=onePertPercent)){return(FALSE)}
    if(any(100*(rSums/ncol(cmapFold))>=oneCellPercent)){return(FALSE)}
  }
  return(TRUE)
}

#' Function to add missingness
#' @param percentMissing percent of missingness to add
#' @return matrix of data with missingness added
addMissing<-function(cmapDS,cmapFold,percentMissing ){
  
  return(res)
}
