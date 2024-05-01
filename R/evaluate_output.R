#Same as evaluateoutput1 and 2, includes functions from cmap_query.R all in one script

evaluate_results = function(trueSigs, imputedSigs, foldAssignMat, RWF=.99){
  cells = rownames(foldAssignMat)
  drugs = colnames(foldAssignMat)
  names(imputedSigs)=cells
  evalMethods = c("corScore","corCmapScorePos.99","corCmapScoreNeg.99")
  scoreMat = array(NA,dim=c(length(cells),length(drugs),length(evalMethods)),dimnames=list(cells,drugs,evalMethods))
  for(cell in 1:length(cells)){
    existing_drugs = names(which(foldAssignMat[cell,]!=0))
    for(drug in 1:length(drugs)){
      if(foldAssignMat[cell,drug]!=0){
        trueSigs[[cell]]=data.frame(trueSigs[[cell]])
        row.names(trueSigs[[cell]])=drugs
        true_sig = as.numeric(trueSigs[[cell]][drug,])
        names(true_sig)=names(trueSigs[[1]][1,])
        imputed_sig = imputedSigs[[cell]][drug,]
        names(imputed_sig)=names(trueSigs[[1]][1,])
        scoreMat[cell,drug,"corScore"] = cor(true_sig,imputed_sig,method="spearman")

        toQueryMat = trueSigs[[cell]][-drug,]

        query_result_true = run_regular_cmap_query_consensus_single_cell(true_sig,toQueryMat)
        query_result_imputed = run_regular_cmap_query_consensus_single_cell(imputed_sig,toQueryMat)
        scoreMat[cell,drug,"corCmapScorePos.99"] = weighted_rank_cor(query_result_true,query_result_imputed,TRUE,RWF)
        scoreMat[cell,drug,"corCmapScoreNeg.99"] = weighted_rank_cor(query_result_true,query_result_imputed,FALSE,RWF)

      }
    }
    save(scoreMat,file="scoreMat.rda")
  }
  return(scoreMat)
}


run_regular_cmap_query_consensus_single_cell <-function(signature=c(), exp.matrix, num_up = 50, q.up=NA,q.down=NA) {
  if(length(signature)>0){
    #Need to select a set of genes to be the query upregulated and query downregulated sets.
    q.up = intersect(names(signature)[order(signature,decreasing=TRUE)][1:num_up],names(which(signature>0)))
    q.down = intersect(names(signature)[order(signature)][1:num_up],names(which(signature<0)))
  }

  if(length(q.up)==0 & length(q.down)==0){return(NA)}

  all_drugs = rownames(exp.matrix)

  query_results = c()
  for(drug in all_drugs){
    ref_sig = as.numeric(exp.matrix[drug,])
    names(ref_sig)=names(trueSigs[[1]][1,])

    if(any(is.na(ref_sig))){
      query_results[drug,cell] = NA
    } else {
      query_results[drug] = WTCS(q.up,q.down,ref_sig)
    }
  }
  return(query_results)
}

# Given two vectors of CMAP query scores for a each signature (where positive scores indicate positive connectivity and vice-versa), calculate a weighted rank correlation.
weighted_rank_cor = function(query_result_true,query_result_imputed,positive=TRUE,RWF=.99){
  x = query_result_true
  y = query_result_imputed
  if(positive){ #the call to rw uses "rank", which will rank from low to high.  Switching signs here makes the the thing with the highest connectivity score have rank 1.
    x = 0-x
    y = 0-y
  }
  wfunc=function(x,n) weight_norm(x,n,f=RWF)
  return(rw(x,y,wfunc)/rw(x,x,wfunc))
}


# Weighted Spearman Rank Correlation
# From Shieh et al "Rank Tests for Independence" (2000)
#
rw <- function(x, y, wfunc=(function (x,n) 1/n)) {
  xr <- rank(x)
  yr <- rank(y)
  n <- length(xr)
  s <- 0
  for (i in 1:n) {
    weight <- wfunc(xr[i],n)
    d <- weight * (xr[i] - ((n+1)/2)) * (yr[i] - ((n+1)/2))
    s <- s+d
  }
  return(s)
}

weight_norm <- function(r,n,f=0.99) {
  return(2*dnorm(r-0.5,0,n-(n*f)))
}


#From cmap_query.R, compute weighted connectivity score for query gene set pair (q.up, q.down) and reference signature r
WTCS <- function(q.up, q.down, r){
  r.ordered = r[order(r, decreasing = TRUE)]
  ES.up = ES(r.ordered, q.up)
  ES.down = ES(r.ordered, q.down)
  if(sign(ES.up) != sign(ES.down)){
    return ((ES.up - ES.down)/ (2.0))
  }else{
    return (0)
  }
}

#From cmap_query.R
ES <- function(gene.list, gene.set, weighted.score.type = 0){
  correl.vector <- gene.list
  gene.list <- names(gene.list)
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
  no.tag.indicator <- 1 - tag.indicator
  correl.vector <- abs(correl.vector)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  if(sum.correl.tag != 0){
    norm.tag <- 1.0/sum.correl.tag
  }else{
    norm.tag <- 0
  }
  norm.no.tag <- 1.0/(length(gene.list) - length(gene.set))
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES){
    ES <- signif(max.ES, digits = 5)
  } else{
    ES <- signif(min.ES, digits = 5)
  }
  return (ES)
}
