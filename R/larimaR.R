#################################################
#LarimaR: Functions to help analyze l1000 data
#
#
#################################################


#' Function to get broad ids from perturbagen name
#' @param pertName the name of the perturbagen
#' @param sigInfo sig file data
#' @param phase which phase of data are we using
#' @return a vector of broad ids
#' @export
getBRDS<-function(pertName,sigInfo,phase=3){
  if(phase==3){
    temp=subset(sigInfo,sigInfo$cmap_name==pertName)
  }else{
    temp=subset(sigInfo,sigInfo$pert_iname==pertName)
  }
  brds=unique(temp$pert_id)
  brds
}


#' Function to generate data list from metadata objects
#' @param cell cmap cell
#' @param pert cmap perturbation
#' @param sigInfo data for the gctx file to subset
#' @return dataframe of genesxperturbations
#' @export
cmapDoseExtract<-function(cell,pert,sigInfo,gctxFile,phase=3){
  require(cmapR)
  require(cmap2020)
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
  colDat= read_gctx_meta(gctxFile, dim="col")
  rowDat= read_gctx_meta(gctxFile, dim="row")
  res=list()
  if(phase==3){
    cids=sigInfo[sigInfo$cell_iname==cell&sigInfo$pert_id==pert,"sig_id"]
  }
  else{
    cids=sigInfo[sigInfo$cell_id==cell&sigInfo$pert_id==pert,"sig_id"]
  }
  if(length(cids)>0){
    gctx=parse_gctx(gctxFile,cid=cids,rid=rids)
    #res[[cellIDs[i]]]=mat(gctx)
    temp=mat(gctx)
    if(phase==3){
      row.names(temp)=cmapGene$gene_symbol}else{
        row.names(temp)=cmapGene$pr_gene_symbol
      }
  }else{temp=NA}
  return(temp)
}

#' Function to get broad ids from perturbagen name
#' @param geneP3 the name of the gene in phase3
#' @return a phase1/2 gene name
#' @export
harmonizeGene<-function(geneP3){
  if(geneP3=="ELP1"){return("IKBKAP")}#https://medlineplus.gov/genetics/gene/elp1/#synonyms
  if(geneP3=="TOMM70"){return("TOMM70A")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=TOMM70
  if(geneP3=="WASHC5"){return("KIAA0196")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=WASHC5
  if(geneP3=="RXYLT1"){return("TMEM5")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=RXYLT1
  if(geneP3=="TENT4A"){return("PAPD7")}#http://www.proteinatlas.org/ENSG00000112941-TENT4A
  if(geneP3=="KHDC4"){return("KIAA0907")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=KHDC4
  if(geneP3=="WASHC4"){return("KIAA1033")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=WASHC4
  if(geneP3=="CEMIP2"){return("TMEM2")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=CEMIP2
  if(geneP3=="DMAC2L"){return("ATP5S")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=DMAC2L
  if(geneP3=="HDGFL3"){return("HDGFRP3")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=HDGFL3
  if(geneP3=="CARMIL1"){return("LRRC16A")}#https://www.proteinatlas.org/ENSG00000079691-CARMIL1
  if(geneP3=="MINDY1"){return("FAM63A")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=MINDY1
  if(geneP3=="COQ8A"){return("ADCK3")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=COQ8A
  if(geneP3=="SQOR"){return("SQRDL")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=SQOR
  if(geneP3=="PRUNE1"){return("PRUNE")}
  if(geneP3=="CIAO3"){return("NARFL")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=CIAO3
  if(geneP3=="JPT2"){return("HN1L")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=JPT2
  if(geneP3=="STIMATE"){return("TMEM110")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=STIMATE
  if(geneP3=="DIPK1A"){return("FAM63A")}#https://www.genecards.org/cgi-bin/carddisp.pl?gene=DIPK1A
  return(NA)
}

#' Function to get a measure of cell viability based on the achilles model
#' @param myvec a vector of gene expression from l1000
#' @return a measure of cell viability based on the achilles model
#' @export
cmapTox<-function(myvec){
  data(achModList)
  cur=achModList[["intercept"]]
  for(i in 1:length(myvec)){
    coef=NA
    coef=achModList[[names(myvec)[i]]]
    if(is.null(coef)){coef=achModList[[harmonizeGene(names(myvec)[i])]]}
    if(is.na(coef)){error("missing coef")}
    cur=cur+(coef*myvec[i])
  }
  return(as.numeric(cur))
}

#' Function to get harmonize genenames
#' @param genes vector of phase3 gene names
#' @return a vector of phase1/2 gene names
#' @export
harmonizeGenes<-function(genes){
  newnames=c()
  for(i in 1:nrow(genes)){
    test=achModList[[rownames(genes)[i]]]
    if(is.null(test)){newnames[i]=harmonizeGene(rownames(genes)[i])}else{
      newnames[i]=rownames(genes)[i]
    }
  }
  rownames(genes)=newnames
  genes
}

#' Function to combine data across broad ids
#' @param broadIDS character vector of Broad ids
#' @return a dataframe of gene expression for 1 cell mult Broad ids
#' @export
mergeBroad<-function(broadIDS,cell,sigInfo,gctxFile,phase=3){
  dat=cmapDoseExtract(cell,broadIDS[1],sigInfo,gctxFile,phase=phase)
  if(length(broadIDS)>1){
    for(i in 2:length(broadIDS)){
      tdat=cmapDoseExtract(cell,broadIDS[i],sigInfo,gctxFile,phase)
      if(!is.null(nrow(tdat))){dat=cbind(dat,tdat)}
    }
  }
  dat
}

#' Function to get broad ids from perturbagen name
#' @param sigID the name of the perturbagen
#' @param sigInfo sig file data
#' @return a vector of broad ids
#' @export
getDosage<-function(sigID,sigInfo,phase=3){
  temp=subset(sigInfo,sigInfo$sig_id==sigID)
  if(nrow(temp)>1){warning("multiple sig id records found\n")}
  if(phase>=2){
    return(list(unit=temp[["pert_dose_unit"]][1],
      dose=temp[["pert_dose"]][1],
      doseTime=temp[["pert_time"]][1],
      doseTimeUnit=temp[["pert_time_unit"]][1]
    ))}else{
      #split on space take first argument for now, lets see what happens
      return(list(unit=strsplit(temp$pert_idose[1]," ")[[1]][2],
        dose=strsplit(temp$pert_idose[1]," ")[[1]][1],
        doseTime=strsplit(temp$pert_itime[1]," ")[[1]][1],
        doseTimeUnit=strsplit(temp$pert_itime[1]," ")[[1]][2]
      ))
    }
}

#' Function to run dosage pipeline
#' @return a list of data
#' @export
getDosagePipeline<-function(pertName,sigInfo,gctxFileLocation,cell,curphase=3){
  broadIDS=getBRDS(pertName,sigInfo,phase=curphase)
  if(length(broadIDS)==0){return(NA)}
  expData=mergeBroad(broadIDS,cell,sigInfo,gctxFileLocation,curphase)
  dosages=NA
  achilles=NA
  phase=NA
  if(length(colnames(expData))>0){
    dosages=sapply(colnames(expData),function(x){getDosage(x,sigInfo,curphase)})
    achilles=apply(expData,2,cmapTox)
    phase=rep(curphase,ncol(expData))
  }
  return(list(pertName=pertName,broadIDS=broadIDS,
     cell=cell, expData=expData,dosages=dosages,achilles=achilles,phase=phase))
}

#' Function to merge 2 data extracts
#' @return a list of data
#' @param p1 output from dosage pipeline
#' @param p2 output from dosage pipeline
#' @return list of data
#' @export
mergeDosageExtracts<-function(p1,p2){
  if(class(p1)=="logical" & class(p2)=="logical"){return(NA)}
  if(class(p1)=="logical"){return(p2)}
  if(class(p2)=="logical"){return(p1)}
  if(class(p1$expData)[1]=="logical"){return(p2)}
  if(class(p2$expData)[1]=="logical"){return(p1)}

  p3=p1
  p3$broadIDS=c(p1$broadIDS,p2$broadIDS)
  p3$expData=cbind(p1$expData,p2$expData)
  p3$dosages=cbind(p1$dosages,p2$dosages)
  p3$achilles=c(p1$achilles,p2$achilles)
  p3$phase=c(p1$phase,p2$phase)
  p3
}


#' Function to run dosage pipeline on phase 1 and 2
#' @return a list of data
#' @param sigInfo1 sig info for phase 1
#' @param sigInfo2 sig info for phase 1
#' @param gctxFile1Location gctx for phase 1
#' @param gctxFile2Location gctxfor phase 2
#' @export
getPhase1and2<-function(pertName,sigInfo1,sigInfo2,gctxFile1Location, gctxFile2Location,cell){
  p1=getDosagePipeline(pertName,sigInfo1,gctxFile1Location,cell,curphase=1)
  p2=getDosagePipeline(pertName,sigInfo2,gctxFile2Location,cell,curphase=2)
  return(mergeDosageExtracts(p1,p2))
}



#' Function to Plot Cell viability by dosage
#' @return plot of cell viability by dosage
#' @export
plotCmapDosage<-function(dataOb){
  t2 <- list(
  size = 36,
  color = "black")

  t3 <- list(
  size = 24,
  color = "black")

  prepLab<-function(lab){
    paste(lab$dose,lab$unit, lab$doseTime,lab$doseTimeUnit, sep="<br>")
  }
  require(plotly)
  mydat=data.frame(
    dosageAmt=as.numeric(as.character(dataOb$dosage["dose",])),
    dosageTime=as.factor(as.character(dataOb$dosage["doseTime",])),
    dosageTimeNum=as.numeric(as.character(dataOb$dosage["doseTime",])),

    phase=as.factor(as.character(dataOb$phase)),
    viability=dataOb$achilles,
    text=sapply(1:ncol(dataOb$dosage), function(x){prepLab(dataOb$dosage[,x])})
  )
  fig=plot_ly(data = mydat, x = ~dosageAmt,
    y = ~viability, color=~dosageTime, symbol=~phase,text=~text,
    type="scatter", mode="markers",
     marker = list(size = 20)
  )%>%layout(title = list(text=paste('Viability vs Dose:',dataOb$pertName)
    ,y=.9,font=t2),
  xaxis = list(title = list(text ='Dose Amount', font = t2)),
  yaxis = list(title = list(text ='Cell Viability', font = t2)),
  legend=list(font=t3,title=list(text="Color is Time\n Symbol is Phase"))
)

  times=unique(mydat$dosageTimeNum)
  for(i in 1:length(times)){
    temp=mydat[mydat$dosageTimeNum==times[i],]
    tlm=lm(viability~dosageAmt,data=temp)
    fig=fig%>%add_trace(x= ~dosageAmt, y = fitted(tlm),
      data=temp,line=list(width=4),mode = 'lines',marker=list(opacity=0))
    }
  fig
}
