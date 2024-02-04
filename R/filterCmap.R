#' Function to classify points
#' @return dataframe discribing points
#' @param quorumPercent the percent of data above tox threshhold to keep any
#' @export
labelSubset<-function(dataOb,quorumPercent=50){
  prepLab<-function(lab){
    paste(lab$dose,lab$unit, lab$doseTime,lab$doseTimeUnit, sep="<br>")
  }

  if(quorumPercent>1){quorumPercent=quorumPercent/100}
  #rules
  mydat=data.frame(
    dosageAmt=as.numeric(as.character(dataOb$dosage["dose",])),
    dosageTime=as.factor(as.character(dataOb$dosage["doseTime",])),
    phase=as.factor(as.character(dataOb$phase)),
    viability=dataOb$achilles,
    text=sapply(1:ncol(dataOb$dosage), function(x){prepLab(dataOb$dosage[,x])})
  )

  times=sort(as.numeric(as.character(unique(mydat$dosageTime))),decreasing=TRUE)
  doses=sort(unique(mydat$dosageAmt),decreasing=TRUE)

  slopes=c()
  pvals=c()

  for(i in 1:length(times)){
    slopes[i]=NA
    pvals[i]=NA
    temp=mydat[mydat$dosageTime==times[i]&!is.na(mydat$dosageTime),]
    if(nrow(temp)>1&length(unique(temp$dosageAmt))>1){
      tlm=lm(viability~dosageAmt,data=temp)
      slopes[i]=summary(tlm)$coefficients["dosageAmt",1]
      pvals[i]=summary(tlm)$coefficients["dosageAmt",4]
    }
  }

if(length(times)==1){
  #one
  mydat$mult=0
}else{mydat$mult=1}

#tox screen
mydat$nonToxic=mydat$viability>=-3

#keep
mydat$keep=TRUE

#apply tox screen
mydat$keep=ifelse(mydat$nonToxic,mydat$keep, FALSE )

#apply NA screen
mydat$keep=ifelse(!is.na(mydat$viability),mydat$keep, FALSE )

targetTime=NA
targetDose=NA

#start at bottom right -times sorted highest prioritize up over left
for(j in 1:length(doses)){
  for(i in 1:length(times)){
    temp=subset(mydat,mydat$dosageAmt==doses[j]&mydat$dosageTime==times[i])
    if(nrow(temp)>0){
      if((sum(temp$nonToxic)/nrow(temp))>quorumPercent){
        targetTime=times[i]
        targetDose=doses[j]
        break
      }
    }
  }
  if(!is.na(targetTime)){break}
}

if(!is.na(targetTime)&!is.na(targetDose)){
  mydat$keep=(mydat$keep & mydat$dosageTime==targetTime & mydat$dosageAmt==targetDose)
}
#if no acceptable group, keep everything
if(is.na(targetTime)|is.na(targetDose)){
  mydat$keep=TRUE
  #apply NA screen
  mydat$keep=ifelse(!is.na(mydat$viability),mydat$keep, FALSE )
  print(paste("Default for :",dataOb$pertName))
}
return(mydat)
}

#' Function to Plot Cell viability by dosage Filtered
#' @return plot of cell viability by dosage
#' @export
plotCmapFilter<-function(dataOb, quorumPercent=.5){
  require(plotly)
  mydat=labelSubset(dataOb,quorumPercent)
  fig=plot_ly(data = mydat, x = ~dosageAmt,
    y = ~viability, color=~dosageTime, symbol=~keep,text=~text,
    type="scatter", mode="markers",
     marker = list(size = 20)
  )%>%layout(title = paste('Viability vs Dose:',dataOb$pertName))

  times=unique(mydat$dosageTime)
  for(i in 1:length(times)){
    temp=mydat[mydat$dosageTime==times[i],]
    tlm=lm(viability~dosageAmt,data=temp)
    fig=fig%>%add_lines(x= ~dosageAmt, y = fitted(tlm),data=temp)
  }
  fig
}

#' Function to extract filtered data
#' @return dataframe
#' @export
procFilter<-function(dataOb, quorumPercent=.5){
  mydat=labelSubset(dataOb,quorumPercent)
  expNames=row.names(mydat)[mydat$keep]
  if(length(expNames>0)){
    expDat=data.frame(dataOb$expData[,expNames])
    require(cmap2020)
    res=data.frame(stouffer(expDat,weight_signature(expDat)))
  }else{res=data.frame(matrix(ncol=1,nrow=978))}
  names(res)=dataOb$pertName
  row.names(res)=row.names(dataOb$expDat)
  res
}

#' Function to classify points
#' @return dataframe discribing points
#' @export
labelSubsetMax<-function(dataOb){
  prepLab<-function(lab){
    paste(lab$dose,lab$unit, lab$doseTime,lab$doseTimeUnit, sep="<br>")
  }

  #rules
  mydat=data.frame(
    dosageAmt=as.numeric(as.character(dataOb$dosage["dose",])),
    dosageTime=as.factor(as.character(dataOb$dosage["doseTime",])),
    phase=as.factor(as.character(dataOb$phase)),
    viability=dataOb$achilles,
    text=sapply(1:ncol(dataOb$dosage), function(x){prepLab(dataOb$dosage[,x])})
  )

  times=sort(as.numeric(as.character(unique(mydat$dosageTime))),decreasing=TRUE)
  doses=sort(unique(mydat$dosageAmt),decreasing=TRUE)

#keep
mydat$keep=TRUE

#apply NA screen
mydat$keep=ifelse(!is.na(mydat$viability),mydat$keep, FALSE )

targetTime=NA
targetDose=NA


targetTime=times[1]
targetDose=doses[1]

mydat$keep=(mydat$keep & mydat$dosageTime==targetTime & mydat$dosageAmt==targetDose)

return(mydat)
}
