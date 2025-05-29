####################################
#Name: cmapFilter
#Goal: To help filter files to find data
#      based on ontologies
#
####################################

#'find data based on cell information
#'@import shiny
#'@export
runCmapCellFilter<-function(cellFile){
	myfile=system.file(package="larimaR","shinyApps/cmapCell.R")
	source(myfile)
	cmapCellFilter(cellFile)
}

#'find data based on perturbation information
#'@import shiny
#'@export
runCmapPertFilter<-function(pertFile){
	myfile=system.file(package="larimaR","shinyApps/cmapPert.R")
	source(myfile)
	cmapPertFilter(pertFile)
}

#'show available functions
#'@export
showAll=function(){
	cat(paste("Functions:","runCmapPertFilter", "runCmapCellFilter",sep="\n"))
}
