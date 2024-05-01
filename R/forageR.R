####################################
#Name: forageR
#Goal: To help filter files to find data
#      based on ontologies
#
####################################

#'find data based on cell information
#'@import shiny
#'@export
runforageRCell<-function(){
	myfile=system.file(package="larimaR","shinyApps/forageRAppCell.R")
	source(myfile)
	forageRCell()
}

#'find data based on perturbation information
#'@import shiny
#'@export
runforageRPert<-function(){
	myfile=system.file(package="larimaR","shinyApps/forageRAppPert.R")
	source(myfile)
	forageRPert()
}

#'show available functions
#'@export
showAll=function(){
	cat(paste("Functions:","runforageRCell", "runforageRPert",sep="\n"))
}
