####################################
#Name: cmapPertFilter
#Goal: To help filter files to find data
#      based on ontologies
#
####################################

########################
#load libraries
########################
require(shiny)
require(ggplot2)
require(DT)
require(ontologyIndex)
########################


########################
#functions
########################

relFunc<-function(mydf,myvar,myrel,mypattern){
  if(myrel%in%relations){
    inInd=NA
    if(myrel=="contains"){
      inInd=regexpr(mydf[[myvar]],pattern=mypattern,ignore.case = TRUE)>0
    }
    if(myrel=="does not contain"){
      inInd=regexpr(mydf[[myvar]],pattern=mypattern,ignore.case = TRUE)<0
    }
    if(myrel=="=="){
      inInd=(mydf[[myvar]]==mypattern)
    }
    if(myrel=="!="){
      inInd=(mydf[[myvar]]!=mypattern)
    }
    if(myrel==">"){
      inInd=(mydf[[myvar]]>mypattern)
    }
    if(myrel==">="){
      inInd=(mydf[[myvar]]>=mypattern)
    }
    if(myrel=="<"){
      inInd=(mydf[[myvar]]<mypattern)
    }
    if(myrel=="<="){
      inInd=(mydf[[myvar]]<=mypattern)
    }
    if(myrel=="not NA"){
      inInd=(!is.na(mydf[[myvar]]))
    }

    res=list();
    res[["mydfIn"]]=mydf[inInd,]
  	res[["mydfOut"]]=mydf[!inInd,]
  	return(res)
  }else{stop("That is not a valid relation")}
}

relText=list(
  "contains"='inInd=regexpr(mydf[["myvar"]],pattern="mypattern",ignore.case = TRUE)>0',
  "does not contain"='inInd=regexpr(mydf[["myvar"]],pattern="mypattern",ignore.case = TRUE)<0',
  "=="='inInd=(mydf[["myvar"]]=="mypattern")',
  "!="='inInd=(mydf[["myvar"]]!="mypattern")',
  ">"='inInd=(mydf[["myvar"]]>"mypattern")',
  ">="='inInd=(mydf[["myvar"]]>="mypattern")',
  "<"='inInd=(mydf[["myvar"]]<"mypattern")',
  "<="='inInd=(mydf[["myvar"]]<="mypattern")',
  "not NA"='inInd=(!is.na(mydf[["myvar"]]))'
  )

relTextFunc<-function(myreltext,mydf,myvar,mypattern,myrel){
  tempText=relText[[myrel]]
  tempText=sub("myvar",myvar,tempText)
  tempText=sub("mypattern",mypattern,tempText)
  tempText=sub("myrel",myrel,tempText)
  tempText=sub("mydf",deparse(substitute(mydf)),tempText)
  tempText
}
#text mining help from here
#http://www.sthda.com/english/wiki/text-mining-and-word-cloud-fundamentals-in-r-5-simple-steps-you-should-know

mkcloud<-function(mydf, myvar,mypalette){
	require(wordcloud)
	require(tm)
	mycorp=Corpus(VectorSource(mydf[[myvar]]))
	mydocs=tm_map(mycorp, removeWords, stopwords("english"))

	dtm=TermDocumentMatrix(mydocs)
	m=as.matrix(dtm)
	v=sort(rowSums(m),decreasing=TRUE)
	d=data.frame(word = names(v),freq=v)
	wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=50, random.order=FALSE, rot.per=0.35,scale=c(6,1))
          #,
          #colors=brewer.pal(8, mypalette)[3:8])
}


filterCmap<-function(mydf,myvar,mypattern){
	res=list();
	inInd=regexpr(mydf[[myvar]],pattern=mypattern,ignore.case = TRUE)>0
	res[["mydfIn"]]=mydf[inInd,]
	res[["mydfOut"]]=mydf[!inInd,]
	res
}

mkcloudFilter<-function(mylist,myvar){
  par(mfrow=c(1,2))
  mkcloud(mylist[[1]],myvar,"Blues")
  mkcloud(mylist[[2]],myvar,"Reds")
}

listElementBullet<-function(myListElement, myListElementName){
  myhtml=paste0("<br>",myListElementName,"<ul>")
  if(length(myListElement)>0){
    for(i in 1:length(myListElement)){
      myhtml=paste0(myhtml,"<li>",myListElement[[i]],"</li>")
    }
  }
  myhtml=paste0(myhtml,"</ul><br>")
  return(myhtml)
}

listBullet<-function(mylist,mylistnames){
  myhtml=""
  for(i in 1:length(mylistnames)){
    myhtml=paste0(myhtml,listElementBullet(mylist[[i]],mylistnames[i]))
  }
  return(myhtml)
}


########################
#shiny
########################
cmapPertFilter = function(PertDataFile,oboFile=NA) {
  #Setup:
    cmapPertdata<<-read.csv(PertDataFile,sep="\t")
    origdata<<-cmapPertdata
    relations<<-c("contains","does not contain","==","!=",">",">=","<","<=","not NA")

########################
#load data - for testing now,
# may query databases and such later
########################

  #load ontology
  if(!is.na(oboFile)){
    myont=get_ontology(oboFile, propagate_relationships = "is_a",
	   extract_tags = "everything")
   }
  mylog=c()

    #start shiny app config
    shinyApp(
        ##########
        #Start UI Config
        ##########
        ui = fluidPage(
            titlePanel("cmapPertFilter"),
            sidebarLayout(position = "left",
                  sidebarPanel(width=2,
                    selectInput("filterVar","Filter Column",names(cmapPertdata)),
                    selectInput("filterRel","Relation",relations),
                    textInput("filterPattern", "Filter Pattern"),
                    selectInput("plotVar","Plot Column",names(cmapPertdata)),
                    actionButton("filterButton", "Filter"),
                    actionButton("resetButton", "Reset"),
                    textInput("downloadName", "Download Name"),
                    downloadButton('downloadData', 'Download Data')
                                       ),
              mainPanel("",width=10,
               tabsetPanel(
                  tabPanel("Data",h3("Current Data"),
                  dataTableOutput('mytable')),

                  tabPanel("Word Clouds",h3("Word Clouds: Keep and Remove"),
                  plotOutput('wcPlot', height = "600px")),

                  tabPanel("Synonyms",h3("Synonyms from Disease Ontology"),
                  htmlOutput('syns')),

                  tabPanel("Log",h3("Log:"),
                  h3(htmlOutput('myLog')))

                  )
                                    )
                          )
            )
        ,

        ####################
        #Start Server Config
        ####################
        server = function(input, output, session) {
            #set up output
         observeEvent(input$filterButton, {
               #subs=filterForage(cmapPertdata,input$filterVar,input$filterPattern)
               subs=relFunc(cmapPertdata,input$filterVar,input$filterRel,input$filterPattern)
               #table
               output$mytable = renderDataTable({subs[[1]]},options = list(scrollX = TRUE,pageLength=5))
               #plots
               output$wcPlot<-renderPlot(mkcloudFilter(subs,input$plotVar))
               if(!is.na(oboFile)){
	                #Synonyms
                  syns=myont$synonym[grep(x=myont$name, pattern=input$filterPattern)]
                  synsNames=myont$name[grep(x=myont$name, pattern=input$filterPattern)]
	                names(syns)=synsNames
                  output$syns= renderUI({HTML(listBullet(syns,synsNames))})
                }
               cmapPertdata<<-subs[[1]]
               #logtext=paste("cmapPertdata", input$filterVar, input$filterRel, input$filterPattern)
               logtext="inInd=NA"
               mylog<<-c(mylog,logtext)
               logtext=relTextFunc(relText[[input$filterRel]],cmapPertdata,input$filterVar,input$filterPattern,input$filterRel)
               mylog<<-c(mylog,logtext)
               logtext="cmapPertdata=cmapPertdata[inInd,]"
               mylog<<-c(mylog,logtext)
               output$myLog=renderUI({HTML(paste(mylog,collapse="<br>"))})
          })

          observeEvent(input$resetButton, {
            cmapPertdata<<-origdata
            output$mytable = renderDataTable({cmapPertdata},options = list(scrollX = TRUE,pageLength=5))
            output$wcPlot<-renderPlot(plot.new())
            mylog<<-c()
            output$myLog=renderUI({HTML(paste(mylog,collapse="<br>"))})


          })
          output$downloadData <- downloadHandler(
             filename = function() { paste(input$downloadName, Sys.Date(), ".csv", sep="")},
             content = function(file) {
                 write.csv(cmapPertdata, file,row.names=F)
             }
          )

        }
        )
}
