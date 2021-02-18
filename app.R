library(shiny)
library(tibble)
library(VennDiagram)
library(RJSONIO)
library(STRINGdb)
library(feather)
library(grid)
library(futile.logger)

#some code from stackoverflow, which allows the app to deploy 
options(repos = BiocManager::repositories())
getOption("repos")

#Download string DB
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="")
#Download protein names for string IDs
terms<-string_db$get_aliases()
#Load list of diseases on OpenTargets
diseasePath <- "20.11_disease_list.csv"
dt<-read.table(diseasePath, sep=",", header=TRUE)
diseases<-dt[,2]

#Load CBD data
CBDPath <- "net.csv"
CBD<-read.table(CBDPath, sep=";", header=TRUE)
CBDtargets<-unique(CBD$ensembl, na.rm=TRUE)
CBDtarget <- CBDtargets[!is.na(CBDtargets)]
CBDtargets<-c(CBDtarget)
CBDstring<-as.data.frame(CBDtarget)

#Map CBD targets for StringIDs
myTargets_mapped <- string_db$map(CBDstring, "CBDtarget", removeUnmappedRows = TRUE)
myHits <- myTargets_mapped$STRING_id

#UI
ui <- fluidPage(title = "Dropdown example",
                   selectizeInput('foo', label = NULL, choices = NULL,selected = NULL,
                     options = list(maxOptions = 5, maxItems = 1)),
                   dataTableOutput("enrich"),
                   dataTableOutput("test"),
                   plotOutput("venn"),
                   plotOutput("string"))

#SERVER
server <- function(input, output, session) {
  
  #Create dropdown list of diseases
  updateSelectizeInput(session, 'foo', choices = diseases, server = TRUE, selected=NULL)
  
  #create link to diseases
  index<-reactive({
    index<-which(dt[,2] == input$foo)
    URL<-paste("https://platform-api.opentargets.io/v3/platform/public/association/filter?disease=", dt[index,1],"&size=10000&fields=target.id", sep="")
  })
  
  #create table with targets of interest
  diseaseTargets<- reactive({
    foodMarketsRaw<-fromJSON(index())
    foodMarkets<-foodMarketsRaw[['data']]
    foodMarkets[[1]][[1]][[1]]
    targetID<-sapply(foodMarkets, function(x) x[[1]][[1]])
    targetID<-c(targetID)
  })
  
  #Venn analysis
  VennTable<-reactive({
    x<-list(diseaseTargets(),CBDtargets)
    y<-calculate.overlap(x)
    VennTable<-y$a3
  })

  #Venn plot
  output$venn<-renderPlot({
    x<-list(diseaseTargets(),CBDtargets)
    y<-calculate.overlap(x)
    area1<-length(y$a1)
    area2<-length(y$a2)
    cross.area<-length(y$a3)
    grid.newpage()
    draw.pairwise.venn(area1, area2, cross.area, category = c(input$foo, "CBD"), scaled = FALSE,lty = rep("blank", 2), 
                       fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2),)
  })
  
  #String analysis
  string<-reactive({
    allTargets<-as.data.frame(diseaseTargets())
    colnames(allTargets) <- "target"
    example1_mapped <- string_db$map(allTargets, "target", removeUnmappedRows = TRUE)
    hits <- example1_mapped$STRING_id
    interactions<-string_db$get_interactions(hits)
    #get interactions of interest
    myInteractions=c()
    for (i in myHits){
      tempTable <- interactions[interactions$from == i, ]
      myInteractions <- rbind(myInteractions, tempTable) 
    } 
    
    #Find the node with maximum number of edges
    tempTable <- table(myInteractions$from)
    tempTable<-as.data.frame(table(myInteractions$from), stringsAsFactors = FALSE)
    attach(tempTable)
    newdata <- tempTable[order(-Freq),] 
    
    #rename string IDs to protein names
    x=1
    for(i in newdata$Var1){
      index<-which(terms$STRING_id==i)
      index<-index[1]
      newdata$Var1[x]=terms$alias[index]
      x=x+1
    }
    newdata
  })
  
  #Render enrichment results
  output$enrich<-renderDataTable({
    enrichment <- string_db$get_enrichment(VennTable(), category = "KEGG")
    x<-head(enrichment, n=20)
    y<-subset(x, select = c(description, preferredNames, p_value, fdr))})
    #plot Venn 
    output$string<-renderPlot({string_db$plot_network(VennTable())})
    #Plot string
    output$test<-renderDataTable(string())

}
  

  shinyApp(ui = ui, server = server)
