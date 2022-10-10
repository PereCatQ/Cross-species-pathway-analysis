library(RCy3)
library(tibble)

# set working environment to source file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# check cytoscape connection
RCy3::cytoscapePing()

#RCy3::installApp("enhancedGraphics")
#RCy3::installApp("WikiPathways")

data <- read.table("data/normalized-data.tsv", sep = "\t", header=TRUE)
data$HumanAverage<-apply(data[,3:7], 1, mean)
data$RabbitAverage<-apply(data[,8:12], 1, mean)
data$SheepAverage<-apply(data[,13:17], 1, mean)

# Notch signaling & TFG beta
pathways <- c("WP268","WP560")

for(p in pathways) {
  command <- paste0("wikipathways import-as-pathway id=",p)
  RCy3::commandsRun(command) 
  RCy3::loadTableData(data, data.key.column = "ENSGID", table.key.column = "Ensembl")
  RCy3::setNodeCustomHeatMapChart(c("HumanAverage","SheepAverage","RabbitAverage"), slot = 2, style.name = "WikiPathways")
}

data[data$ENSGID == "ENSG00000168214",]
