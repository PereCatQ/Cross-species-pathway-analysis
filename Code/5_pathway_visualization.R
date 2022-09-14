library(RCy3)
library(tibble)

# set working environment to source file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# check cytoscape connection
RCy3::cytoscapePing()

#RCy3::installApp("enhancedGraphics")
#RCy3::installApp("WikiPathways")

DataWithAverages <- read.delim("DataWithAverages.txt")
DataWithAverages <- tibble::rownames_to_column(DataWithAverages, "Ensembl")
pathways <- c("WP560","WP2884", "WP4357", "WP3", "WP268","WP4874")

for(p in pathways) {
  command <- paste0("wikipathways import-as-pathway id=",p)
  RCy3::commandsRun(command) 
  RCy3::loadTableData(DataWithAverages, data.key.column = "Ensembl", table.key.column = "Ensembl")
  RCy3::setNodeCustomHeatMapChart(c("HumanAverage","SheepAverage","RabbitAverage"), slot = 2, style.name = "WikiPathways")
}
