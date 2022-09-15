# INSTALL PACKAGES
if(!require(BiocManager)){
  install.packages("BiocManager")
  library(BiocManager)}

if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)}

if(!require(rWikiPathways)){
  BiocManager::install("rWikiPathways")
  library(rWikiPathways)}

if(!require(dplyr)){
  BiocManager::install("dplyr")
  library(dplyr)}

if(!require(tidyr)){
  BiocManager::install("tidyr")
  library(tidyr)}

if(!require(heatmaply)){
  BiocManager::install("heatmaply")
  library(heatmaply)}

###########################################################

downloadPathwayArchive(date="20220810", organism="Homo sapiens", format = "gmt")
wp <- list.files(pattern = "wikipathways")
wp <- clusterProfiler::read.gmt(wp)
colnames(wp)[c(1,2)] <- c("pathway", "entrezgene")
wp <- wp[wp$pathway %in% names(table(wp$pathway))[table(wp$pathway) >= 1],]
colnames(wp)[2] <- "entrezgene_id"
wp <- separate(data = wp, col = pathway, into = c("name", "version", "id", "species"), sep = "%")

#######################################################################
# TODO: replace with all differently active pathways (read in from file)
#######################################################################
pathways2<-read.table("HumanVRabbit_SheepVRabbit.txt", sep = "\t", header = T)
pathways<-pathways2[,c(1)]
#pathways <- c("WP5049","WP384", "WP5192", "WP3844", "WP268","WP2526", "WP4874", "WP4357" )

wp.filt <- wp[wp$id %in% pathways,]

jaccard <- matrix(ncol=length(pathways), nrow=length(pathways))
overlap <- matrix(ncol=length(pathways), nrow=length(pathways))
for(i in 1:length(pathways)) {
  for(j in i:length(pathways)) {
    genes1 <- wp.filt[wp.filt$id == pathways[i],"entrezgene_id"]
    genes2 <- wp.filt[wp.filt$id == pathways[j],"entrezgene_id"]
    intersection = length(intersect(genes1, genes2))
    union = length(genes1) + length(genes2) - intersection
    score <- (intersection/union)
    jaccard[i,j] <- score
    jaccard[j,i] <- score
    overlap[i,j] <- intersection
    overlap[j,i] <- intersection
  }
}

colnames(jaccard) <- pathways
rownames(jaccard) <- pathways

colnames(overlap) <- pathways
rownames(overlap) <- pathways
diag(overlap) <- 0

heatmaply(jaccard,main="Similarity plot",plot_method = "plotly", colors = viridis(n = 256,  option = "magma"))
heatmaply(overlap,main="Shared genes (diagonal = 0)",plot_method = "plotly", colors = viridis(n = 256,  option = "magma"), dendrogram = "none")

# SIMILARITY NETWORK
#######################################################################
# TODO: adapt cutoff to higher similarity score when using all significant pathways
#######################################################################
sim.cutoff <- 0.05 #between 0 and 1
adj.matrix <- jaccard
adj.matrix[adj.matrix < sim.cutoff] <- 0 
diag(adj.matrix) <- 0

graph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, mode="undirected", weighted = TRUE)
suid <- RCy3::createNetworkFromIgraph(graph, title = "Similarity network", collection = "Similarity")
RCy3::loadTableData(wp.filt, data.key.column = "id", table.key.column = "name")
RCy3::analyzeNetwork(directed=FALSE)
RCy3::toggleGraphicsDetails()
RCy3::copyVisualStyle("default","similarity")
RCy3::setVisualStyle("similarity")
RCy3::setNodeShapeDefault("ellipse", style.name="similarity")
RCy3::lockNodeDimensions(TRUE, style.name="similarity")
RCy3::setNodeSizeMapping('Degree', c(1,15), c(30,80), mapping.type = "c", style.name = "similarity")

# GENE OVERLAP NETWORK
overlap.cutoff <- 0 # increase if you want to remove pathways only sharing few genes
adj.matrix2 <- overlap
adj.matrix2[adj.matrix2 < overlap.cutoff] <- 0 

graph2 <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix2, mode="undirected", weighted = TRUE)
suid <- RCy3::createNetworkFromIgraph(graph2, title = "Shared Genes network", collection = "Overlap")
RCy3::loadTableData(wp.filt, data.key.column = "id", table.key.column = "name")
RCy3::analyzeNetwork(directed=FALSE)
RCy3::toggleGraphicsDetails()
RCy3::copyVisualStyle("default","overlap")
RCy3::setVisualStyle("overlap")
RCy3::setNodeShapeDefault("ellipse", style.name="overlap")
RCy3::lockNodeDimensions(TRUE, style.name="overlap")
RCy3::setNodeSizeMapping('Degree', c(1,15), c(30,80), mapping.type = "c", style.name = "overlap")
RCy3::setEdgeLineWidthMapping('weight', c(1,max(overlap)), c(1,5), mapping.type="c", style.name = "overlap")
