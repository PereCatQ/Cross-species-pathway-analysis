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

if(!require(qusage)){
  BiocManager::install("qusage")
  library(qusage)}

if(!require(plyr)){
  BiocManager::install("plyr")
  library(plyr)}

if(!require(dplyr)){
  BiocManager::install("dplyr")
  library(dplyr)}

if(!require(biomaRt)){
  BiocManager::install("biomaRt")
  library(biomaRt)}

if(!require(ggvenn)){
  BiocManager::install("ggvenn")
  library(ggvenn)}

if(!require(EnhancedVolcano)){
  BiocManager::install("EnhancedVolcano")
  library(EnhancedVolcano)}

if(!require(tidyr)){
  BiocManager::install("tidyr")
  library(tidyr)}

if(!require(readxl)){
  BiocManager::install("readxl")
  library(readxl)}

# SET DIRECTORY
DATA.DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(DATA.DIR)

# OPEN DATA_MARKERS FILE (THE CORRECTED COUNTS DATA)
data_markers <- read.table("data_markers.txt", sep = "\t")

# DOWNLOAD WIKIPATHWAYS DATABASE
downloadPathwayArchive(date="20220810", organism="Homo sapiens", format = "gmt")
wp <- list.files(pattern = "wikipathways")
wp <- clusterProfiler::read.gmt(wp)
colnames(wp)[c(1,2)] <- c("pathway", "entrezgene")
wp <- wp[wp$pathway %in% names(table(wp$pathway))[table(wp$pathway) >= 1],]
colnames(wp)[2] <- "entrezgene_id"
wp.count <- wp %>% group_by(pathway) %>% mutate(count = n())
wp.count <- wp.count[,c(1,3)]
wp.count <- unique(wp.count)

# CONVERT ENTREZ TO ENSEMBL GENE ID
grch38.gene <- useMart(biomart = "ensembl", 
                        dataset = "hsapiens_gene_ensembl")
geneNames <- getBM(attributes = c("hgnc_symbol",
                                   "ensembl_gene_id","entrezgene_id"), 
                    filters = "ensembl_gene_id", 
                    values = unique(rownames(data_markers)), 
                    mart = grch38.gene)

# check filtering > loosing 12k genes is a lot!

wp <- wp[wp$entrezgene_id %in% geneNames$entrezgene_id,]
geneNames$entrezgene_id <- as.character(geneNames$entrezgene_id)
wp <- full_join(wp, geneNames)
wp <- wp[!is.na(wp$pathway),]
wp <- merge(wp, wp.count, by = "pathway")
wp <- separate(data = wp, col = pathway, into = c("name", "version", "id", "species"), sep = "%")

# disease and size filtering for pathways
pathway_disease_filtering <- read_excel("pathway disease filtering.xlsx", sheet = "Diseases", col_names = FALSE)
colnames(pathway_disease_filtering) <- c("num","pathway")
pathway_disease_filtering <- separate(data = pathway_disease_filtering, col = pathway, into = c("name", "version", "id", "species"), sep = "%")

disease.pwy <- pathway_disease_filtering$id
wp <- subset(wp, !(id %in% disease.pwy))
wp<-subset(wp, count>10)

#=============================================================
# CREATE META-DATA FOR GROUPING OF SPECIES SAMPLES
groupLabels <- t(data_markers)
groupLabels <- as.data.frame(groupLabels)
groupLabels$sample <- rownames(groupLabels)
groupLabels$labels <- NA
ncol(groupLabels)
groupLabels <- groupLabels[,c(9670,9671)]
groupLabels$labels <- groupLabels$sample 
groupLabels$labels <- gsub("_.*", "", groupLabels$labels)

# LOOP OVER PATHWAYs AND CALCULATE STATISTIC SCORE PER COMPARISON PER PATHWAY GENELIST
#   AND ALSO AN EFFECT SIZE BASED ON MEAN 1 MINUS MEAN 2
# First make some empty list to put the results in later
mwVals <- list()
mwEffectVals <- list()

# Take one pathway from the total set of pathways
for (p in unique(wp$id)) {
  # Take the genes in that pathway
  pathGenes <- wp[wp$id==p,"ensembl_gene_id"]
  # Now subset the data_markers to only those genes
  test_markers <- data_markers[rownames(data_markers) %in% pathGenes,]
  
  # Split the test_markers based on the labels, so values grouped per species
  testData_human <- matrix(t(as.matrix(test_markers[1:5])), nrow = 1)
  testData_rabbit <- matrix(t(as.matrix(test_markers[6:10])), nrow = 1)
  testData_sheep <- matrix(t(as.matrix(test_markers[11:15])), nrow = 1)
  testData <- list (Human = testData_human, Rabbit = testData_rabbit, Sheep = testData_sheep)
  #testData <- split(as.vector(as.matrix(test_markers)), groupLabels[, "labels"])
  
  # FIRST compare HUMAN and RABBIT
  mw <- NA
  try ({
    mw <- wilcox.test(testData$Human, 
                      testData$Rabbit, na.rm = T)$p.value
  }, silent = TRUE)
  # Effect size
  eff <- NA
  try ({eff <- (mean(testData[["Human"]], na.rm = T) - mean(testData[["Rabbit"]], 
                                                            na.rm = T))
  }, silent = TRUE)
  # Send output to results list for Human vs Rabbit under name of the pathway
  mwVals$Human_Rabbit[[paste0(p)]] <- mw
  mwEffectVals$Human_Rabbit[[paste0(p)]] <- eff
  
  # SECOND compare HUMAN and SHEEP
  mw <- NA
  try ({
    mw <- wilcox.test(testData$Human, 
                      testData$Sheep, na.rm = T)$p.value
  }, silent = TRUE)
  # Effect size
  eff <- NA
  try ({eff <- (mean(testData[["Human"]], na.rm = T) - mean(testData[["Sheep"]], 
                                                            na.rm = T))
  }, silent = TRUE)
  # Send output to results list for Human vs Sheep under name of the pathway
  mwVals$Human_Sheep[[paste0(p)]] <- mw
  mwEffectVals$Human_Sheep[[paste0(p)]] <- eff
  
  # THIRD compare SHEEP and RABBIT
  mw <- NA
  try ({
    mw <- wilcox.test(testData$Sheep, 
                      testData$Rabbit, na.rm = T)$p.value
  }, silent = TRUE)
  # Effect size
  eff <- NA
  try ({eff <- (mean(testData[["Sheep"]], na.rm = T) - mean(testData[["Rabbit"]], 
                                                             na.rm = T))
  }, silent = TRUE)
  # Send output to results list for Rabbit vs Sheep under name of the pathway
  mwVals$Sheep_Rabbit[[paste0(p)]] <- mw
  mwEffectVals$Sheep_Rabbit[[paste0(p)]] <- eff
  
  # Now the loop restarts and goes on to the next pathway name, and calculates
  # the mean counts for those genes per column/sample
}

# COMBINE AND STRUCTURE THE RESULTS NICELY
resultsHR <- as.data.frame(as.matrix(mwVals$Human_Rabbit))
resultsHS <- as.data.frame(as.matrix(mwVals$Human_Sheep))
resultsSR <- as.data.frame(as.matrix(mwVals$Sheep_Rabbit))

colnames(resultsHR) <- "Human_Rabbit"
colnames(resultsHS) <- "Human_Sheep"
colnames(resultsSR) <- "Sheep_Rabbit"

resultsHR$pathway <- rownames(resultsHR)
resultsHS$pathway <- rownames(resultsHS)
resultsSR$pathway <- rownames(resultsSR)

resultsEffectHR <- as.data.frame(as.matrix(mwEffectVals$Human_Rabbit))
resultsEffectHS <- as.data.frame(as.matrix(mwEffectVals$Human_Sheep))
resultsEffectSR <- as.data.frame(as.matrix(mwEffectVals$Sheep_Rabbit))

colnames(resultsEffectHR) <- "Human_Rabbit_Effect"
colnames(resultsEffectHS) <- "Human_Sheep_Effect"
colnames(resultsEffectSR) <- "Sheep_Rabbit_Effect"

resultsEffectHR$pathway <- rownames(resultsEffectHR)
resultsEffectHS$pathway <- rownames(resultsEffectHS)
resultsEffectSR$pathway <- rownames(resultsEffectSR)

resultsTotal <- full_join(resultsHR,resultsHS)
resultsTotal <- full_join(resultsTotal, resultsSR)
resultsTotal <- full_join(resultsTotal, resultsEffectHR)
resultsTotal <- full_join(resultsTotal, resultsEffectHS)
resultsTotal <- full_join(resultsTotal, resultsEffectSR)

resultsTotal <- resultsTotal[,c(2,1,5,3,6,4,7)]
resultsTotal$Human_Rabbit <- as.numeric(resultsTotal$Human_Rabbit)
resultsTotal$Human_Sheep <- as.numeric(resultsTotal$Human_Sheep)
resultsTotal$Sheep_Rabbit <- as.numeric(resultsTotal$Sheep_Rabbit)
resultsTotal$Human_Rabbit_Effect <- as.numeric(resultsTotal$Human_Rabbit_Effect)
resultsTotal$Human_Sheep_Effect <- as.numeric(resultsTotal$Human_Sheep_Effect)
resultsTotal$Sheep_Rabbit_Effect <- as.numeric(resultsTotal$Sheep_Rabbit_Effect)
rownames(resultsTotal) <- resultsTotal$pathway

write.table(resultsTotal, "pathway_analysis.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#============================================================== Filter pathways with pval> 0.05 per each species comparison

HumanVsRabbitSorted<-subset(resultsTotal, Human_Rabbit<0.05)
HumanVsRabbitSorted<-subset(HumanVsRabbitSorted, abs(Human_Rabbit_Effect)>1)
#add step if HumanSheep is <0.05 then Human_Sheep_Effect must be >1?
HumanVsRabbitSorted<-HumanVsRabbitSorted%>%
  arrange(Human_Rabbit)

HumanVsSheepSorted<-subset(resultsTotal, Human_Sheep<0.05)
HumanVsSheepSorted<-subset(HumanVsSheepSorted, abs(Human_Sheep_Effect)>1)
HumanVsSheepSorted<-HumanVsSheepSorted%>%
  arrange(Human_Sheep)

SheepVsRabbitSorted<-subset(resultsTotal, Sheep_Rabbit<0.05)
SheepVsRabbitSorted<-subset(SheepVsRabbitSorted, abs(Sheep_Rabbit_Effect)>1)
#add step if HumanSheep is <0.05 then Human_Sheep_Effect must be >1?
SheepVsRabbitSorted<-SheepVsRabbitSorted%>%
  arrange(Sheep_Rabbit)

#Take note: In each comparison (i.e. Rabbit_Sheep) the effect is the mean of Rabbit minus mean of Sheep. This means that a negative value would mean that the Rabbit pathway is less expressed than the Sheep one.

pathways<-list(HumanSheep=HumanVsSheepSorted$pathway, HumanRabbit=HumanVsRabbitSorted$pathway, SheepRabbit=SheepVsRabbitSorted$pathway)
pathways2<-list(HumanRabbit=HumanVsRabbitSorted$pathway, SheepRabbit=SheepVsRabbitSorted$pathway)

ggvenn(pathways, stroke_size = 1, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),set_name_size = 10)
ggvenn(pathways2, stroke_size = 1, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),set_name_size = 10)


#MISTAKE HERE, WE ARE NOT SORTING BASED ON EFFECT IN SHEEP DUE TO HOW WE MERGE THE DATA
#add in highlight effect size sorting
HighlightDef <- subset(HumanVsRabbitSorted, pathway %in% SheepVsRabbitSorted$pathway)
HighlightDef <- merge(HighlightDef, unique(wp[,c(1,3)]), by.x="pathway", by.y ="id")


EnhancedVolcano(resultsTotal, lab = "", x = "Human_Rabbit_Effect", y = "Human_Rabbit", pCutoff = 0.05, FCcutoff = 1)
EnhancedVolcano(resultsTotal, lab = "", x = "Sheep_Rabbit_Effect", y = "Sheep_Rabbit", pCutoff = 0.05, FCcutoff = 1)
EnhancedVolcano(resultsTotal, lab = "", x = "Human_Sheep_Effect", y = "Human_Sheep", pCutoff = 0.05, FCcutoff = 1)

write.csv2(resultsTotal, "results_pathway_analysis.csv", row.names=TRUE, col.names = TRUE)
write.csv2(HighlightDef, "results_SheepVRabbit_HumanVRabbit.csv",row.names=TRUE, col.names = TRUE )



#to write table with Highlight pathways
write.table(HighlightDef, "HumanVRabbit_SheepVRabbit.txt", sep = "\t", quote = F, row.names = F, col.names = T)
