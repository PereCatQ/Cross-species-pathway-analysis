# Pathway activity analysis
# Authors: PereCatQ and mkutmon
# Description: Assessing the activity of the pathways and comparing the activity in the different species

#-----------------------------
# Setup and required libraries 
#-----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", update=FALSE)
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi", update=FALSE)
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways", update=FALSE)
if(!"biomaRt" %in% installed.packages()) BiocManager::install("biomaRt", update=FALSE)
if(!"readxl" %in% installed.packages()) BiocManager::install("readxl", update=FALSE)
if(!"EnhancedVolcano" %in% installed.packages()) BiocManager::install("EnhancedVolcano", update=FALSE)
if(!"ggvenn" %in% installed.packages()) BiocManager::install("ggvenn", update=FALSE)
if(!"tidyr" %in% installed.packages()) BiocManager::install("tidyr", update=FALSE)
if(!"data.table" %in% installed.packages()) BiocManager::install("data.table", update=FALSE)
if(!"pheatmap" %in% installed.packages()) BiocManager::install("pheatmap", update=FALSE)

library(rstudioapi)
library(rWikiPathways)
library(readxl)
library(biomaRt)
library(ggvenn)
library(EnhancedVolcano)
library(tidyr)
library(data.table)
library(pheatmap)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#------------------------------------
# Import combined normalized data
#------------------------------------
data <- read.table("data/normalized-data.tsv", sep = "\t", header=TRUE)

# identifier mapping 
biomart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chro <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","MT","X","Y")
geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id","entrezgene_id"), 
                    filters = c("ensembl_gene_id", "chromosome_name"), 
                    values = list(ensembl_gene_id=unique(data$ENSGID),chromosome_name=chro), 
                    mart = biomart)
geneNames <- geneNames[!is.na(geneNames$entrezgene_id),]

# 9757 identifiers in homology mapped data - 9717 with mapping to NCBI Gene - make sure you don't loose too many genes due to ID mapping
length(unique(geneNames$ensembl_gene_id))

#----------------------------------------------------
# Retrieve pathway collection from WikiPathways
#----------------------------------------------------
gmt <- rWikiPathways::downloadPathwayArchive(date="20220810", organism="Homo sapiens", format = "gmt")
wp <- clusterProfiler::read.gmt(gmt)
colnames(wp)[c(1,2)] <- c("pathway", "entrezgene")

# number of genes in pathway
wp.count <- wp %>% group_by(pathway) %>% mutate(count = n())
wp.count <- wp.count[,c(1,3)]
wp.count <- unique(wp.count)

wp.filt <- wp[wp$entrezgene %in% geneNames$entrezgene_id,]
geneNames$entrezgene_id <- as.character(geneNames$entrezgene_id)
wp.filt <- merge(wp.filt, geneNames, by.x="entrezgene", by.y="entrezgene_id", all.x=TRUE)
wp.filt <- merge(wp.filt, wp.count, by = "pathway")
wp.filt <- tidyr::separate(data = wp.filt, col = pathway, into = c("name", "version", "id", "species"), sep = "%")

#----------------------------------------------------
# Filter pathways (disease & size)
#----------------------------------------------------

pathway_disease_filtering <- readxl::read_excel("data/disease-pathways.xlsx", sheet = "Diseases", col_names = TRUE)
pathway_disease_filtering <- separate(data = pathway_disease_filtering, col = Pathway, into = c("name", "version", "id", "species"), sep = "%")

disease.pwy <- pathway_disease_filtering$id
wp.filt <- subset(wp.filt, !(id %in% disease.pwy))
wp.filt <- subset(wp.filt, count>10)

#----------------------------------------------------
# Pathway activity comparison between species
#----------------------------------------------------

# groups
samples <- as.data.frame(colnames(data[,3:17]))
samples$labels <- c("Human","Human","Human","Human","Human","Rabbit","Rabbit","Rabbit","Rabbit","Rabbit","Sheep","Sheep","Sheep","Sheep","Sheep")
colnames(samples) <- c("sample","labels")
rownames(samples) <- samples$sample

# pathway activity statisical comprison (wilcox test)
# LOOP OVER PATHWAYs AND CALCULATE STATISTIC SCORE PER COMPARISON PER PATHWAY GENELIST
#   AND ALSO AN EFFECT SIZE BASED ON MEAN 1 MINUS MEAN 2
# First make some empty list to put the results in later
mwVals <- list()
mwEffectVals <- list()

# Take one pathway from the total set of pathways
for (p in unique(wp.filt$id)) {
  # Take the genes in that pathway
  pathGenes <- wp.filt[wp.filt$id==p,"ensembl_gene_id"]
  # Now subset the data_markers to only those genes
  data.selected <- data[data$ENSGID %in% pathGenes,3:17]
  
  # Split the data based on the labels, so values grouped per species
  testData_human <- matrix(t(as.matrix(data.selected[1:5])), nrow = 1)
  testData_rabbit <- matrix(t(as.matrix(data.selected[6:10])), nrow = 1)
  testData_sheep <- matrix(t(as.matrix(data.selected[11:15])), nrow = 1)
  testData <- list (Human = testData_human, Rabbit = testData_rabbit, Sheep = testData_sheep)

  # FIRST compare HUMAN and RABBIT
  mw <- wilcox.test(testData$Human, testData$Rabbit, na.rm = T)$p.value # significance
  eff <- (mean(testData[["Human"]], na.rm = T) - mean(testData[["Rabbit"]], na.rm = T)) # effect size
  mwVals$Human_Rabbit[[paste0(p)]] <- mw
  mwEffectVals$Human_Rabbit[[paste0(p)]] <- eff
  
  # SECOND compare HUMAN and SHEEP
  mw <- wilcox.test(testData$Human, testData$Sheep, na.rm = T)$p.value # significance
  eff <- (mean(testData[["Human"]], na.rm = T) - mean(testData[["Sheep"]], na.rm = T)) # effect size
  mwVals$Human_Sheep[[paste0(p)]] <- mw
  mwEffectVals$Human_Sheep[[paste0(p)]] <- eff
  
  # THIRD compare SHEEP and RABBIT
  mw <- wilcox.test(testData$Sheep, testData$Rabbit, na.rm = T)$p.value # significance
  eff <- (mean(testData[["Sheep"]], na.rm = T) - mean(testData[["Rabbit"]], na.rm = T)) # effect size
  mwVals$Sheep_Rabbit[[paste0(p)]] <- mw
  mwEffectVals$Sheep_Rabbit[[paste0(p)]] <- eff
}

# COMBINE AND STRUCTURE THE RESULTS NICELY
res.pval <- as.data.frame(t(do.call(rbind, mwVals)))
res.effect <- as.data.frame(t(do.call(rbind, mwEffectVals)))

colnames(res.pval) <- paste("Pval", colnames(res.pval), sep = "_")
colnames(res.effect) <- paste("Effect", colnames(res.effect), sep = "_")

results.total <- cbind(res.pval, res.effect)
results.total$Pathway <- rownames(results.total)
results.total <- results.total[,c(7,4,1,6,3,5,2)]
results.total[,2] <- as.numeric(results.total[,2])
results.total[,3] <- as.numeric(results.total[,3])
results.total[,4] <- as.numeric(results.total[,4])
results.total[,5] <- as.numeric(results.total[,5])
results.total[,6] <- as.numeric(results.total[,6])
results.total[,7] <- as.numeric(results.total[,7])

results.total <- merge(results.total, unique(wp.filt[,c(1,3)]), by.x="Pathway", by.y="id", all.x=TRUE) 

write.table(results.total, "data/pathway-analysis.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#----------------------------------------------------
# Pathway overlap analysis
#----------------------------------------------------

png('figures/volcano_HR.png', width = 2300, height = 2000, res = 300)
EnhancedVolcano(results.total, lab = "", x = "Effect_Human_Rabbit", y = "Pval_Human_Rabbit", pCutoff = 0.05, FCcutoff = 1, xlab = "effect size", ylab = bquote(~-log[10] ~ "p-value"),legendLabels = c("NS", "effect size", "p-value", "p-value and effect size"), caption = "", subtitle="Volcano plot pathway activity effect size (Human vs Rabbit)", title="", titleLabSize = 0, captionLabSize = 0, xlim=c(-5,5), ylim= c(0,25))
dev.off()

png('figures/volcano_SR.png', width = 2300, height = 2000, res = 300)
EnhancedVolcano(results.total, lab = "", x = "Effect_Sheep_Rabbit", y = "Pval_Sheep_Rabbit", pCutoff = 0.05, FCcutoff = 1, xlab = "effect size", ylab = bquote(~-log[10] ~ "p-value"),legendLabels = c("NS", "effect size", "p-value", "p-value and effect size"), caption = "", subtitle="Volcano plot pathway activity effect size (Sheep vs Rabbit)", title="", titleLabSize = 0, captionLabSize = 0, xlim=c(-5,5), ylim= c(0,25))
dev.off()

png('figures/volcano_HS.png', width = 2300, height = 2000, res=300)
EnhancedVolcano(results.total, lab = "", x = "Effect_Human_Sheep", y = "Pval_Human_Sheep", pCutoff = 0.05, FCcutoff = 1, xlab = "effect size", ylab = bquote(~-log[10] ~ "p-value"),legendLabels = c("NS", "effect size", "p-value", "p-value and effect size"), caption = "", subtitle="Volcano plot pathway activity effect size (Human vs Sheep)", title="", titleLabSize = 0, captionLabSize = 0, xlim=c(-5,5), ylim= c(0,25))
dev.off()

sign.pwy.Human_Rabbit <- subset(results.total, Pval_Human_Rabbit < 0.05 & abs(Effect_Human_Rabbit) > 1)
sign.pwy.Sheep_Rabbit <- subset(results.total, Pval_Sheep_Rabbit < 0.05 & abs(Effect_Sheep_Rabbit) > 1)
sign.pwy.Human_Sheep <- subset(results.total, Pval_Human_Sheep < 0.05 & abs(Effect_Human_Sheep) > 1)

pathways <- list(HumanSheep=sign.pwy.Human_Sheep$Pathway, HumanRabbit=sign.pwy.Human_Rabbit$Pathway, SheepRabbit=sign.pwy.Sheep_Rabbit$Pathway)
ggvenn(pathways, stroke_size = 1, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),set_name_size = 5, show_percentage = FALSE)

selected.pwys <- subset(sign.pwy.Human_Rabbit, Pathway %in% sign.pwy.Sheep_Rabbit$Pathway)
rownames(selected.pwys) <- selected.pwys$name

data.heatmap <- selected.pwys[,c(2,4)]
colnames(data.heatmap) <- c("Human vs Rabbit", "Sheep vs Rabbit")

png('figures/heatmap-pathways.png', width = 2500, height = 3800, res=300)
pheatmap(data.heatmap, color = viridisLite::viridis(256,option="D"),cluster_rows = T ,cluster_cols = F, scale="none", show_rownames = T, angle_col = 315,fontsize_col = 16, fontsize_row = 10)
dev.off()

write.table(selected.pwys, "data/selected-pathways.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
