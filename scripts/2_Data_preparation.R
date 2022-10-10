# Data preparation and exploration
# Authors: PereCatQ and mkutmon
# Description: Merging of data from three species, data normalization, data exploration

#-----------------------------
# Setup and required libraries 
#-----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", update=FALSE)
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi", update=FALSE)
if(!"splitstackshape" %in% installed.packages()) BiocManager::install("splitstackshape", update=FALSE)
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr", update=FALSE)
if(!"pheatmap" %in% installed.packages()) BiocManager::install("pheatmap", update=FALSE)
if(!"viridisLite" %in% installed.packages()) BiocManager::install("viridisLite", update=FALSE)
if(!"tidyr" %in% installed.packages()) BiocManager::install("tidyr", update=FALSE)
if(!"limma" %in% installed.packages()) BiocManager::install("limma", update=FALSE)
if(!"reshape2" %in% installed.packages()) BiocManager::install("reshape2", update=FALSE)
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2", update=FALSE)

library(rstudioapi)
library(splitstackshape)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridisLite)
library(limma)
library(ggplot2)
library(reshape2)

# set working environment to source file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#----------------------------------------
# Read homology mapping file 
#----------------------------------------

# read homology mapping file (make sure script 1_Homology_mapping.R has been run)
mapping <- read.delim("data/homology-mapping.tsv")

#----------------------------------------
# Prepare human data 
#----------------------------------------

# read human files
human <- read.table("count-data/human/MER-PC-b002_total.ReadCounts.tsv", header=TRUE, sep="\t")
human.split <- splitstackshape::cSplit(human, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]
human.filt <- human.split[human.split$X_3 == "ProteinCoding",]
human.filt <- human.filt[,c(1,2,4,5,6,7,8)]

#----------------------------------------
# Prepare rabbit data
#----------------------------------------

# read rabbit files
rabbit1 <- read.table("count-data/rabbit/MER-PC-b004_total.ReadCounts.tsv", header=TRUE, sep="\t")
rabbit2 <- read.table("count-data/rabbit/MER-PC-b005_total.ReadCounts.tsv", header=TRUE, sep="\t")

rabbit.full <- dplyr::full_join(rabbit1, rabbit2, by = "X")

# Heatmap to check for missing genes accross rabbit (detected in one batch but not the other)
rabbit.full2<-rabbit.full[,c(2:6)]
rabbit.full2<-log2(rabbit.full2+1)
pheatmap(rabbit.full2, color = viridisLite::viridis(256,option="D"),cluster_rows = F,cluster_cols = F, scale="none", show_rownames = F )

# Change NA to 0 in rabbit full join (change from NA to 0 = not expressed)
rabbit.full[is.na(rabbit.full)] <- 0

rabbit.split <- splitstackshape::cSplit(rabbit.full, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]
rabbit.filt <- rabbit.split[rabbit.split$X_3 == "ProteinCoding",]
rabbit.filt <- rabbit.filt[,c(1,2,4,5,6,7,8)]

#----------------------------------------
# Prepare rabbit data
#----------------------------------------

sheep <- read.table("count-data/sheep/MER-PC-b003_total.ReadCounts.tsv", header=TRUE, sep="\t")
sheep.split <- splitstackshape::cSplit(sheep, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]

sheep.filt <- sheep.split[sheep.split$X_03 == "ProteinCoding",]
sheep.filt <- sheep.filt[,c(1,2,4,5,6,7,8)]

rm(human.split, rabbit.split,rabbit1,rabbit2, sheep.split, rabbit.full2, rabbit.full)

#----------------------------------------
# Combine data of the three species
#----------------------------------------

# map rabbit to human
human.rabbit <- unique(merge(rabbit.filt, mapping[,c(1,5)], by.x = "X_1", by.y = "Rabbit.GeneID", all.y = TRUE))

# map sheep to human
human.sheep <- unique(merge(sheep.filt, mapping[,c(1,4)], by.x= "X_01", by.y = "Sheep.Name", all.y=TRUE))

# genes not mapped to human are removed 
combined.data <- merge(human.filt, human.rabbit, by.x = "X_1", by.y = "Human.GeneID", all.y = TRUE)
colnames(combined.data) <- c("ENSGID","GeneName","Human_1","Human_2","Human_3","Human_4","Human_5","ENSOCUGID","OcGeneName","Rabbit_1","Rabbit_2","Rabbit_3","Rabbit_4","Rabbit_5")

combined.data <- merge(combined.data, human.sheep, by.x = "ENSGID", by.y = "Human.GeneID", all.y = TRUE)
combined.data<- merge(combined.data, mapping [,c(1,3)], by.x="ENSGID", by.y="Human.GeneID", all.x=TRUE)
combined.data.final <- combined.data[,c(1:7,10:14,17:21)]

# genes to measured in all species are dropped
combined.data.final <- combined.data.final%>% tidyr::drop_na()

#----------------------------------------
# Data normalization (quantile normalization)
#----------------------------------------

combined.data.final[,c(3:17)] <- log2(combined.data.final[,c(3:17)]+1)
combined.data.norm <- combined.data.final
combined.data.norm[,c(3:17)] <- limma::normalizeQuantiles(combined.data.norm[,c(3:17)], ties=TRUE)

write.table(combined.data.norm, file="data/normalized-data.tsv", quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

#--------------------------------------------------------------------------
# Exploration: visualization of data before and after normalization
#--------------------------------------------------------------------------

# before normalization
# =========================
# density
data.dens <- apply(combined.data.final[,3:17], 2, density)
plot(NA, xlim=range(sapply(data.dens, "[", "x")), ylim=range(sapply(data.dens, "[", "y")))
mapply(lines, data.dens, col=1:length(data.dens))
legend("topright", legend=names(data.dens), fill=1:length(data.dens), cex = 0.7, title="Density plot (before norm)")

# boxplot
melt.dat <-data.frame(reshape2::melt(combined.data.final[,c(1,3:17)], id.vars = "ENSGID"))
ggplot(melt.dat, aes(x = as.factor(variable), y = value)) + geom_boxplot(aes(fill = variable), position = position_dodge(0.9)) + theme_bw()

# heatmap
pheatmap(combined.data.final[,c(3:17)], color = viridisLite::viridis(256,option="D"),cluster_rows = T ,cluster_cols = T, scale="none", show_rownames = F )

# PCA
PCA<-prcomp(t(combined.data.final[,c(3:17)]),scale=TRUE)
PCA.Data<-data.frame(Sample=rownames(PCA$x), X=PCA$x[,1],Y=PCA$x[,2])
ggplot(data=PCA.Data, aes(x=X, y=Y, color=Sample, label=Sample))+
  geom_point(size=5)+geom_text(size=4,nudge_x = 0.25, nudge_y = -0.08, check_overlap = T)+
  xlab("PCA1")+ylab("PCA2")+
  theme_light()


# after normalization
# =========================
# density
data.norm.dens <- apply(combined.data.norm[,3:17], 2, density)
plot(NA, xlim=range(sapply(data.norm.dens, "[", "x")), ylim=range(sapply(data.norm.dens, "[", "y")))
mapply(lines, data.norm.dens, col=1:length(data.norm.dens))
legend("topright", legend=names(data.norm.dens), fill=1:length(data.norm.dens), cex = 0.7, title="Density plot (after norm)")

# boxplot
melt.dat.norm <-data.frame(reshape2::melt(combined.data.norm[,c(1,3:17)], id.vars = "ENSGID"))
ggplot(melt.dat.norm, aes(x = as.factor(variable), y = value)) + geom_boxplot(aes(fill = variable), position = position_dodge(0.9)) + theme_bw()

# heatmap
pheatmap(combined.data.norm[,c(3:17)], color = viridisLite::viridis(256,option="D"),cluster_rows = T ,cluster_cols = T, scale="none", show_rownames = F )

# PCA
pca.norm <- prcomp(t(combined.data.norm[,c(3:17)]),scale=TRUE)
pca.norm.data <- data.frame(Sample=rownames(pca.norm$x), X=pca.norm$x[,1],Y=pca.norm$x[,2])
ggplot(data=pca.norm.data, aes(x=X, y=Y, color=Sample, label=Sample))+
  geom_point(size=5)+geom_text(size=4,nudge_x = 0.25, nudge_y = -0.08, check_overlap = T)+
  xlab("PCA1")+ylab("PCA2")+
  theme_light()

#--------------------------------------------------------------------------
# Check: visualization of known endothelial, epithelium, stromal, immune and blood marker genes
#--------------------------------------------------------------------------
#(ALCAM, TJP1, NCAM1, SLC4A4, ATP1A1, KERA, CD34, LUM, WNT7a, PAX6, KRT5, ITGAM, CCR7, CD19)
#Absent epithelium markers (not detected): KRT12 (not in human, rabbit 0,0,1,38,0, sheep 0,0,0,0,1), GJB6, KRT14 (absent rabbit and sheep, human reads 0,0,0,0,4 -> lacking homology)
#Absent stromal markers (not detected in rabbit): PTGDS, ALDH3A1
#Absent immune markers: CD16a, CD3, CD11b
#Absent blood markers CCL21, CLDN5, GYPA
#Absent Endo marker due to homology PRDX6 (present in all samples)

markers <- combined.data.norm[c("ENSG00000170017","ENSG00000104067","ENSG00000149294","ENSG00000080493", "ENSG00000163399", "ENSG00000139330", "ENSG00000174059", "ENSG00000139329", "ENSG00000154764", "ENSG00000007372", "ENSG00000169896", "ENSG00000126353", "ENSG00000177455"),c(3:17)]
rownames(markers)<-c("ALCAM", "TJP1", "NCAM1", "SLC4A4", "ATP1A1", "KERA", "CD34", "LUM","WNT7A", "PAX6","ITGAM", "CCR7", "CD19")
pheatmap(markers, color = viridisLite::viridis(256,option="D"),cluster_rows = F,cluster_cols = F, scale="none", show_rownames = T, na_col = "grey", fontsize = 15 )
