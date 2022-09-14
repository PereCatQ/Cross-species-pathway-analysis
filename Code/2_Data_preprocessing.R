if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"plyr" %in% installed.packages()) BiocManager::install("plyr")
if(!"splitstackshape" %in% installed.packages()) BiocManager::install("splitstackshape")
if(!"limma" %in% installed.packages()) BiocManager::install("limma")


library(rstudioapi)
library(plyr)
library(splitstackshape)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridisLite)
library(limma)
library(ggplot2)
library(tibble)


# set working environment to source file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read homology mapping file
mapping <- read.delim("homology-mapping/homology-mapping.tsv")

# read data files
human <- read.table("data-2022-02-21/human/MER-PC-b002_total.ReadCounts.tsv", header=TRUE, sep="\t")
human.split <- cSplit(human, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]
human.filt <- human.split[human.split$X_3 == "ProteinCoding",]
human.filt[,c(4:8)]

# read rabbit files
rabbit1 <- read.table("data-2022-02-21/rabbit/MER-PC-b004_total.ReadCounts.tsv", header=TRUE, sep="\t")
rabbit2 <- read.table("data-2022-02-21/rabbit/MER-PC-b005_total.ReadCounts.tsv", header=TRUE, sep="\t")

#Rabbit data full join
rabbit.full<-full_join(rabbit1, rabbit2, by = "X")
#Heatmap to check for missing genes accross rabbit
rabbit.full2<-rabbit.full[,c(2:6)]
rabbit.full2<-log2(rabbit.full2+1)
pheatmap(rabbit.full2, color =viridis(256,option="D"),cluster_rows = F,cluster_cols = F, scale="none", show_rownames = F )
#Change NA to 0 in rabbit full join
rabbit.full[is.na(rabbit.full)] <- 0

rabbit.split <- cSplit(rabbit.full, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]
rabbit.filt <- rabbit.split[rabbit.split$X_3 == "ProteinCoding",]
rabbit.filt[,c(4:8)]

# read sheep files
sheep <- read.table("data-2022-02-21/sheep/MER-PC-b003_total.ReadCounts.tsv", header=TRUE, sep="\t")
sheep.split <- cSplit(sheep, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]

sheep.filt <- sheep.split[sheep.split$X_03 == "ProteinCoding",]
sheep.filt[,c(4:8)]

rm(human.split, rabbit.split,rabbit1,rabbit2, sheep.split, rabbit.full2, rabbit.full)


#################################THIS STEP WILL ADD ROW MEANS, AND THE COMBINED DATA STEP WILL GIVE AN ERROR IF NOT REMOVED, ONLY RUN IN CASE YOU WANT TO SEE GENE DISTRIBUTIONS FOR EACH SPECIES

# Visualize distributions
sheep.filt[,c(4:8)] <- log2(sheep.filt[,c(4:8)]+1)
human.filt[,c(4:8)] <- log2(human.filt[,c(4:8)]+1)
rabbit.filt[,c(4:8)] <- log2(rabbit.filt[,c(4:8)]+1)

human.dens <- apply(human.filt[,4:8], 2, density)
plot(NA, xlim=range(sapply(human.dens, "[", "x")), ylim=range(sapply(human.dens, "[", "y")))
mapply(lines, human.dens, col=1:length(human.dens))
legend("topright", legend=names(human.dens), fill=1:length(human.dens))
human.filt$rowMean <- rowMeans(human.filt[,4:8])

rabbit.dens <- apply(rabbit.filt[,4:8], 2, density)
plot(NA, xlim=range(sapply(rabbit.dens, "[", "x")), ylim=range(sapply(rabbit.dens, "[", "y")))
mapply(lines, rabbit.dens, col=1:length(rabbit.dens))
legend("topright", legend=names(rabbit.dens), fill=1:length(rabbit.dens))
rabbit.filt$rowMean <- rowMeans(rabbit.filt[,4:8])

sheep.dens <- apply(sheep.filt[,4:8], 2, density)
plot(NA, xlim=range(sapply(sheep.dens, "[", "x")), ylim=range(sapply(sheep.dens, "[", "y")))
mapply(lines, sheep.dens, col=1:length(sheep.dens))
legend("topright", legend=names(sheep.dens), fill=1:length(sheep.dens))
sheep.filt$rowMean <- rowMeans(sheep.filt[,4:8])

rm(human.dens, rabbit.dens, sheep.dens)


#========================================================================================================================================
# Combine data

# map rabbit to human
rabbit.hs <- unique(merge(rabbit.filt, mapping[,c(1,5)], by.x = "X_1", by.y = "Rabbit.GeneID", all.y = TRUE))

# map sheep to human
sheep.hs<-unique(merge(sheep.filt, mapping[,c(1,4)], by.x= "X_01", by.y = "Sheep.GeneName", all.y=TRUE))


#==============================================================================================================================================================
# merge all datasets
# TODO: decide what to do with genes not mapped to human (currently they are removed)
combined.data <- merge(human.filt, rabbit.hs, by.x = "X_1", by.y = "Human.GeneID", all.y = TRUE)
colnames(combined.data) <- c("ENSGID","GeneName","Type","Human_1","Human_2","Human_3","Human_4","Human_5","ENSOCUGID","OcGeneName","OcType","Rabbit_1","Rabbit_2","Rabbit_3","Rabbit_4","Rabbit_5")

#Adding sheep data
combined.data <- merge(combined.data, sheep.hs, by.x = "ENSGID", by.y = "Human.GeneID", all.y = TRUE)
combined.data<- merge(combined.data, mapping [,c(1,3)], by.x="ENSGID", by.y="Human.GeneID", all.x=TRUE)
combined.data<-combined.data[,c(1,2, 4:10, 12:16, 25,17,20:24)]


combined.data<-combined.data%>% drop_na()
combined.data[,c(3:7, 10:14, 17:21)] <- log2(combined.data[,c(3:7, 10:14, 17:21)]+1)
combined.data[,c(3:7, 10:14, 17:21)] <- normalizeQuantiles(combined.data[,c(3:7, 10:14, 17:21)], ties=TRUE)


write.table(combined.data, file="MergedDataFull.tsv", quote=FALSE, sep="\t", col.names = TRUE, row.names = FALSE)


#======================================================== Visualize distributions after homology mapping THIS STEP WILL ADD ROW MEANS, AND THE COMBINED DATA STEP WILL GIVE AN ERROR IF NOT REMOVED, ONLY RUN IN CASE YOU WANT TO SEE GENE DISTRIBUTIONS AFTER HOMOLOGY MAPPING

# Visualize distributions
human.dens <- apply(combined.data[,3:7], 2, density)
plot(NA, xlim=range(sapply(human.dens, "[", "x")), ylim=range(c(0, 0.165)))
mapply(lines, human.dens, col=1:length(human.dens))
legend("topright", legend=names(human.dens), fill=1:length(human.dens))
human.filt$rowMean <- rowMeans(combined.data[,3:7])


rabbit.dens <- apply(combined.data[,10:14], 2, density)
plot(NA, xlim=range(sapply(rabbit.dens, "[", "x")), ylim=range(c(0,0.165)))
mapply(lines, rabbit.dens, col=1:length(rabbit.dens))
legend("topright", legend=names(rabbit.dens), fill=1:length(rabbit.dens))
rabbit.filt$rowMean <- rowMeans(rabbit.filt[,4:8])

sheep.dens <- apply(combined.data[,17:21], 2, density)
plot(NA, xlim=range(sapply(sheep.dens, "[", "x")), ylim=range(c(0,0.165)))
mapply(lines, sheep.dens, col=1:length(sheep.dens))
legend("topright", legend=names(sheep.dens), fill=1:length(sheep.dens))
sheep.filt$rowMean <- rowMeans(sheep.filt[,4:8])
