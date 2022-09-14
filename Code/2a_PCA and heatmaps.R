library(rstudioapi)
library(plyr)
library(splitstackshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(viridisLite)
library(DESeq2)
library(ggfortify)

# set working environment to source file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# read merged data file
DataAll <- read.delim("MergedDataFull.tsv")

Data<-DataAll%>% drop_na()
Data<-Data[,c(1,3:7,10:14,17:21)]

#Removing duplicate gene names
Data<-Data %>% distinct(ENSGID, .keep_all = TRUE)


#############################################################################

#Format table, Gene ID column should be removed as variable and introduced as row names and log transform data
Data2<-Data[,-1]
rownames(Data2) <- Data[,1]


Data3<-Data2[,c(1,2,3,4,5,11,12,13,14,15,6,7,8,9,10)]
colnames(Data3)<-c("Human_1", "Human_2", "Human_3", "Human_4", "Human_5", "Sheep_1", "Sheep_2", "Sheep_3", "Sheep_4", "Sheep_5", "Rabbit_1", "Rabbit_2", "Rabbit_3", "Rabbit_4", "Rabbit_5")



#Heatmap all
pheatmap(Data3, color =viridis(256,option="D"),cluster_rows = T ,cluster_cols = T, scale="none", show_rownames = F )

#PCA
PCA<-prcomp(t(Data3),scale=TRUE)

autoplot(PCA)

PCA.Data<-data.frame(Sample=rownames(PCA$x), X=PCA$x[,1],Y=PCA$x[,2])
ggplot(data=PCA.Data, aes(x=X, y=Y, color=Sample, label=Sample))+
  geom_point(size=5)+geom_text(size=4,nudge_x = 0.25, nudge_y = -0.08, check_overlap = T)+
  xlab("PCA1")+ylab("PCA2")+
  theme_light()

ggplot(data=PCA.Data, aes(x=X, y=Y, color=Sample, label=Sample, position="jitter"))+
  geom_jitter(width=3.3,height=-0.8, size=5)+
  xlab("PCA1: 41.64% variance")+ylab("PCA2: 27.6% variance")+
  theme_bw(base_size = 20)

#Correlation heatmap
correlation<-cor(Data3)
pheatmap(correlation, color=viridis(256,option="D"), cluster_rows=FALSE, cluster_cols = FALSE, fontsize = 15)




write.table(Data2,"data_markers.txt", col.names = T, row.names = T, quote = F, sep = "\t")

############################################################################################### Selection of layer specific markers

#(ALCAM, TJP1, NCAM1, SLC4A4, ATP1A1, KERA, CD34, LUM, WNT7a, PAX6, KRT5, ITGAM, CCR7, CD19)
#Absent epithelium markers (not detected): KRT12 (not in human, rabbit 0,0,1,38,0, sheep 0,0,0,0,1), GJB6, KRT14 (absent rabbit and sheep, human reads 0,0,0,0,4 -> lacking homology)
#Absent stromal markers (not detected in rabbit): PTGDS, ALDH3A1
#Absent immune markers: CD16a, CD3, CD11b
#Absent blood markers CCL21, CLDN5, GYPA
#Absent Endo marker due to homology PRDX6 (present in all samples)


Markers<-Data3[c("ENSG00000170017","ENSG00000104067","ENSG00000149294","ENSG00000080493", "ENSG00000163399", "ENSG00000139330", "ENSG00000174059", "ENSG00000139329", "ENSG00000154764", "ENSG00000007372", "ENSG00000169896", "ENSG00000126353", "ENSG00000177455"),]

rownames(Markers)<-c("ALCAM", "TJP1", "NCAM1", "SLC4A4", "ATP1A1", "KERA", "CD34", "LUM","WNT7A", "PAX6","ITGAM", "CCR7", "CD19")


pheatmap(Markers, color =viridis(256,option="D"),cluster_rows = F,cluster_cols = F, scale="none", show_rownames = T, na_col = "grey", fontsize = 15 )





Data2$HumanAverage<-apply(Data2[,1:5], 1, mean)
Data2$RabbitAverage<-apply(Data2[,6:10], 1, mean)
Data2$SheepAverage<-apply(Data2[,11:15], 1, mean)


write.table(Data2[,16:18], "DataWithAverages.txt", col.names = T, row.names = T, quote = F, sep="\t")

