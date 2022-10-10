# Homology mapping
# Authors: PereCatQ and mkutmon
# Description: Mapping of all measured human genes to homologs in rabbit and sheep using orthogene R-package

#-----------------------------
# Setup and required libraries 
#-----------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", update=TRUE)
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi", update=TRUE)
if(!"orthogene" %in% installed.packages()) BiocManager::install("orthogene", update=TRUE)
if(!"splitstackshape" %in% installed.packages()) BiocManager::install("splitstackshape", update=TRUE)
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr", update=TRUE)

library(rstudioapi)
library(orthogene)
library(splitstackshape)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#---------------------------------------------------------------------
# Get list of human identifiers in dataset (protein coding genes only)
#---------------------------------------------------------------------

human <- read.table("count-data/human/MER-PC-b002_total.ReadCounts.tsv", header=TRUE, sep="\t")
human.split <- splitstackshape::cSplit(human, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]
human.filt <- human.split[human.split$X_3 == "ProteinCoding",]
colnames(human.filt)[1:3] <- c("GeneID", "GeneName", "Type")

human.ids <- human.filt[,c(1,2)]
rm(human,human.split,human.filt)

#--------------------------------------------
# Map human genes to sheep (gprofiler method)
#--------------------------------------------
method <- "gprofiler"  
mapping.sheep <- orthogene::convert_orthologs(gene_df = human.ids$GeneID,
                                        gene_input = "Human.GeneID", 
                                        gene_output = "columns", 
                                        input_species = "human",
                                        output_species = "oarambouillet",
                                        non121_strategy = "kbs",
                                        method = method) 
mapping.sheep <- mapping.sheep[,c(2,3)]
colnames(mapping.sheep) <- c("Human.GeneID","Sheep.GeneName")

mapping.sheep.ids <- orthogene::map_genes(genes = mapping.sheep$Sheep.GeneName, species = "oarambouillet")

mapping.sheep <- merge(mapping.sheep, mapping.sheep.ids[,c(2,4)], by.x="Sheep.GeneName", by.y = "input",  all.x=TRUE)
mapping.sheep <- merge(mapping.sheep, human.ids, by.x="Human.GeneID", by.y = "GeneID",  all.x=TRUE)
mapping.sheep <- mapping.sheep[,c(1,4,3,2)]
colnames(mapping.sheep) <- c("Human.GeneID","Human.Name","Sheep.GeneID","Sheep.Name")
mapping <- mapping.sheep[startsWith(mapping.sheep$Sheep.GeneID,"ENSOAR"),]
mapping <- dplyr::distinct(mapping)

rm(mapping.sheep,mapping.sheep.ids)

#--------------------------------------------
# Map human genes to rabbit (gprofiler method)
#--------------------------------------------

# convert to rabbit
mapping.rabbit <- orthogene::convert_orthologs(gene_df = human.ids$GeneID,
                                              gene_input = "Human.GeneID", 
                                              gene_output = "columns", 
                                              input_species = "human",
                                              output_species = "rabbit",
                                              non121_strategy = "kbs",
                                              method = method) 

mapping.rabbit <- mapping.rabbit[,c(2,3)]
colnames(mapping.rabbit) <- c("Human.GeneID","Rabbit.GeneName")

mapping.rabbit.ids <- orthogene::map_genes(genes = mapping.rabbit$Rabbit.GeneName, species = "rabbit")

mapping.rabbit <- merge(mapping.rabbit, mapping.rabbit.ids[,c(2,4)], by.x="Rabbit.GeneName", by.y = "input",  all.x=TRUE)
mapping.rabbit <- mapping.rabbit[,c(2,3,1)]
colnames(mapping.rabbit) <- c("Human.GeneID","Rabbit.GeneID","Rabbit.Name")
mapping.rabbit <- mapping.rabbit[startsWith(mapping.rabbit$Rabbit.GeneID,"ENSOCUG"),]
mapping.rabbit <- dplyr::distinct(mapping.rabbit)
mapping <- merge(mapping, mapping.rabbit, by.x="Human.GeneID", by.y = "Human.GeneID")


rm(mapping.rabbit, mapping.rabbit.ids)

#--------------------------------------------
# Export mapping file
#--------------------------------------------

write.table(mapping, file="data/homology-mapping.tsv", quote=FALSE, sep="\t", row.names=FALSE)
