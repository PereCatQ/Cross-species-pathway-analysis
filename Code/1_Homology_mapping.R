if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"plyr" %in% installed.packages()) BiocManager::install("plyr")
if(!"splitstackshape" %in% installed.packages()) BiocManager::install("splitstackshape")
if(!"orthogene" %in% installed.packages()) BiocManager::install("orthogene")
if(!"org.Mm.eg.db" %in% installed.packages()) BiocManager::install("org.Mm.eg.db")
if(!"biomaRt" %in% installed.packages()) BiocManager::install("biomaRt")

library(biomaRt)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(rstudioapi)
library(plyr)
library(dplyr)

library(splitstackshape)
library(orthogene)
library(clusterProfiler)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# set up identifier mapping
ensembl <- useEnsembl(biomart = "genes")
ensembl.sheep <- useDataset(dataset = "oarambouillet_gene_ensembl", mart = ensembl)
ensembl.rabbit <- useDataset(dataset = "ocuniculus_gene_ensembl", mart = ensembl)

# read human identifiers
human <- read.table("data-2022-02-21/human/MER-PC-b002_total.ReadCounts.tsv", header=TRUE, sep="\t")
human.split <- cSplit(human, "X" ,"_")[,c(6,7,8,1,2,3,4,5)]
human.filt <- human.split[human.split$X_3 == "ProteinCoding",]
human.filt[,c(4:8)] <- log10(human.filt[,c(4:8)]+1)
colnames(human.filt)[1:3] <- c("GeneID", "GeneName", "Type")

human.ids <- human.filt[,c(1,2)]
rm(human,human.split,human.filt)

# convert to sheep
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

test <- mapping.sheep[startsWith(mapping.sheep$Sheep.GeneName,"ENSOAR"),]
colnames(test)[2] <- "Sheep.GeneID"
test <- cbind(test, test$Sheep.GeneID)
colnames(test)[3] <- "Sheep.GeneName"

res <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name', values = mapping.sheep$Sheep.GeneName, mart = ensembl.sheep)

mapping <- merge(x=mapping.sheep,y=res, by.x = "Sheep.GeneName", by.y = "external_gene_name")
colnames(mapping)[3] <- "Sheep.GeneID"
mapping <- mapping[,c(2,3,1)]
mapping <- rbind(mapping,test)
mapping <- distinct(mapping)

##############################################################

# convert to rabbit
method <- "gprofiler"  
mapping.rabbit <- orthogene::convert_orthologs(gene_df = human.ids$GeneID,
                                              gene_input = "Human.GeneID", 
                                              gene_output = "columns", 
                                              input_species = "human",
                                              output_species = "rabbit",
                                              non121_strategy = "kbs",
                                              method = method) 

mapping.rabbit <- mapping.rabbit[,c(2,3)]
colnames(mapping.rabbit) <- c("Human.GeneID","Rabbit.GeneName")

test <- mapping.rabbit[startsWith(mapping.rabbit$Rabbit.GeneName,"ENSOCU"),]
colnames(test)[2] <- "Rabbit.GeneID"
test <- cbind(test, test$Rabbit.GeneID)
colnames(test)[3] <- "Rabbit.GeneName"

res <- getBM(attributes = c('ensembl_gene_id','external_gene_name'), filters = 'external_gene_name', values = mapping.rabbit$Rabbit.GeneName, mart = ensembl.rabbit)
mapping.rabbit <- merge(x=mapping.rabbit,y=res, by.x = "Rabbit.GeneName", by.y = "external_gene_name")
colnames(mapping.rabbit)[3] <- "Rabbit.GeneID"
mapping.rabbit <- mapping.rabbit[,c(2,3,1)]

mapping.rabbit <- rbind(mapping.rabbit,test)
mapping.rabbit <- distinct(mapping.rabbit)

mapping <- merge(x = mapping,y = mapping.rabbit, by = "Human.GeneID")

mapping <- merge (x=mapping, y= human.ids, by.x="Human.GeneID", by.y = "GeneID")
mapping <- mapping[,c(1,6,2:5)]
colnames(mapping)[2] <- "Human.GeneName"

write.table(mapping, file="homology-mapping.tsv", quote=FALSE, sep="\t", row.names=FALSE)


# add gene names for sheep

