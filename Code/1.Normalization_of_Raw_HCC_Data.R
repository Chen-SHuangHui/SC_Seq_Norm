#To load raw counts
raw_HCC_data <- read.table(file='GSE98638_HCC.TCell.S5063.count.txt.gz',
                       header = TRUE,stringsAsFactors = FALSE,
                       check.names = F);
raw_HCC_data$geneID <- NULL;
na_index <- is.na(raw_HCC_data$symbol)
raw_HCC_data$symbol[na_index] <- paste("na",1:sum(na_index),sep = "-")
####This line can be removed if there is no NA in raw expression matrix
rownames(raw_HCC_data) <- raw_HCC_data$symbol;                    
raw_HCC_data$symbol <- NULL;

#1. Normalization by DESeq
library("DESeq")
condition <- factor(rep("untreated",5063))
cds <- newCountDataSet(raw_HCC_data, condition)
cds <- estimateSizeFactors(cds)
#To extract the normalized counts
DESeq_norm <- counts(cds, normalized=TRUE)

#2. Normalization by Scran
raw_counts <- as.matrix(raw_HCC_data)

library("scran")
larger.sce <- SingleCellExperiment(list(counts=raw_counts))
clusters <- quickCluster(larger.sce, min.size=100)
larger.sce <- computeSumFactors(larger.sce, cluster=clusters)
larger.sce <- normalize(larger.sce)
##To extract the normalized counts: method_1
##Scran_norm <- matrix(ncol=ncol(raw_counts),nrow=nrow(raw_counts))
##size_factor <- larger.sce@int_colData$size_factor
##for(i in 1:ncol(counts(larger.sce))){
##  Scran_norm[,i] <- raw_counts[,i]/size_factor[i]
##}
##To extract the normalized counts: method_2 (recommended)
Scran_norm <- logcounts(larger.sce)^2-1

#3. Normalization by ISnorm
source("ISnorm_function.R")
mat <- as.matrix(raw_HCC_data)                                 
gene_dis<-calculate.dis(mat=mat,detection_rate=0.9,ncore=4)
spike_candidate<-dbscan.pick(dis=gene_dis,ngene=(1:floor(nrow(gene_dis)/25))*5,solution=100)
candidate_res<-candidate.norm(mat=mat,spike_candidate=spike_candidate,ncore=4)
candidate_res$inst["PTS89-0205",1] <- 0
####PTS89-0205 might cause error when ISnorm calculates instability scores
ISnorm_res<-opt.candidate(mat=mat,candidate_res=candidate_res,threshold=0.1,switch_check=2)
#To extract the normalized counts
ISnorm_norm <- ISnorm_res$normalized
