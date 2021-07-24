#To remove gene with zero-counts in all cells
zero_index <- apply(SC3_mat,1,sum) != 0
filtered_gene_mat <- SC3_mat[zero_index,]

#Marker genes
T_marker_genes <- c("LEF1","CCR7","SELL","TCF7")
B_marker_genes <- c("BMI1","CD24","CD69","CR1","ENTPD1","FCER2",
                    "FCGR2B","FCGRT","IGHD")

#To divide cells into clusters
##Method1
###library(XLConnect)
###connect <- loadWorkbook('SC3_clusters.xls')
###clst <- readWorksheet(connect,'Cells')
##Method2
library(readxl)
clst <- as.data.frame(read_excel("SC3_clusters.xls",sheet = "Cells"))

clst <- clst[,c("...1","sc3_16_clusters")]
PBMC_clusters <- paste("C",clst[,2],sep = "0")
PBMC_clusters <- gsub("C01","C1",PBMC_clusters)
PBMC_clusters[PBMC_clusters == "C1"] <- "C01"
names(PBMC_clusters) <- clst[,1]

#To construct matrices containing marker gene expression data
#                              or  surface protein TotalSeq
library(reshape2)
##Marker genes
naive_markers_mat <- filtered_gene_mat[c(T_marker_genes,B_marker_genes),]
naive_marker_content <- marker_collect(naive_markers_mat,PBMC_clusters)
##Protein TotalSeq
pr_content <- marker_collect(raw_pr_mat,PBMC_clusters)
##To melt the information
merged_content <- cbind(pr_content,
                             naive_marker_content[rownames(pr_content),
                                                  c(T_marker_genes,
                                                    B_marker_genes)])
merged_content$Cluster <- as.character(merged_content$Cluster)
log_merged_content <- cbind(merged_content[,1:2],
                            log2(merged_content[,-(1:2)]+1))
log_T_content <- log_merged_content[,c("Barcode","Cluster",T_marker_genes)]
log_T_content <- melt(log_T_content,id.vars = c("Barcode","Cluster"))
colnames(log_T_content)[3:4] <- c("Gene","Log_UMI")
##RNA abundance of each cell
abd_clst_PBMC <- abd_clst_cal(filtered_gene_mat,PBMC_clusters)
log_abd_clst_PBMC <- transform(abd_clst_PBMC, 
                               log_RNA_abundance = log2(RNA_abundance+1))

#Significant Tests
####For example,to compare the content of gene LEF1 in cluster C02 and C16
with(log_merged_content,t.test(log_merged_content[Cluster == "C02","LEF1"],
                               log_merged_content[Cluster == "C16","LEF1"]))

#To draw graphs
library(ggplot2)
##CD3
ggplot(log_merged_content, aes(x = Cluster, y = CD3_TotalSeqC,
                              fill = Cluster)) +
  geom_violin(trim = T) +
  labs(title = "CD3")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(limits=c(paste("C0",c(1:5,8:9),sep = ""),
                            paste("C",c(10,11,13,16),sep = "")
                            ))+
  scale_fill_manual(values = rep(c("#E69F00", "#56B4E9","#E69F00"),
                                 times = c(5,5,1)))+
  xlab("Cluster")+
  ylab("Log CD3 Content")+
  stat_summary(fun.y=mean, geom="point", size=1, color="black")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

##CD19
ggplot(log_merged_content, aes(x = Cluster, y = CD19_TotalSeqC,
                               fill = Cluster)) +
  geom_violin(trim = T) +
  labs(title = "CD19")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(limits=c(paste("C0",c(1:5,8:9),sep = ""),
                            paste("C",c(10,11,13,16),sep = "")
  ))+
  scale_fill_manual(values = rep(c("#56B4E9", "#E69F00","#56B4E9"),
                                times = c(5,2,4)))+
  xlab("Cluster")+
  ylab("Log CD19 Content")+
  stat_summary(fun.y=mean, geom="point", size=1, color="black")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

##CD4
ggplot(log_merged_content, aes(x = Cluster, y = CD4_TotalSeqC,
                               fill = Cluster)) +
  geom_violin(trim = T) +
  labs(title = "CD4")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(limits=c(paste("C0",c(1:5),sep = ""),"C16"))+
  scale_fill_manual(values = rep(c("#E69F00","#56B4E9",
                                   "#E69F00","#56B4E9"),
                                 times = c(1,1,3,1)))+
  xlab("Cluster")+
  ylab("Log CD4 Content")+
  theme(legend.position = "none")+
  stat_summary(fun.y=mean, geom="point", size=1, color="black")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

##CD8A
ggplot(log_merged_content, aes(x = Cluster, y = CD8a_TotalSeqC,
                               fill = Cluster)) +
  geom_violin(trim = T) +
  scale_x_discrete(limits=c(paste("C0",c(1:5),sep = ""),"C16"))+
  scale_fill_manual(values = rep(c("#56B4E9","#E69F00",
                                   "#56B4E9","#E69F00"),
                                times = c(1,1,3,1)))+
  labs(title = "CD8A")+
  xlab("Cluster")+
  ylab("Log CD8A Content")+
  theme(legend.position = "none")+
  stat_summary(fun.y=mean, geom="point", size=1, color="black")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") 

##T markers
ggplot(log_T_content, aes(x = Cluster, y = Log_UMI,
                               fill = Gene)) +
  geom_boxplot(outlier.size = 1) +
  labs(title = "Naive T markers in T Clusters")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Cluster",
                   limits=c("C01","C02","C03",
                            "C04","C05","C16"))+
  scale_y_continuous(name = "Log UMI")+
  stat_summary(fun.y=mean, geom="line", size=1, color="black")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

##RNA abundance
ggplot(log_abd_clst_PBMC, aes(x = cluster, y = log_RNA_abundance,
                          fill = cluster)) +
  geom_boxplot(color = "black") +
  labs(title = "RNA_abundance")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(name = "Cluster",
                   limits=c("C01","C02","C03",
                            "C04","C05","C16"))+
  scale_y_continuous(name = "Log RNA abundance")+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
