source("2.Calculation_of_RNA_abundance.R")
HCC_Cells_Table <- read.table(file = "Individual_Cells_Data_Table.txt",
                          stringsAsFactors = F,header = T,check.names = F)
HCC_clusters <- HCC_Cells_Table[,"majorCluster"]
names(HCC_clusters) <- HCC_Cells_Table[,"UniqueCell_ID"]
HCC_clusters <- HCC_clusters[HCC_clusters != "unknown"]
cellnum <- as.numeric(table(HCC_clusters))
DESeq_clst_abd <- abd_clst_cal(DESeq_norm,HCC_clusters)
ISnorm_clst_abd <- abd_clst_cal(ISnorm_norm,HCC_clusters)
Scran_clst_abd <- abd_clst_cal(Scran_norm,HCC_clusters)

boxplot(RNA_abundance~cluster,data = ISnorm_clst_abd,range = 0,
        names = paste("C",1:11,"(n=",cellnum,")",sep=""),
        width=cellnum,col = c(rep("yellow",5),rep("green",6)),
        ylab = "RNA_abundance",xlab = "cell_clusters",cex.axis = 0.65)

boxplot(RNA_abundance~cluster,data = DESeq_clst_abd,range = 0,
        names = paste("C",1:11,"(n=",cellnum,")",sep=""),
        width=cellnum,col = c(rep("yellow",5),rep("green",6)),
        ylab = "RNA_abundance",xlab = "cell_clusters",cex.axis = 0.65)

boxplot(RNA_abundance~cluster,data = Scran_clst_abd,range = 0,
        names = paste("C",1:11,"(n=",cellnum,")",sep=""),
        width=cellnum,col = c(rep("yellow",5),rep("green",6)),
        ylab = "RNA_abundance",xlab = "cell_clusters",cex.axis = 0.65)