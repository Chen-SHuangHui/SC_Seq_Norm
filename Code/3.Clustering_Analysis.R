library(SingleCellExperiment)
library(SC3)
library(scater)

#To filter low_quality cells
filter <- function(norm_data){
  ##The variable "norm data" refers to normalized expression data or UMI counts
  ##To recognize cells with too high mitochondrial gene expression
  ##                   or too low nuclear expression
  rna_abd1 <- apply(norm_data,2,sum)
  A_few_MT_genes <- intersect(rownames(norm_data),MT_genes)
  MT_mat <- norm_data[A_few_MT_genes[,2],]
  MT_abd <- apply(MT_mat,2,sum)
  filtered_cells1 <- names(rna_abd1[rna_abd1 > 500])
  filtered_cells2 <- names(rna_abd1[MT_abd/rna_abd1 < 0.05])
  #To recognize cells expressing CD14
  filtered_cells3 <- colnames(norm_data[,which(norm_data["CD14",] == 0)])
  #To select proper cells
  filtered_cells <- intersect(intersect(filtered_cells1,filtered_cells2),
                              filtered_cells3)
  norm_data <- norm_data[,filtered_cells]
  return(norm_data)
}

sc3_clustering <-function(norm_data,cell_types,kc){
  ##To create a SingleCellExperiment object
  ##kc: customerized k, a numeric vector
  sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(norm_data),
                  logcounts = log2(as.matrix(norm_data) + 1)), 
    colData = cell_types)
  ##To define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  ##To remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  sce <- sc3(sce, ks = kc, biology = TRUE)
  return(sce)
}

#To draw graphs
Graphing <- function(kc){
  ## kc: customized k
  jpeg(filename = "SC3_plotPCA.jpg")
  plotPCA(SC3_clusters, 
          colour_by = paste("sc3",kc,"clusters",sep = "_"), 
          size_by = paste("sc3",kc,"log2_outlier_score",sep = "_")
          )
  dev.off()
  
  jpeg(filename = "SC3_DE_genes.jpg")
  sc3_plot_de_genes(
    SC3_clusters, k = kc, 
    show_pdata = c(
      "cell_type", 
      "log10_total_features",
      paste("sc3_",kc,"_clusters",sep = ""), 
      paste("sc3_",kc,"_log2_outlier_score",sep = "")
    )
  )
  dev.off()
  
  jpeg(filename = "SC3_consensus.jpg")
  sc3_plot_consensus(
    SC3_clusters, k = kc, 
    show_pdata = c(
      "cell_type", 
      "log10_total_features",
      paste("sc3_",kc,"_clusters",sep = ""), 
      paste("sc3_",kc,"_log2_outlier_score",sep = "")
    )
  )
  dev.off()
  
  jpeg(filename = "SC3_markers.jpg")
  sc3_plot_markers(
    SC3_clusters, k = kc, 
    show_pdata = c(
      "cell_type", 
      "log10_total_features",
      paste("sc3_",kc,"_clusters",sep = ""), 
      paste("sc3_",kc,"_log2_outlier_score",sep = "")
    )
  )
  dev.off()
}

#Data pre-treatment
##To load raw gene data
library(Matrix)
###To get raw matrix, which contains gene expression data and protein data 
untar("vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.tar")
raw_mat <- as.matrix(readMM(file="filtered_feature_bc_matrix/matrix.mtx.gz"))
###To set colnames
barcodes <- read.table(file="filtered_feature_bc_matrix/barcodes.tsv.gz",
                       header=F,check.names = F)
colnames(raw_mat) <- barcodes[,1]
###To set rownames
features <- read.table(file="filtered_feature_bc_matrix/features.tsv.gz",
                       header=F,check.names = F,stringsAsFactors = F)
rownames(raw_mat) <- features[,2]
###To divide the raw matrix into 2 parts
raw_gene_mat <- raw_mat[features[V3 == "Gene",2],]
raw_pr_mat <- raw_mat[features[V3 == "Antibody",2],]
###Filtering low_quality cells
SC3_mat <- filter(raw_gene_mat)
###Setting cell types(cell types can be customerized. Here is an example)
cell_type <- data.frame(cell_type = rep("None",times = ncol(SC3_mat)),
                        row.names = colnames(SC3_mat))

##To get mitochondiral genes
MT_genes <- read.table("MT_genes.txt",stringsAsFactors = F)

#Clustering analysis and export the results
SC3_clusters <- sc3_clustering(SC3_mat,cell_type,9:16)
save(SC3_clusters,file = "SC3_clusters.RData")
Graphing(16)
sc3_export_results_xls(SC3_clusters,filename = "SC3_clusters.xls")

