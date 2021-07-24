####The variable rna_mat is a gene expression matrix, or part of that.
####The variable clusters is a named character vector.
####Names of this vector are cell IDs/barcodes, 
####and elements in this vector are cluster names.
####eg. name: AAACCTGAGACCACGA-1 element: C01

#To calculate the mean, median and sd value of RNA abundance of each cell cluster
abundance_cal <- function(rna_mat,clusters){
  abundance_mean <- numeric()
  abundance_sd <- numeric()
  abundance_median <- numeric()
  rna_percell <- apply(rna_mat,2,sum)
  types <- unique(clusters)
  for(i in types){
    abundance_cell <- rna_percell[clusters[names(rna_percell)] == i]
    abundance_mean[i] <- mean(abundance_cell)
    abundance_median[i] <- median(abundance_cell)
    abundance_sd[i] <- sd(abundance_cell)
  }
  abundance <- data.frame(mean = abundance_mean,sd = abundance_sd,
                          median = abundance_median,
                          row.names=names(abundance_mean))
  return(abundance)
}  

#To calculate RNA abundance of each cell
abd_clst_cal <- function(rna_mat,clusters){
  abundance_all_cells <- data.frame()
  rna_percell <- apply(rna_mat,2,sum)
  types <- unique(clusters)
  for(i in types){
    abundance_cell <- rna_percell[clusters[names(rna_percell)] == i]
    abundance_cluster <- data.frame(names(abundance_cell),
                                    abundance_cell,
                                    rep(i,length(abundance_cell)),
                                    stringsAsFactors = F)
    abundance_all_cells <- rbind(abundance_all_cells,abundance_cluster)
  }
  colnames(abundance_all_cells) <- c("barcodes","RNA_abundance","cluster")
  abundance_all_cells$cluster <- factor(abundance_all_cells$cluster,
                                        levels = types)
  return(abundance_all_cells)
}

#To collect cell markers in the expression matrix
marker_collect <- function(rna_mat,clusters){
  cluster_names <- unique(clusters)
  marker_allcells <- data.frame()
  for(i in cluster_names){
    marker_perclstuer <- rna_mat[,clusters[colnames(rna_mat)] == i]
    marker_perclstuer <- t(marker_perclstuer)
    marker_percluster <- data.frame(rownames(marker_perclstuer),
                                    rep(i,nrow(marker_perclstuer)),
                                    marker_perclstuer,
                                    stringsAsFactors = F)
    marker_allcells <- rbind(marker_allcells,marker_percluster)
  }
  colnames(marker_allcells) <- c("Barcode","Cluster",rownames(rna_mat))
  marker_allcells$Cluster <- factor(marker_allcells$Cluster,
                                    levels = cluster_names)
  return(marker_allcells)
}  #calculate marker genes of each cell cluster


#To calculate content of markers in each cell
marker_abd_cal <- function(rna_mat,clusters){
  abundance_mean <- numeric()
  abundance_sd <- numeric()
  abundance_median <- numeric()
  types <- unique(clusters)
  for(i in types){
    mat_percluster <- rna_mat[rna_mat$Cluster == i,]
    rna_percell <- apply(rna_mat,2,sum)
    abundance_cell <- rna_percell[clusters[names(rna_percell)] == i]
    abundance_mean[i] <- mean(abundance_cell)
    abundance_median[i] <- median(abundance_cell)
    abundance_sd[i] <- sd(abundance_cell)
  }
  abundance <- data.frame(mean = abundance_mean,sd = abundance_sd,
                          median = abundance_median,
                          row.names=names(abundance_mean))
  return(abundance)
}

