library(reticulate)
library(Seurat)
library(ggplot2)

#set reticulate python path to correct python path
Sys.setenv(RETICULATE_PYTHON = '/opt/anaconda3/envs/r-reticulate/bin/python')
#import the scGene2Vec library (after installed in r-reticulate environment) as g2v
g2v <- import('scGene2Vec.scGene2Vec')



#First, generate a SeuratObject, and find variable genes 
#or use SCTransform to get variable genes under scale.data slot

pbmc <- Read10X('/Users/beth/Downloads/filtered_gene_bc_matrices 3/hg19/')
pbmc <- CreateSeuratObject(pbmc)
MT.genes <- grep("^MT-", rownames(pbmc), value=TRUE)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = MT.genes)
mt.threshold <- 10
Count.threshold <- 6000
pbmc <- subset(pbmc, subset = percent.mt < mt.threshold & 
                      nCount_RNA < Count.threshold)
pbmc <- SCTransform(pbmc)

#make a correlation table for the variable genes
cor_table <- cor(t(pbmc@assays$SCT@scale.data))

#use g2v to make a table of vectors for each gene from a gene2vec model
#window = 1/2 * the number of most highly correlated genes to include in each gene's group (sentence in word2vec)
#size = the length of the desired vector for each gene (also size in word2vec)
vector_table <- g2v$make_vector_table(as.data.frame(cor_table), window=as.integer(10), size=as.integer(10))

#function to normalize a vector
norm_vec <- function(x) sqrt(sum(x^2))

#function to get cell embedding for each cell and normalize (to use cosine similarity)
g2v_cell_embeddings_fast <- function(SeuratObj, gene2vec){
  counts <- SeuratObj@assays$SCT@scale.data
  cell_embeddings <- t(counts) %*% t(as.matrix(gene2vec))
  cell_embeddings <- apply(cell_embeddings, 1, function (x) x/norm_vec(x))
  colnames(cell_embeddings) <- colnames(counts)
  return(t(cell_embeddings))
}

#get cell embeddings based on gene vectors
cell_embeddings <- g2v_cell_embeddings_fast(pbmc, vector_table)

#add this reduction data object to the SeuratObject
reduction.data <- CreateDimReducObject(embeddings = cell_embeddings, 
                                       assay = "SCT", key = "g2v")
pbmc[["g2v"]] <- reduction.data 

#get pca reduction for comparison
pbmc <- RunPCA(pbmc)

#continue with analysis (Umap, Neighbors, Clustering, ect. using new g2v reduction)

#Finding optimal resolution for clusters by maximum mean silhoutte score
#Shilhoutte score and number of clusters vs. resolution plots -------------------------------------------
#the 3 functions below can be used to make the silhoutte score and number of clusters vs. resolution plots 
#they require cluster and purrr packages to be loaded
mean_sil <- function(SeuratObj, r, reduction, dims) {
  SeuratObj <-  FindClusters(SeuratObj, verbose = FALSE, resolution = r)  
  dist.matrix <- dist(x = Embeddings(object = SeuratObj[[reduction]])[, 1:dims])
  clusters <- SeuratObj$seurat_clusters
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  return(mean(sil[,3]))
}

number_clusters <- function(SeuratObj, r) {
  SeuratObj <-  FindClusters(SeuratObj, verbose = FALSE, resolution = r)  
  return(length(levels(SeuratObj@active.ident)))
}

make_res_plot <- function(SeuratObj, reduction, dims){
  SeuratObj<-FindNeighbors(SeuratObj,reduction,dims=1:dims)
  mean_sil_score <- map(seq(0.1,2,0.05), mean_sil, SeuratObj=SeuratObj, reduction=reduction, dims=dims)
  num_clusters <-  map(seq(0.1,2,0.05), number_clusters, SeuratObj=SeuratObj)
  
  num_clusters <- unlist(num_clusters)
  mean_sil_score <- unlist(mean_sil_score)
  
  df <- data.frame(seq(0.1,2,0.05),num_clusters,mean_sil_score)
  colnames(df) <- c('resolution','num_clusters', 'mean_sil_score')
  
  p = ggplot(df, aes(x=resolution)) +
    
    geom_line( aes(y=num_clusters), color = 'blue') +
    geom_line( aes(y=mean_sil_score*25), color = 'red') +
    
    scale_y_continuous(
      limits = c(-0,25),
      # Features of the first axis
      name = "Number of Clusters",
      
      # Add a second axis and specify its features
      sec.axis = sec_axis(~./25, name="Mean Sil Score") ) +
    
    theme_classic() +
    
    theme(axis.title.y = element_text(color = 'blue', size=12),
          axis.title.y.right = element_text(color = 'red', size=12)
    )
  
  return(p)
}

#Function to make a pdf (only way to plot here) of a silhoullete plot for all the spots/cells
#with the current clustering result (uses cluster package)
make_sil_plot_pdf <- function(file_name, SeuratObj, reduction_method_used, dims_used){
  dist.matrix <- dist(x = Embeddings(object = se[[reduction_method_used]])[, 1:dims_used])
  clusters <- SeuratObj$seurat_clusters
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  pdf(file=file_name)
  plot(sil)
  dev.off()
}

#Function for plotting the heatmap for the top x number of genes in each cluster ---------------
#usess Seurat's DoHeatmap function
plot_top_markers_heatmap <- function(SeuratObj, all_markers, num_genes_per_cluster){
  num_cluster = length(unique(all_markers$cluster))
  temp <- head(all_markers[all_markers$cluster==0  & (all_markers$avg_log2FC>0),],num_genes_per_cluster) 
  for (i in seq(1,num_cluster)){
    temp <- rbind(temp, head(all_markers[all_markers$cluster==i  & (all_markers$avg_log2FC>0),],num_genes_per_cluster))
  }
  p1 <- DoHeatmap(SeuratObj, features = temp$gene)
  return(p1)
}


pca_res_plot <- make_res_plot(pbmc,'pca',10)
pca_res_plot

g2v_res_plot <- make_res_plot(pbmc,'g2v',10)
g2v_res_plot

pbmc <- FindNeighbors(pbmc, 'pca', dims=1:10)
pbmc <- FindClusters(object = pbmc, verbose = TRUE, resolution=1.8)
make_sil_plot_pdf('pca_pbmc',pbmc,'pca',10)

allmarkers <- FindAllMarkers(pbmc, only.pos=T, logfc.threshold = 0.25, min.pct = 0.1)
pca_heat <- plot_top_markers_heatmap(pbmc, allmarkers, 10)

pbmc <- RunUMAP(pbmc,'pca',dims = 1:10)
pca_umap <- DimPlot(pbmc, label = TRUE, reduction = "umap")

pbmc <- FindNeighbors(pbmc, 'g2v', dims=1:10)
pbmc <- FindClusters(object = pbmc, verbose = TRUE, resolution=1)
make_sil_plot_pdf('g2v_pbmc',pbmc,'g2v',10)

pbmc <- RunUMAP(pbmc,'g2v',dims = 1:10)
g2v_umap <- DimPlot(pbmc, label = TRUE, reduction = "umap")

allmarkers <- FindAllMarkers(pbmc, only.pos=T, logfc.threshold = 0.25, min.pct = 0.1)
g2v_heat <- plot_top_markers_heatmap(pbmc, allmarkers, 10)

library(caret)

cm <- confusionMatrix(pbmc@meta.data$SCT_snn_res.1,pbmc@meta.data$SCT_snn_res.1.8,dnn=c('g2v','pca'))

test <- FindMarkers(pbmc, ident.1=0, ident.2 = c(3,5,6,7,8,9,14), only.pos=T)

library(SeuratDisk)
library(SeuratData)
