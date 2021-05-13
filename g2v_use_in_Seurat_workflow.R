library(reticulate)
library(Seurat)
library(ggplot2)
library(SeuratData)
library(cluster)
library(purrr)

v#set reticulate python path to correct python path
Sys.setenv(RETICULATE_PYTHON = '/opt/anaconda3/envs/r-reticulate/bin/python')
#import the scGene2Vec library (after installed in r-reticulate environment) as g2v
g2v <- import('scGene2Vec.scGene2Vec')

#First, generate a SeuratObject, and find variable genes 
#or use SCTransform to get variable genes under scale.data slot
InstallData("pbmc3k")
data("pbmc3k")

#used this previous download earlier
#pbmc <- Read10X('/Users/beth/Downloads/filtered_gene_bc_matrices 3/hg19/')

pbmc3k <- CreateSeuratObject(pbmc3k)

MT.genes <- grep("^MT-", rownames(pbmc3k), value=TRUE)
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, features = MT.genes)
mt.threshold <- 10
Count.threshold <- 6000
pbmc3k <- subset(pbmc3k, subset = percent.mt < mt.threshold & 
                      nCount_RNA < Count.threshold)
pbmc3k <- SCTransform(pbmc3k)

#make a correlation table for the variable genes
cor_table <- cor(t(pbmc3k@assays$SCT@scale.data))

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
cell_embeddings <- g2v_cell_embeddings_fast(pbmc3k, vector_table)

#add this reduction data object to the SeuratObject
reduction.data <- CreateDimReducObject(embeddings = cell_embeddings, 
                                       assay = "SCT", key = "g2v_")
pbmc3k[["g2v"]] <- reduction.data 

#get pca reduction for comparison
pbmc3k <- RunPCA(pbmc3k)

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
  dist.matrix <- dist(x = Embeddings(object = SeuratObj[[reduction_method_used]])[, 1:dims_used])
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

#-----------------------------------------------compare res plots
pca_res_plot <- make_res_plot(pbmc3k,'pca',10)
pca_res_plot

g2v_res_plot <- make_res_plot(pbmc3k,'g2v',10)
g2v_res_plot

#-----------------------------------------------pca with 22 clusters
pbmc3k <- FindNeighbors(pbmc3k, 'pca', dims=1:10)
pbmc3k <- FindClusters(object = pbmc3k, verbose = TRUE, resolution=3)
make_sil_plot_pdf('pca_pbmc',pbmc3k,'pca',10)

allmarkers <- FindAllMarkers(pbmc3k, only.pos=T, logfc.threshold = 0.25, min.pct = 0.1)
pca_heat <- plot_top_markers_heatmap(pbmc3k, allmarkers, 10)

pbmc3k <- RunUMAP(pbmc3k,'pca',dims = 1:10)
pca_umap <- DimPlot(pbmc3k, label = TRUE, reduction = "umap")

#--------------------------------------------------g2v with 22 clusters
pbmc3k <- FindNeighbors(pbmc3k, 'g2v', dims=1:10)
pbmc3k <- FindClusters(object = pbmc3k, verbose = TRUE, resolution=2)
make_sil_plot_pdf('g2v_pbmc',pbmc3k,'g2v',10)

pbmc3k <- RunUMAP(pbmc3k,'g2v',dims = 1:10)
g2v_umap <- DimPlot(pbmc3k, label = TRUE, reduction = "umap")

allmarkers <- FindAllMarkers(pbmc3k, only.pos=T, logfc.threshold = 0.25, min.pct = 0.1)
g2v_heat <- plot_top_markers_heatmap(pbmc3k, allmarkers, 10)

#--------------------------------------------------cluster relationships heatmap
library(caret)

cm <- confusionMatrix(pbmc3k@meta.data$SCT_snn_res.2,pbmc3k@meta.data$SCT_snn_res.3,dnn=c('g2v','pca'))

heatmap(cm$table, xlab='pca',ylab='g2v')
#--------------------------------------------------Azimuth web portal integration results
predictions <- read.delim('azimuth_pred.tsv', row.names = 1)
pbmc3k <- AddMetaData(
       object = pbmc3k,
       metadata = predictions)

par(mar=c(10,1,1,1))

g2v_cm <- table(pbmc3k@meta.data$SCT_snn_res.2,pbmc3k@meta.data$predicted.celltype.l2)
g2v_cm <- prop.table(g2v_cm, margin=2)
heatmap(g2v_cm,ylab='g2v')


pca_cm <- table(pbmc3k@meta.data$SCT_snn_res.3,pbmc3k@meta.data$predicted.celltype.l2)
pca_cm <- prop.table(pca_cm, margin=2)
heatmap(pca_cm,ylab='pca')
