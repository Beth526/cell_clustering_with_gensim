library(reticulate)
library(Seurat)

#set reticulate python path to correct python path
Sys.setenv(RETICULATE_PYTHON = '/opt/anaconda3/envs/r-reticulate/bin/python')
#import the scGene2Vec library (after installed in r-reticulate environment) as g2v
g2v <- import('scGene2Vec.scGene2Vec')

#First, generate a SeuratObject, and find variable genes 
#or use SCTransform to get variable genes under scale.data slot

#make a correlation table for the variable genes
cor_table <- cor(t(SeuratObject@assays$SCT@scale.data))

#use g2v to make a table of vectors for each gene from a gene2vec model
#window = 1/2 * the number of most highly correlated genes to include in each gene's group (sentence in word2vec)
#size = the length of the desired vector for each gene (also size in word2vec)
vector_table <- g2v$make_vector_table(as.data.frame(cor_table), window=10, size=10)

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
cell_embeddings <- g2v_cell_embeddings_fast(SeuratObject, vector_table)

#add this reduction data object to the SeuratObject
reduction.data <- CreateDimReducObject(embeddings = cell_embeddings, 
                                       assay = "SCT", key = "g2v")
SeuratObject[["g2v"]] <- reduction.data 

#continue with analysis (Umap, Neighbors, Clustering, ect. using new g2v reduction)
