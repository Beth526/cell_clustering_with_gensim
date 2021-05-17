import pandas as pd 
from gensim.models import Word2Vec
import numpy as np

def make_sentences(gene_corr, window):
    '''
    Makes the "sentences" for Word2Vec, which are in this case group of genes that are most highly correlated with the gene of interest (center gene)
    
    Parameters:
    window = 1/2 the number of most highly correlated genes to include in the group with the gene of interest
    
    Returns:
    A list of "sentences" or groups of genes to use as input for Word2Vec
    '''
    sentences = []
    for i in range(gene_corr.shape[1]):
        top_corr = gene_corr.iloc[i].sort_values()[-(window*2):-1].index
        sentence = list(top_corr[:window]) + [gene_corr.columns[i]] + list(top_corr[window:])
        sentences.append(sentence)
    return sentences



def make_vector_table(gene_corr,window=10,size=10):
    '''
    Gets the gene vectors using gensim's Word2Vec implementation
    
    Parameters:
    window = 1/2 the number of most highly correlated genes to include in the group with the gene of interest
    size = same as size parameter in Word2Vec (corresponds to length of gene vectors returned)
    
    Returns:
    Dataframe of gene vectors
    '''
    sentences = make_sentences(gene_corr,window)
    gene2vec_model = Word2Vec(sentences,
                    min_count=1,
                    window=window,
                    size=size,
                    workers=1,
                    seed=1)
    df = np.zeros((size,gene_corr.shape[1]))
    df = pd.DataFrame(df)
    df.columns = gene_corr.columns
    for i in list(gene_corr.columns):
        df[i] = gene2vec_model[i]
        
    return df




