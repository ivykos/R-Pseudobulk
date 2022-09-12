library(Seurat)
library(Matrix)
library(tidyverse)



#Create pseudo-bulk counts file
pseudobulk <- function(obj){
  matrix <- GetAssayData(obj, slot="counts")
  df <- data.frame(rownames(matrix))
  samples <- as.character(unique(obj@meta.data$projid))
  for (i in samples){
    matrix <- GetAssayData(subset(obj, subset = projid == i), slot="counts")
    bulk <- data.frame(Matrix::rowSums(matrix))
    df[as.character(i)] <- bulk[1]
    
  }
  df <- data.frame(df, row.names = 1)
  names(df) <- sub("^X","",names(df)) 
  return(df)
  write.csv(df, "bulk.csv")
}

#Generate Deseq ColData
coldata <- function(obj){
  meta <- obj@meta.data
  meta <- meta %>% distinct(projid, .keep_all=TRUE)
  col_data <- data.frame(meta)
  col_data <- col_data %>% remove_rownames %>% column_to_rownames(var="projid")
  return(col_data)
  write.csv(col_data, "coldata.csv")
}