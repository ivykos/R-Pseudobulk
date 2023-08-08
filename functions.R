library(Seurat)
library(Matrix)
library(tidyverse)


#Define pseudobulk function
pseudobulk <- function(obj){
  matrix <- GetAssayData(obj, slot="counts")
  df <- data.frame(rownames(matrix))
  for (i in unique(obj@meta.data$projid_region)){
    matrix <- GetAssayData(subset(obj, subset = projid_region == as.character(i)))
    bulk <- data.frame(Matrix::rowSums(matrix))
    df[as.character(i)] <- bulk[1]
    
  }
  df <- data.frame(df, row.names = 1)
  return(df)
}

#Define function to make DESeq design matrix
coldata <- function(obj){
  meta <- obj@meta.data
  meta <- meta %>% distinct(projid_region, .keep_all=TRUE)
  col_data <- data.frame(meta)
  col_data <- col_data %>% remove_rownames %>% column_to_rownames(var="projid_region")
  return(col_data)
  write.csv(col_data, "coldata.csv")
}