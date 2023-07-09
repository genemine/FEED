# 1. Description

FEED (Feature Extraction based on gene Expression Decomposition) is a feature selection method based on gene expression decomposition for single-cell clustering. Feature selection is a critical step to determine a subset of genes to improve clustering accuracy.

# 2. Input data

A gene expression matrix is used as input, where rows represent genes and columns represent cells. 

# 3. Usage

`genes <- FEED(dataset)`

The genes selected by FEED can be used to update gene expression matrix and cluster cells.
