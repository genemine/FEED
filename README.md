# 1. Description

FEED (Feature Extraction based on gene Expression Decomposition) is a feature selection method based on gene expression decomposition for single-cell clustering. Feature selection is a critical step to determine a subset of genes to improve clustering accuracy. 

![overview of FEED](/overview.jpg)

FEED consists of four parts: input, gene decomposition, gene correlation and gene selection. 

The input data is gene expression matrix where rows represent genes and columns represent cells. 
FEED filters out genes that are expressed in less than 6% of cells and expressed in more than 94% of cells. In addition, the filtered gene expression matrix X is log2-transformed after adding a pseudo-count 1.

First, FEED decomposes the distribution that expression of each gene follows into multiple Gaussian distributions called components. It then filters out genes with less than two components.
Second, genes are grouped into bins based on their average expression. The subsequent analysis is performed within each bin. Jensen-Shannon divergence is applied to calculate correlations between genes.
Finally, based on the gene correlation matrix, coefficient of variation (CVs) are calculated as gene importance. A permutation-based threshold calculation method is used to determine a cutoff of gene importance and thus identify marker genes.

# 2. Implementation
FEED is implemented in R. It is tested on  Linux operating system. They are freely available for non-commercial use.

# 3. Usage

`genes <- FEED(dataset)`

The genes selected by FEED can be used to update gene expression matrix and cluster cells.The demo input data are provided in the folder 'dataset'.

# 4. Contact

If any questions, please do not hesitate to contact us at: hongdong@csu.edu.cn
