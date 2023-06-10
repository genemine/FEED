library(SingleCellExperiment)
library(tidyverse)
library(Rtsne)
library(dbscan)
library(fpc)
library(Seurat)
library(nortest)
library(mclust)
library(parallel)
library(doRNG)

library(scater)
library(RANN)
library(igraph)
library(mclust)
library(cluster)

library(SCMarker)
library(M3Drop)
library(FEAST)
library(DUBStepR)

setwd("/home/zhangchao_204712149/my_idea/R/")

####################################################################
#                     clustering selection method                  #
####################################################################
HC_clustering <- function(dataset, K) {
  distance <- as.dist(1 - cor(dataset))
  h_tree <- hclust(distance, method = "complete")
  cluster_result <- cutree(h_tree, k = K)
  return(cluster_result)
}

Louvain_clustering <- function(dataset,gene_list) {
  SeuObj = CreateSeuratObject(counts=data)
  SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
  all.gene = rownames(SeuObj)
  SeuObj <- ScaleData(SeuObj, features = all.gene, verbose = FALSE)
  if(ncol(data)<50){
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE,npcs=20)
  }else{
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE)
  }
  if(length(gene_list)<10){
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE, dims = 1:(length(gene_list)-1))
  }else{
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE)
  }
  SeuObj <- FindClusters(SeuObj, verbose = FALSE)
  cluster_ident=as.numeric(Idents(SeuObj))
  return(cluster_ident)
}

GMM_clustering <- function(dataset) {
  # dataset <- log2(dataset+1)
  pc_res <-  prcomp(t(dataset))$x
  tmp_pca_mat = pc_res[, 1:10]
  res <- Mclust(tmp_pca_mat)
  # cluster_result <-  apply(res$z, 1, which.max)
  cluster_result <- res$classification
  return(cluster_result)
}


####################################################################
#                     sota gene selection method                   #
####################################################################
SCMarker <- function(data){
  Res <- ModalFilter(data = data, geneK = 10, cellK = 10)
  Res <- GeneFilter(obj = Res)
  Res <- getMarker(obj = Res, k = 300, n = 30)
  scMarker_genes <- Res$marker
  return(scMarker_genes)
}

#HVG
HVG <- function(data){
  s <- 2000
  se_obj <- CreateSeuratObject(data)
  se_obj <- NormalizeData(object = se_obj)
  se_obj <- FindVariableFeatures(object = se_obj,nfeatures = s)
  Seurat_genes=se_obj@assays[["RNA"]]@var.features
  return(Seurat_genes)
}

Expr <- function(data){
  meanExp = apply(data,1,mean)
  meanExp = meanExp[order(meanExp,decreasing = TRUE)]
  s <- 2000
  Expr_genes = names(meanExp)[1:s]
  return(Expr_genes)
}

# FEAST
feast <- function(data,K){
  s <- 2000
  Y <- process_Y(data,thre = 2)
  ixs  = FEAST_fast(Y, k=K)
  return(rownames(Y)[ixs][1:s])
}

Gini_Filtering <- function(ExprM.RawCounts) {
  minCellNum <-  3
  minGeneNum <- 2000
  expressed_cutoff <- 1
  #expressed_cutoff:value lower than this value could be just noise.
  ExpressedinCell_per_gene = apply(ExprM.RawCounts,
                                   1,
                                   function(x)
                                     length(x[x > expressed_cutoff]))
  nonMir = grep("MIR|Mir", rownames(ExprM.RawCounts), invert = T)  # because Mir gene is usually not accurate
  Genelist = intersect(rownames(ExprM.RawCounts)[nonMir], 
                       rownames(ExprM.RawCounts)[ExpressedinCell_per_gene >= minCellNum])
  ExpressedGene_per_cell = apply(ExprM.RawCounts[Genelist, ], 2, function(x)
    length(x[x > 0]))
  ExprM.RawCounts.filter = ExprM.RawCounts[Genelist, ExpressedGene_per_cell >= minGeneNum]
  return(list('raw' = ExprM.RawCounts.filter))
}

GiniClust_Fitting <- function(ExprM.RawCounts.filter){
  
  # Calculate gini index
  calcul.gini = function(x, unbiased = TRUE, na.rm = FALSE){
    if (!is.numeric(x)){
      warning("'x' is not numeric; returning NA")
      return(NA)
    }
    if (!na.rm && any(na.ind = is.na(x)))
      stop("'x' contain NAs")
    if (na.rm)
      x = x[!na.ind]
    n = length(x)
    mu = mean(x)
    N = if (unbiased) n * (n - 1) else n * n
    ox = x[order(x)]
    dsum = drop(crossprod(2 * 1:n - n - 1,  ox))
    dsum / (mu * N)
  }
  
  gini = apply(as.data.frame(ExprM.RawCounts.filter), 1, function(x){calcul.gini(as.numeric(x)) } )    #theoretically, gini have very low chance to have a 1 value
  GiniIndex = as.data.frame(cbind(1:dim(ExprM.RawCounts.filter)[1], gini))
  
  Maxs          = apply(ExprM.RawCounts.filter,1,max)
  Means         = apply(ExprM.RawCounts.filter,1,mean)
  log2.Maxs     = log2(Maxs+0.1)
  ExprM.Stat1   = as.data.frame(cbind(Maxs,GiniIndex$gini,log2.Maxs))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs")
  ExprM.Stat1 = ExprM.Stat1[ExprM.Stat1$log2.Maxs>0 & ExprM.Stat1$log2.Maxs<=20 ,] 
  log2.Maxs = ExprM.Stat1$log2.Maxs
  Gini      = ExprM.Stat1$Gini
  Maxs      = ExprM.Stat1$Maxs
  
  # fitting in max-gini space 
  Gini.loess.fit        = loess(Gini~log2.Maxs, span=0.9, degree=1)
  Normlized.Gini.Score  = Gini.loess.fit$residuals   #residuals = Gini - Gini.fitted
  Gini.fitted           = Gini.loess.fit$fitted    
  ExprM.Stat1           = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs")], Normlized.Gini.Score, Gini.fitted))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs", "Norm.Gini", "Gini.fitted")
  ### remove 25% of first round outlier genes, do second round loess
  Gini.loess.fit.residual = residuals(Gini.loess.fit)                               
  thresh.outlier = quantile(Gini.loess.fit.residual[Gini.loess.fit.residual>0], 0.75) 
  id.genes.loess.fit = which(Gini.loess.fit.residual < thresh.outlier)               
  id.outliers.loess.fit = which(Gini.loess.fit.residual >= thresh.outlier)          
  log2.Maxs.genes = log2.Maxs[id.genes.loess.fit]                                   
  log2.Maxs.outliers = log2.Maxs[id.outliers.loess.fit]                            
  Gini.loess.fit.2 = loess(Gini[id.genes.loess.fit]~log2.Maxs[id.genes.loess.fit], span=0.9, degree = 1)
  Gini.loess.fit.2.predict = predict(Gini.loess.fit.2)  
  
  #plot second round fit
  Gini.loess.fit.2.x.y = cbind(log2.Maxs.genes,Gini.loess.fit.2.predict)
  Gini.loess.fit.2.x.y.uniq = as.data.frame(unique(Gini.loess.fit.2.x.y))
  Gini.loess.fit.2.x.y.uniq = Gini.loess.fit.2.x.y.uniq[order(Gini.loess.fit.2.x.y.uniq[,1]),]
  log2.Maxs.genes.sorted = log2.Maxs.genes[order(log2.Maxs.genes)]                   
  Gini.loess.fit.2.predict.sorted = Gini.loess.fit.2.predict[order(log2.Maxs.genes)] 
  #using Gini.loess.fit.2 as model, predict gini value for those outlier which are not used for build model.
  #for each max in outliers set, find the id of max value which is most close in fitted data set
  loc.outliers = apply(matrix(log2.Maxs.outliers),1,function(x){
    if(x<max(log2.Maxs.genes.sorted)){
      return(which(log2.Maxs.genes.sorted>=x)[1])
    }else{
      return(which.max(log2.Maxs.genes.sorted))
    }})                
  #check the results
  outlier_max_in_fit <- cbind(log2.Maxs.outliers, loc.outliers, log2.Maxs.genes.sorted[loc.outliers])
  
  #based on Gini.loess.fit.2, predict outliers which was not used for fitting
  Gini.outliers.predict = apply(cbind(seq(length(log2.Maxs.outliers)),log2.Maxs.outliers),1,function(x){
    id = x[1]
    value = x[2]
    if(value == log2.Maxs.genes.sorted[loc.outliers[id]]){
      return(as.numeric(Gini.loess.fit.2.x.y.uniq[which(Gini.loess.fit.2.x.y.uniq$log2.Maxs.genes>=value)[1],2]))
    }else{
      if(loc.outliers[id]>1){
        return(Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1]+(Gini.loess.fit.2.predict.sorted[loc.outliers[id]]-Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1])*(value-log2.Maxs.genes.sorted[loc.outliers[id]-1])/(log2.Maxs.genes.sorted[loc.outliers[id]]-log2.Maxs.genes.sorted[loc.outliers[id]-1]))
      }else{
        return(Gini.loess.fit.2.predict.sorted[2]-(Gini.loess.fit.2.predict.sorted[2]-Gini.loess.fit.2.predict.sorted[1])*(log2.Maxs.genes.sorted[2]-value)/(log2.Maxs.genes.sorted[2]-log2.Maxs.genes.sorted[1]))
      }
    }
  })
  
  #plot outliers predict results
  outliers.precit.x.y.uniq = as.data.frame(unique(cbind(log2.Maxs.outliers, Gini.outliers.predict)))
  colnames(outliers.precit.x.y.uniq) = colnames(Gini.loess.fit.2.x.y.uniq)
  Gini.loess.fit.2.full.x.y.uniq = rbind(Gini.loess.fit.2.x.y.uniq, outliers.precit.x.y.uniq)
  
  #calcualte Normlized.Gini.Score2
  Normlized.Gini.Score2                        = rep(0,length(Gini.loess.fit.residual))               
  Normlized.Gini.Score2[id.genes.loess.fit]    = residuals(Gini.loess.fit.2)                         
  Normlized.Gini.Score2[id.outliers.loess.fit] = Gini[id.outliers.loess.fit] - Gini.outliers.predict 
  
  Gini.fitted2           = Gini - Normlized.Gini.Score2         
  ExprM.Stat1            = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs", "Gini.fitted", "Norm.Gini" )], Gini.fitted2, Normlized.Gini.Score2))
  colnames(ExprM.Stat1)  = c("Maxs","Gini","log2.Maxs", "Gini.fitted","Norm.Gini",  "Gini.fitted2", "Norm.Gini2")
  Gini.pvalue            = pnorm(-abs(scale(ExprM.Stat1$Norm.Gini2, center=TRUE,scale=TRUE)))
  ExprM.Stat2            = cbind(ExprM.Stat1, Gini.pvalue)  #first time use ExprM.Stat2
  
  #for each measurement, first ranked by themself.
  # dentify High Gini Genes with Norm.Gini
  ExprM.Stat2         = ExprM.Stat2[rev(order(ExprM.Stat2$Norm.Gini2)),]
  Genelist.HighNormGini = rownames(ExprM.Stat2[ExprM.Stat2$Norm.Gini2 > 1,])  # cutoff approach, use a cutoff to choose top genes. 
  
  # identify High Gini Genes with pvalue
  ExprM.Stat2         = ExprM.Stat2[order(ExprM.Stat2$Gini.pvalue),]
  Genelist.top_pvalue = rownames(ExprM.Stat2[ExprM.Stat2$Gini.pvalue < 0.05 & ExprM.Stat2$Norm.Gini2 > 0,])
  
  return(Genelist.top_pvalue)
}




####################################################################
#                                 FEED                             #
####################################################################
#' @description calculate NMI and ARI
#' @param truelabel reference cell labels
#' @param predlabel cluster cell labels
evaluation <- function(truelabel, predlabel) {
  if (length(truelabel) != length(predlabel))
    stop("truelabel and predlabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(predlabel)
  #Mutual information
  MI = 0.0
  for (idx in x_ids) {
    for (idy in y_ids) {
      idxOccur = which(truelabel == idx)
      idyOccur = which(predlabel == idy)
      idxyOccur = intersect(idxOccur, idyOccur)
      if (length(idxyOccur) > 0) {
        MI = MI + (length(idxyOccur) / total) * log2((length(idxyOccur) * total) /
                                                       (length(idxOccur) * length(idyOccur)))
      }
    }
  }
  #Normalized Mutual information
  Hx = 0
  #Entropies
  for (idx in x_ids) {
    idxOccurCount = length(which(truelabel == idx))
    
    Hx = Hx - (idxOccurCount / total) * log2(idxOccurCount / total)
  }
  Hy = 0
  #Entropies
  for (idy in y_ids) {
    idyOccurCount = length(which(predlabel == idy))
    Hy = Hy - (idyOccurCount / total) * log2(idyOccurCount / total)
  }
  nmi = 2 * MI / (Hx + Hy)
  
  #(adjusted) Rand Index
  tab = table(truelabel, predlabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab ^ 2) - (sum(ni ^ 2) + sum(nj ^ 2)) / 2) / n2
  ari = c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2) / n2) / ((nis2 + njs2) /
                                                                  2 - (nis2 * njs2) / n2)
  out = c(nmi,ari)
  names(out)=c("NMI","ARI")
  return(out)
}

#' @description load csv file and pre-processing it
#' @param file_path csv file path
#' @return gene expression matrix(genes x cells)
load_data <- function(file_path) {
  data_set <- read.csv(file = file_path, header = T)[, -1]
  gnames <- read.csv(file = file_path, header = T)[, 1]
  gnames <- gsub(pattern = "\\|", replacement = "-", gnames)
  gnames <- gsub(pattern = "\\_", replacement = "", gnames)
  rownames(data_set) <- gnames
  return(data_set)
}

#' @description scRNA-seq data processing
#' @param data genes in rows and samples in column
#' @param log2 log2-transform,default is TRUE
#' @return pre-processed data
rnaseq_processing <- function(data, log2 = T) {
  
  pos <-
    which(apply(data, 1, function(x) {
      length(which(x == 0))
    }) < 0.94 * ncol(data))
  data <- data[pos, ]
  
  pos <-
    which(apply(data, 1, function(x) {
      length(which(x == 0))
    }) > 0.06 * ncol(data))
  data <- data[pos, ]
  
  if (log2) {
    data <-  log2(data + 1)
  }
  return(data)
}

#' @description group genes into bins
#' @param dataset gene expression matrix(genes x cells)
#' @param num the size of bins, default 20
#' @return a list including binned genes by mean expression
bin_gene <- function(data_set, num = 20) {
  #the counts of bins
  num_bins = min(num, (nrow(data_set) - 1))
  
  gene_mean <- Matrix::rowMeans(data_set)
  gene_cut <-
    cut(x = gene_mean,
        breaks = num_bins,
        labels = c(1:num_bins))
  names(gene_cut) <- names(gene_mean)
  gene_bins <- list()
  for (i in c(1:num_bins)) {
    gene_bins[[i]] <- names(gene_cut[which(gene_cut == i)])
  }
  return(gene_bins)
}

#' @description gene decomposition
#' @param M expression of single gene
#' @param pvalue 
#' @return components and the corresponding cell labels
NORM_CLUST = function(M, pvalue = 0.05) {
  p = 0
  start = 1
  DFs = tibble()
  DF = tibble(M)
  
  DF['label'] = 1:length(M)
  DF['C'] = 0
  
  NROW = nrow(DF)
  BREAK = F
  while (NROW >= 3) {
    G = 1
    p = 0
    while (p < pvalue) {
      #GMM
      Model = Mclust(DF$M, G = G)
      if (is.null(Model)) {
        BREAK = T
        break
      }
      #select class which is more than 3 samples
      CLASS = table(Model$classification)[table(Model$classification) > 3]
      #generate a zero vector whose length is length of CLASS
      ps = numeric(length = length(CLASS))
      #calculate p-value of every class(ie. every component)
      for (cc in 1:length(CLASS)) {
        ps[cc] = ifelse(length(unique(DF$M[Model$classification ==
                                             as.numeric(names(CLASS))[cc]])) <= 1,
                        0,
                        shapiro.test(DF$M[Model$classification ==
                                            as.numeric(names(CLASS))[cc]])$p.value)
      }
      #find max p-value
      p = ps[which(ps == max(ps))]
      #the component corresponding to the max p-value
      Class = names(CLASS)[which(ps == max(ps))]
      G = G + 1
      
      #the number of components should be less than number of samples
      if (G >= length(unique(DF$M))) {
        BREAK = T
        break
      }
    }
    if (isTRUE(BREAK)) {
      DF$C = -1
      DFs = rbind(DFs, DF %>% filter(C != 0))
      DF = DF %>% filter(C == 0)
      NROW = nrow(DF)
    } else{
      #allocate labels for current component using parameter 'start'
      DF$C[Model$classification == Class] = start
      #get filtered DF
      DFs = rbind(DFs, DF %>% filter(C != 0))
      DF = DF %>% filter(C == 0)
      NROW = nrow(DF)
      start = start + 1
    }
  }
  if (nrow(DF) > 0) {
    DF$C = -1
    DFs = rbind(DFs, DF)
  }
  return(DFs %>% arrange(label) %>% .$C)
}

#' @description 
#' @param data_set pre-processed gene matrix expression matrix
#' @param gene_bin binned genes
#' @return a list of MU and SIGMA that are the mean and standard variance of gene
cal_MU_SIGMA_parallel <- function(data_set,
                                  gene_bin) {
  n_cell <- ncol(data_set)
  n_gene <- nrow(data_set)
  MU <- list()
  SIGMA <- list()
  
  if (length(gene_bin) != 0) {
    binned_data <- data_set[gene_bin, ]
    try({
      tmp_DF <- apply(binned_data,
                      1,
                      function(x) {
                        #gene decomposition
                        DF <- tibble(M = x)
                        DF = DF %>% mutate(C = NORM_CLUST(M) %>% factor())
                      })
    })
    for (l in c(1:length(tmp_DF))) {
      gene_name <- names(tmp_DF[l])
      gene_exp <-
        c(t(data_set[which(rownames(data_set) == gene_name),]))
      
      predict_label <- tmp_DF[[l]]$C
      res <- predict_label[predict_label != -1]
      # remove gene which decompose null or only one
      distri_num <- length(unique(res))
      if (distri_num == 0 | distri_num == 1) {
        next
      } else{
        for (n in c(1:distri_num)) {
          mu <- mean(gene_exp[which(predict_label == n)])
          sigma <- sd(gene_exp[which(predict_label == n)])
          MU[[gene_name]][n] <- mu
          SIGMA[[gene_name]][n] <- sigma
        }
      }
    }
  }
  return(list(MU, SIGMA))
}

#' @description equation of JS correlation
JS_COR <- function(mu1,
                   mu2,
                   sigma1,
                   sigma2) {
  m_mu <- (mu1 + mu2) / 2
  m_sigma <- (sqrt(sigma1 ^ 2 + sigma2 ^ 2)) / 2
  
  # KL_pm <- log(m_sigma/sigma1)+(sigma1^2+(mu1-m_mu)^2)/(2*m_sigma^2)-0.5
  # KL_qm <- log(m_sigma/sigma2)+(sigma2^2+(mu2-m_mu)^2)/(2*m_sigma^2)-0.5
  # return(1 / 2 * KL_pm + 1 / 2 * KL_qm)
  return(0.5 * (log(m_sigma / sigma1) + log(m_sigma / sigma2) - 1) +
           (sigma1 ^ 2 + sigma2 ^ 2 + (mu1 - m_mu) ^ 2 + (mu2 - m_mu) ^
              2) / (4 * m_sigma ^ 2))
}

#' @description calculate JS within each bin
#' @param MU,SIGMA the mean and standard variance of genes
#' @return JS correlation matrix
cal_JS_parallel <- function(MU, SIGMA) {
  group_len <- length(MU)
  if (group_len == 0) {
    JS_init <- NULL
  } else if (group_len == 1) {
    JS_init <- 1
    names(JS_init) <- names(MU)
  } else{
    js <- matrix(0, nrow = group_len, ncol = group_len)
    for (index_gene1 in c(1:(group_len - 1))) {
      for (index_gene2 in c((index_gene1 + 1):group_len)) {
        for (index1 in c(1:length(MU[[index_gene1]]))) {
          for (index2 in c(1:length(MU[[index_gene2]]))) {
            js[index_gene1, index_gene2] <- js[index_gene1, index_gene2] +
              JS_COR(MU[[index_gene1]][index1],
                     MU[[index_gene2]][index2],
                     SIGMA[[index_gene1]][index1],
                     SIGMA[[index_gene2]][index2])
          }
        }
        js[index_gene1, index_gene2] <-
          js[index_gene1, index_gene2] /
          (length(MU[[index_gene1]]) *
             length(MU[[index_gene2]]))
        js[index_gene2, index_gene1] <- js[index_gene1, index_gene2]
      }
    }
    js <- 1 - ((js - min(js)) / (max(js) - min(js)))
    # js <- 1/(1+js)
    rownames(js) <- names(MU)
    colnames(js) <- names(MU)
    JS_init <- js
  }
  # if (length(JS_init) > 1) {
  #   diag(JS_init) <- 1
  # }
  return(JS_init)
}

#' @description calculate gene importance using coefficient of variant(CV)
#' @param JS_init gene-gene Jensen-Shannon divergence correlation
#' @return scores of all genes
cal_score <- function(JS_init) {
  score <- list()
  temp_genes <- list()
  for (i in c(1:length(JS_init))) {
    JS_init[[i]][is.na(JS_init[[i]])] <- 0
    
    # only retain bin containing more than 2 genes
    if (length(JS_init[[i]]) > 9) {
      sorted_cor <-  apply(JS_init[[i]], 1, sort, decreasing = TRUE)
      max_cor <-  sorted_cor[2, ]
      min_cor <-  sorted_cor[nrow(sorted_cor), ]
      cor_range <-  max_cor - 0.75 * min_cor
      co_factor <- apply(JS_init[[i]],
                         1,
                         function(x) {
                           sd(x) / mean(x)
                         })
      # score[[i]] <- cor_range
      score[[i]] <- co_factor
      names(score[[i]]) <- rownames(JS_init[[i]])
    }
  }
  score <- unlist(score)
  return(score)
}

#' @param data_set gene expression matrix
#' @return JS initial correlation matrix
FEED_first <- function(data_set) {
  message("Running bin genes...")
  gene_bins <- bin_gene(data_set = data_set, num = 20)
  message("Done!")
  
  message("Calculate MU and SIGMA...")
  i <- NULL
  n_cores <- 15
  cl <- makeCluster(n_cores)
  doParallel::registerDoParallel(cl, cores = n_cores)
  
  cms <- foreach::foreach(
    i = c(1:20),
    .packages = c("tidyverse", "mclust"),
    .export = c(
      "data_set",
      "cal_MU_SIGMA_parallel",
      "NORM_CLUST"
    )
  ) %dopar% {
    try({
      cal_MU_SIGMA_parallel(data_set,
                            gene_bins[[i]])
    })
  }
  parallel::stopCluster(cl)
  message("Done!")
  
  message("Calculate JS correlation...")
  j <- NULL
  n_cores <- 15
  cl <- makeCluster(n_cores)
  doParallel::registerDoParallel(cl, cores = n_cores)
  
  JS_init <- foreach::foreach(j = c(1:length(cms)),
                              .export = c("cal_JS_parallel", "JS_COR")) %dopar% {
                                try({
                                  cal_JS_parallel(cms[[j]][[1]],
                                                  cms[[j]][[2]])
                                })
                              }
  parallel::stopCluster(cl)
  message("Done!")
  
  return(JS_init)
}

#' main function
#' @param data_set raw gene expression matrix
#' @return final gene subset
FEED <- function(data_set) {
  message("Running data preprocessing...")
  data_set <- rnaseq_processing(data_set)
  message("Done!")
  message("")
  
  #permute data
  message("Permuting gene expression matrix...")
  permuted_data <-
    as.data.frame(apply(data_set,
                        2,
                        function(x) {
                          sample(x, size = nrow(data_set), replace = F)
                        }))
  rownames(permuted_data) <- rownames(data_set)
  colnames(permuted_data) <- colnames(data_set)
  message("Done!")
  message("")
  
  message("Construct two kinds of initial JS correlation matrix...")
  all_data <- list(data_set, permuted_data)
  t <- NULL
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  doParallel::registerDoParallel(cl, cores = n_cores)
  JS_res <- foreach(
    t = c(1:2),
    .packages = c("parallel", "doRNG", "foreach"),
    .export = c(
      "FEED_first",
      "rnaseq_processing",
      "bin_gene",
      "NORM_CLUST",
      "cal_MU_SIGMA_parallel",
      "JS_COR",
      "cal_JS_parallel"
    )
  ) %dopar% {
    try({
      FEED_first(all_data[[t]])
    })
  }
  parallel::stopCluster(cl)
  JS_init <- JS_res[[1]]
  permuted_JS_init <- JS_res[[2]]
  message("Done!")
  message("")
  
  message("calculate threshold of score...")
  score <- cal_score(JS_init)
  permuted_score <- cal_score(permuted_JS_init)
  
  # score.th <- quantile(permuted_score, probs = 0.75)
  score.th <- mean(permuted_score) + sd(permuted_score)
  # score.th <- mean(permuted_score) + 3 * sd(permuted_score)
  
  message("Done!")
  message("")
  
  message("Select final gene subset!")
  final_genes <- names(score[which(score > score.th)])
  marker_scores <- sort(score[final_genes], decreasing = T)
  
  message("End!")
  return(list(final_genes, marker_scores))
}


####################################################################
#                                 test                             #
####################################################################

ALL_D <- c(
  "Biase",
  "Chu",
  "Chung",
  'Goolam',
  "Grover",
  "Karlsson",
  'Kim',
  'Koh',
  'Kolod',
  'Kumar',
  'Pollen',
  "Ramskold",
  "Robert",
  'Yan',
  'Zhou'
)

gene_selection_methods <-  c("DUBStepR","Expr","FEAST","GiniClust","HVG","SCMarker","FEED")

resmat_Kmeans_NMI <- matrix(nrow = length(ALL_D),ncol = length(gene_selection_methods))
rownames(resmat_Kmeans_NMI) <- ALL_D
colnames(resmat_Kmeans_NMI) <- gene_selection_methods
resmat_Kmeans_ARI <- resmat_Kmeans_NMI

# resmat_HC_NMI <- matrix(nrow = length(ALL_D),ncol = length(gene_selection_methods))
# rownames(resmat_HC_NMI) <- ALL_D
# colnames(resmat_HC_NMI) <- gene_selection_methods
# resmat_HC_ARI <- resmat_HC_NMI
# 
# resmat_Clara_NMI <- matrix(nrow = length(ALL_D),ncol = length(gene_selection_methods))
# rownames(resmat_Clara_NMI) <- ALL_D
# colnames(resmat_Clara_NMI) <- gene_selection_methods
# resmat_Clara_ARI <- resmat_Clara_NMI
# 
# resmat_Louvain_NMI <- matrix(nrow = length(ALL_D),ncol = length(gene_selection_methods))
# rownames(resmat_Louvain_NMI) <- ALL_D
# colnames(resmat_Louvain_NMI) <- gene_selection_methods
# resmat_Louvain_ARI <- resmat_Louvain_NMI
# 
# resmat_GMM_NMI <- matrix(nrow = length(ALL_D),ncol = length(gene_selection_methods))
# rownames(resmat_GMM_NMI) <- ALL_D
# colnames(resmat_GMM_NMI) <- gene_selection_methods
# resmat_GMM_ARI <- resmat_GMM_NMI


for (i in c(1:length(ALL_D))) {
  D <- ALL_D[i]
  message("progressing " ,D , "...")
  
  
  load(paste0("./dataset/",D,".Rdata"))
  data_set <- load_data(paste0("./dataset/",D,".csv"))
  if(D=="Chu"){
    label <- label$label
  }else if(D %in% c("Biase","Goolam","Grover","Karlsson","Ramskold")){
    label <- meta$label
  }
  else{
    label <- label
  }
  K <- length(unique(label))
  
  
  if((D %in% c("Biase","Ramskold"))){
    ng <- c("Expr","FEAST","GiniClust","HVG","FEED")
    Expr_genes <- Expr(data_set)
    FEAST_genes <- feast(data_set,K = K)
    Gini_genes <-  GiniClust_Fitting(Gini_Filtering(data_set)$raw)
    HVG_genes <- HVG(data_set)
    FEED_genes <- FEED(data_set)[[1]]
    genes_list <- list(Expr_genes,FEAST_genes,Gini_genes,HVG_genes,FEED_genes)
    
    for(k in c(1:length(genes_list))){
      write.csv(genes_list[[k]],
                paste0("./genelist/",
                       D,"_",ng[k],"_genes.csv"))
    }
    
    # cell cluster
    for (index in c(1:length(genes_list))) {
      method <- ng[index]
      data_new <- data_set[genes_list[[index]], ]
      
      # Kmeans
      message(method, " Kmeans...")
      NMI <- c()
      ARI <- c()
      for (r in c(1:10)) {
        Kmeans_cluster_result <- kmeans(t(data_new), centers = K)$cluster
        NMI <-
          c(NMI, round(evaluation(label, Kmeans_cluster_result)[[1]], 3))
        ARI <-
          c(ARI, round(evaluation(label, Kmeans_cluster_result)[[2]], 3))
      }
      resmat_Kmeans_NMI[D, method] <- round(mean(NMI), 3)
      resmat_Kmeans_ARI[D, method] <- round(mean(ARI), 3)
      
      # # HC
      # message(method, " hierarchical clustering...")
      # HC_cluster_result <- HC_clustering(data_new,K)
      # resmat_HC_NMI[D,method] <- round(evaluation(label,HC_cluster_result)[[1]],3)
      # resmat_HC_ARI[D,method] <- round(evaluation(label,HC_cluster_result)[[2]],3)
      # 
      # # Clara
      # message(method, " Clara...")
      # Clara_cluster_result <- clara(t(data_new),k = K)$clustering
      # resmat_Clara_NMI[D,method] <- round(evaluation(label,Clara_cluster_result)[[1]],3)
      # resmat_Clara_ARI[D,method] <- round(evaluation(label,Clara_cluster_result)[[2]],3)
      # 
      # # Louvain
      # message(method, " Louvain...")
      # Louvain_cluster_result <- Louvain_clustering(data_set,genes_list[[index]])
      # resmat_Louvain_NMI[D,method] <- round(evaluation(label,Louvain_cluster_result)[[1]],3)
      # resmat_Louvain_ARI[D,method] <- round(evaluation(label,Louvain_cluster_result)[[2]],3)
      # 
      # # GMM
      # message(method, " GMM...")
      # GMM_cluster_result <- GMM_clustering(data_new)
      # resmat_GMM_NMI[D,method] <- round(evaluation(label,GMM_cluster_result)[[1]],3)
      # resmat_GMM_ARI[D,method] <- round(evaluation(label,GMM_cluster_result)[[2]],3)
    }
  }else{
    ng <- c("DUBStepR","Expr","FEAST","GiniClust","HVG","SCMarker","FEED")
    
    DUBStepR_genes <- DUBStepR(data_set)$optimal.feature.genes  
    Expr_genes <- Expr(data_set)
    FEAST_genes <- feast(data_set,K = K)
    Gini_genes <-  GiniClust_Fitting(Gini_Filtering(data_set)$raw)
    HVG_genes <- HVG(data_set)
    SCMarker_genes <- SCMarker(data_set)
    FEED_genes <- FEED(data_set)[[1]]
    genes_list <- list(DUBStepR_genes,Expr_genes,FEAST_genes,Gini_genes,HVG_genes,SCMarker_genes,FEED_genes)
    
    for(k in c(1:length(genes_list))){
      write.csv(genes_list[[k]],
                paste0("./genelist/",
                       D,"_",gene_selection_methods[k],"_genes.csv"))
    }
    
    for (index in c(1:length(gene_selection_methods))) {
      method <- gene_selection_methods[index]
      data_new <- data_set[genes_list[[index]], ]
      data_new <- na.omit(data_new)
      
      # Kmeans
      message(method, " Kmeans...")
      NMI <- c()
      ARI <- c()
      for (r in c(1:10)) {
        Kmeans_cluster_result <- kmeans(t(data_new), centers = K)$cluster
        NMI <-
          c(NMI, round(evaluation(label, Kmeans_cluster_result)[[1]], 3))
        ARI <-
          c(ARI, round(evaluation(label, Kmeans_cluster_result)[[2]], 3))
      }
      resmat_Kmeans_NMI[D, method] <- round(mean(NMI), 3)
      resmat_Kmeans_ARI[D, method] <- round(mean(ARI), 3)
      
      # # HC
      # message(method, " hierarchical clustering...")
      # if(D!="Robert"|(D!="Grover" & method!="GiniClust")){
      #   HC_cluster_result <- HC_clustering(data_new,K)
      #   resmat_HC_NMI[D,method] <- round(evaluation(label,HC_cluster_result)[[1]],3)
      #   resmat_HC_ARI[D,method] <- round(evaluation(label,HC_cluster_result)[[2]],3)
      # }
      # 
      # # Clara
      # message(method, " Clara...")
      # Clara_cluster_result <- clara(t(data_new),k = K)$clustering
      # resmat_Clara_NMI[D,method] <- round(evaluation(label,Clara_cluster_result)[[1]],3)
      # resmat_Clara_ARI[D,method] <- round(evaluation(label,Clara_cluster_result)[[2]],3)
      # 
      # # Louvain
      # message(method, " Louvain...")
      # Louvain_cluster_result <- Louvain_clustering(data_set,genes_list[[index]])
      # resmat_Louvain_NMI[D,method] <- round(evaluation(label,Louvain_cluster_result)[[1]],3)
      # resmat_Louvain_ARI[D,method] <- round(evaluation(label,Louvain_cluster_result)[[2]],3)
      # 
      # # GMM
      # message(method, " GMM...")
      # GMM_cluster_result <- GMM_clustering(data_new)
      # resmat_GMM_NMI[D,method] <- round(evaluation(label,GMM_cluster_result)[[1]],3)
      # resmat_GMM_ARI[D,method] <- round(evaluation(label,GMM_cluster_result)[[2]],3)
    }
  }
    message(D," is Done!")

}

write.csv(resmat_Kmeans_NMI,"./results/Kmeans_NMI.csv")
write.csv(resmat_Kmeans_ARI,"./results/Kmeans_ARI.csv")

# write.csv(resmat_HC_NMI,"./results/HC_NMI.csv")
# write.csv(resmat_HC_ARI,"./results/HC_ARI.csv")
# 
# write.csv(resmat_Clara_NMI,"./results/Clara_NMI.csv")
# write.csv(resmat_Clara_ARI,"./results/Clara_ARI.csv")
# 
# write.csv(resmat_Louvain_NMI,"./results/Louvain_NMI.csv")
# write.csv(resmat_Louvain_ARI,"./results/Louvain_ARI.csv")
# 
# write.csv(resmat_GMM_NMI,"./results/GMM_NMI.csv")
# write.csv(resmat_GMM_ARI,"./results/GMM_ARI.csv")






