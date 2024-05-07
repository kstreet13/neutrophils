# dendrograms a few ways

source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))
filt <- which(rowMeans(cpm) > 1)

# 1: Correlation-based distance
corr <- cor(as.matrix(cpm[filt, ]))

require(ape)
dd1 <- as.dist(1-corr) 
hc1 <- hclust(dd1, method = "ward.D2") 

plot(as.phylo(hc1), tip.color = pheno$color, font = 2, cex = 0.5,
     main = 'Correlation dist, CPM > 1 (10,836 genes)')

clus1 <- factor(cutree(hc1, k=5))
plot(as.phylo(hc1), type = "unrooted", tip.color = colorby(clus1), font = 2, cex = 0.5,
     main = 'Correlation dist')


# 2: PCA-based distance
pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=64)

require(ape)
dd2 <- dist(pca$x[,1:12]) # 12 looks about right 
hc2 <- hclust(dd2, method = "ward.D2") 

plot(as.phylo(hc2), type = "unrooted", tip.color = pheno$color, font = 2, cex = 0.5,
     main = 'PCA dist')

clus2 <- factor(cutree(hc2, k=4))
plot(as.phylo(hc2), type = "unrooted", tip.color = colorby(clus2), font = 2, cex = 0.5,
     main = 'PCA dist')


table(clus1, clus2)


# 3: # use most variable genes in surgical samples (SC, SVS, MVS)
rv <- rowVars(as.matrix(cpm[, pheno$group %in% c('SC','SVS','MVS')]))
filt <- which(rv >= sort(rv,decreasing = TRUE)[5000])

corr <- cor(as.matrix(cpm[filt, ]))

require(ape)
dd3 <- as.dist(1-corr) 
hc3 <- hclust(dd3, method = "ward.D2") 

plot(as.phylo(hc3), tip.color = pheno$color, font = 2, cex = 0.5,
     main = 'Correlation dist, Top 5,000 most variable in surgical')

