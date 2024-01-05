source('setup.R')

cpm <- t(t(1e6*counts)/colSums(counts))

boxplot(as.matrix(log1p(cpm)), col = pheno$color)

require(EDASeq)

plotRLE(as.matrix(counts), col = pheno$color)
plotRLE(as.matrix(cpm), col = pheno$color)

# quick PCA
cpm <- cpm[which(rowMeans(cpm) > 1), ]

pca <- BiocSingular::runPCA(t(cpm), rank=64)

plot(pca$x, asp=1, col=pheno$color)

