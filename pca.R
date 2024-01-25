source('setup.R')

cpm <- t(t(1e6*counts)/colSums(counts))

boxplot(as.matrix(log1p(cpm)), col = pheno$color)

require(EDASeq)

plotRLE(as.matrix(counts), col = pheno$color)
plotRLE(as.matrix(cpm), col = pheno$color)

# quick PCA
cpm <- cpm[which(rowMeans(cpm) > 5), ]

pca <- BiocSingular::runPCA(t(log1p(cpm)), rank=64)

plot(pca$x, asp=1, col=pheno$color, cex=3)


plot(pca$x, asp=1, col='white')
text(pca$x, labels = pheno$patient, col=pheno$color)


require(plotly)


plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$color, colors = sort(unique(pheno$color)))




