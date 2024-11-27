source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

#######
# PCA #
#######
rm <- rowMeans(cpm)
rv <- rowVars(cpm)
filt <- which(rm > 1 & rv > 1)
#filt <- which(rowVars(cpm) >= sort(rowVars(cpm), decreasing = TRUE)[10000])

pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=64)

# PC1 & PC2
plot(pca$x, asp=1, col=pheno$color, cex=3)


# 3D
require(plotly)
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$color, colors = sort(unique(pheno$color)))
