source('setup.R')

cpm <- t(t(1e6*counts)/colSums(counts))

boxplot(as.matrix(log1p(cpm)), col = pheno$color)

require(EDASeq)

# plotRLE(as.matrix(counts), col = pheno$color)
# plotRLE(as.matrix(cpm), col = pheno$color)

# quick PCA
#filt <- which(rowMeans(cpm) > 1)
filt <- which(rowVars(cpm) >= sort(rowVars(cpm), decreasing = TRUE)[5000])

pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=64)

plot(pca$x, asp=1, col=pheno$color, cex=3)


plot(pca$x, asp=1, col='white')
text(pca$x, labels = pheno$patient, col=pheno$color)


require(plotly)
cc <- pheno$color
cc[pheno$mmp8 == 'MMP8+'] <- 'darkorange3'
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = cc, colors = sort(unique(cc)))


# PCA on subsets:
# 1. Surgical controls vs healthy controls (all timepoints)
ind <- which(pheno$group %in% c('SC','EC'))

subcpm <- cpm[,ind]
filt <- which(rowMeans(subcpm) > 1)
pca <- BiocSingular::runPCA(t(log1p(subcpm[filt, ])), rank=length(ind))
plot(pca$sdev^2/sum(pca$sdev^2))
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$color[ind], colors = sort(unique(pheno$color[ind])))



# 2. Post-bypass timepoints (B and C) all patients ?
ind <- which(pheno$timepoint %in% c('B','C'))

subcpm <- cpm[,ind]
filt <- which(rowMeans(subcpm) > 1)
pca <- BiocSingular::runPCA(t(log1p(subcpm[filt, ])), rank=length(ind))
plot(pca$sdev^2/sum(pca$sdev^2))
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$color[ind], colors = sort(unique(pheno$color[ind])))


