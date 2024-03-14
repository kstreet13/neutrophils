source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

# Cluster based on expression of all expressed genes only the surgical control
# samples (SC1A, SC1B, SC1C through SC8A, SC8B, SC8C samples, so 24 samples in
# all) to see if samples from the same time point cluster together. I expect
# this to form 3 distinct clusters based on our prior work.  If you can include
# a description of how many genes are expressed in these samples and the
# criteria you use to decide which ones count as "expressed" in a brief blurb,
# that will hopefully help me ask fewer dumb questions when I revisit the
# analysis 😉

idx <- which(pheno$group == 'SC')
pheno <- pheno[idx, ]
counts <- counts[, idx]
cpm <- cpm[, idx]

# PCA
filt <- which(rowMeans(cpm) > 1)
pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=ncol(counts))
plot(c(0,cumsum(pca$sdev^2))/sum(pca$sdev^2), type = 'b')


# clustering
clusMat <- NULL

require(mclust)
mc <- Mclust(pca$x[,1:4], G=1:6)

for(ndim in c(2,4,7)){
    clusMat <- cbind(clusMat, sapply(3:6, function(k){
        kmeans(pca$x[,1:ndim], k)$cluster
    }))
    clusMat <- cbind(clusMat, sapply(c(3,5,7), function(k){
        cutree(hclust(dist(pca$x[,1:ndim])), k = k)
    }))
}


coClus <- apply(clusMat, 1, function(x){
    colMeans(x == t(clusMat))
})


clus <- cutree(hclust(as.dist(1-coClus)), k=3)


heatmap(coClus, ColSideColors = as.character(clus), RowSideColors = pheno$color)
heatmap(as.matrix(dist(pca$x[,1:4])), ColSideColors = as.character((c(4,6,7))[mc$classification]), RowSideColors = pheno$color, symm = TRUE)



require(plotly)
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = factor(clus))
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$color, colors = sort(unique(pheno$color))) |> layout(scene = list(aspectmode='data'))




pheno$cluster <- factor(clus)



# Cluster based on expression of all expressed genes of only the sepsis samples
# (SS1-12, twelve samples total).  Based on our prior analysis, I expect this
# will be kind of a mess.  We saw most MMP8+ samples (SS1-6) in one cluster, but
# not all of them--and the MMP8- clusters were kind of all over the place.  Only
# when we added the surgical samples to the mix did we see any kind of a
# pattern.  But interested to see what we find with better methods.

source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

idx <- which(pheno$group == 'SS')
pheno <- pheno[idx, ]
counts <- counts[, idx]
cpm <- cpm[, idx]

# PCA
filt <- which(rowMeans(cpm) > 1)
pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=ncol(counts))
plot(c(0,cumsum(pca$sdev^2))/sum(pca$sdev^2), type = 'b')
plot((pca$sdev^2)/sum(pca$sdev^2), type = 'b')




# clustering
clusMat <- NULL

require(mclust)
mc <- Mclust(pca$x[,1:4], G=1:6)

for(ndim in c(2,4,6)){
    clusMat <- cbind(clusMat, sapply(2:4, function(k){
        kmeans(pca$x[,1:ndim], k)$cluster
    }))
    clusMat <- cbind(clusMat, sapply(2:4, function(k){
        cutree(hclust(dist(pca$x[,1:ndim])), k = k)
    }))
}


coClus <- apply(clusMat, 1, function(x){
    colMeans(x == t(clusMat))
})


clus <- cutree(hclust(as.dist(1-coClus)), k=3)


heatmap(coClus, ColSideColors = as.character(clus), RowSideColors = pheno$mmp8col, labRow = pheno$mmp8)
heatmap(as.matrix(dist(pca$x[,1:4])), ColSideColors = as.character((c(4,6,7))[mc$classification]), RowSideColors = pheno$mmp8col, symm = TRUE, labRow = pheno$mmp8)


require(plotly)
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = factor(clus))
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$mmp8col, colors = sort(unique(pheno$mmp8col))) |> layout(scene = list(aspectmode='data'))




