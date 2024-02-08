source('setup.R')

cpm <- t(t(1e6*counts)/colSums(counts))
# PCA
filt <- which(rowMeans(cpm) > 1)
pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=64)



# clustering
clusMat <- NULL

for(ndim in c(5,12)){
    clusMat <- cbind(clusMat, sapply(3:8, function(k){
        kmeans(pca$x[,1:ndim], k)$cluster
    }))
    clusMat <- cbind(clusMat, sapply(c(3,5,7), function(k){
        cutree(hclust(dist(pca$x[,1:ndim])), k = k)
    }))
}


coClus <- apply(clusMat, 1, function(x){
    colMeans(x == t(clusMat))
})


clus <- cutree(hclust(as.dist(1-coClus)), k=6)


heatmap(coClus, ColSideColors = as.character(clus), RowSideColors = pheno$color)


require(plotly)
# plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = factor(clus))


pheno$cluster <- factor(clus)

