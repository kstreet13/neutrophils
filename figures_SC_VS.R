source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

genelist <- unique(c(
    read.table('data/DEresults/MVS_timepoint/A_vs_B/down.csv')[,1],
    read.table('data/DEresults/MVS_timepoint/A_vs_B/up.csv')[,1],
    read.table('data/DEresults/MVS_timepoint/A_vs_C/down.csv')[,1],
    read.table('data/DEresults/MVS_timepoint/A_vs_C/up.csv')[,1],
    read.table('data/DEresults/MVS_timepoint/B_vs_C/down.csv')[,1],
    read.table('data/DEresults/MVS_timepoint/B_vs_C/up.csv')[,1],
    
    read.table('data/DEresults/SVS_timepoint/A_vs_B/down.csv')[,1],
    read.table('data/DEresults/SVS_timepoint/A_vs_B/up.csv')[,1],
    read.table('data/DEresults/SVS_timepoint/A_vs_C/down.csv')[,1],
    read.table('data/DEresults/SVS_timepoint/A_vs_C/up.csv')[,1],
    read.table('data/DEresults/SVS_timepoint/B_vs_C/down.csv')[,1],
    read.table('data/DEresults/SVS_timepoint/B_vs_C/up.csv')[,1]
))

heatdat <- as.matrix(cpm[genelist, pheno$group %in% c('SC','MVS','SVS')])
# force ordering
heatdat <- heatdat[,  c(paste0('SC',1:8,'A'), paste0('SC',1:8,'B'), paste0('SC',1:8,'C'),
                        paste0('MVS',1:3,'A'),paste0('MVS',1:3,'B'),paste0('MVS',1:3,'C'),
                        paste0('SVS',1:4,'A'),paste0('SVS',1:4,'B'),paste0('SVS',1:4,'C'))]
subpheno <- pheno[match(colnames(heatdat), pheno$sample), ]
subpheno$annTime <- subpheno$timepoint
subpheno$annTime[subpheno$group != 'SC'] <- paste0(subpheno$annTime[subpheno$group != 'SC'],'_VS')

require(NMF)
anncolors = list(
    annTime = c(A = "green", B = "red", C = "purple",
                  A_VS = "forestgreen", B_VS = "firebrick", C_VS = "slateblue4"))
aheatmap(log1p(heatdat), 
         annCol = subpheno[,"annTime",drop=FALSE],
         Colv = NA, scale = "row",
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         filename = '~/Desktop/heat_SC_VS.png')
# Expression values are normalized by dividing each count by the total for that sample and multiplying by 1 million (Counts Per Million, CPM), then adding a pseudocount of 1 and taking the natural log, to control outliers. For the heatmap, each gene's normalized expression is then standardized to a mean of 0 and standard deviation of 1 to show differences between timepoints.



# PCA
require(Matrix); require(matrixStats); require(sparseMatrixStats)

rv <- rowVars(cpm[, pheno$group %in% c('SC','MVS','SVS')])
ind <- order(rv, decreasing = TRUE)[1:5000]


pca <- BiocSingular::runPCA(t(log1p(cpm[ind, pheno$group %in% c('SC','MVS','SVS')])), rank = sum(pheno$group %in% c('SC','MVS','SVS')))

pctvar <- pca$sdev^2 / sum(pca$sdev^2)

plot(pca$x[,1:2], col=pheno$color[pheno$group %in% c('SC','MVS','SVS')], asp=1,
     xlab = paste0('PC1 (',format(100*pctvar[1], digits = 3),'%)'),
     ylab = paste0('PC2 (',format(100*pctvar[2], digits = 3),'%)'), 
     pch = c(16,17)[1+(pheno$group=='SVS')[pheno$group %in% c('SC','MVS','SVS')]])
# %var explained



require(rgl)
plot3d(pca$x[,1:3], col = pheno$color[pheno$group %in% c('SC','MVS','SVS')],
       aspect = 'iso', size = 10)









###########
# FIGURES #
###########

# Figure 2C: volcano plot max fold change vs baseline expression, all genes
# Dataset: all SVS samples

x <- cbind(deseqSVS$AB_log2FoldChange, deseqSVS$AC_log2FoldChange)
ind <- ifelse(abs(x[,2]) > abs(x[,1]), 2, 1)
ind[is.na(ind)] <- 1
x <- sapply(1:nrow(x),function(i){
    x[i, ind[i]]
})
means <- cbind(
    rowMeans(log1p(cpm[,paste0('SVS',1:4,'A')])),
    rowMeans(log1p(cpm[,paste0('SVS',1:4,'B')])),
    rowMeans(log1p(cpm[,paste0('SVS',1:4,'C')]))
)
y <- means[,1]

# all genes
plot(x,y, col=rgb(0,0,0,.5),
     xlab = 'Max. (abs) fold change from baseline',
     ylab = 'Avg. normalized expression at baseline',
     main = 'All genes - SVS')
abline(v=0,lty=2,col=2)


# Figure 2D: volcano plot max fold change vs baseline expression, DE genes only
# Dataset: all SVS samples

# de genes
DEind <- which(deseqSVS$AB_padj < .05 | deseqSVS$AC_padj < .05 | deseqSVS$BC_padj < .05)
plot(x[DEind], y[DEind], col=rgb(0,0,0,.5),
     xlab = 'Max. (abs) fold change from baseline',
     ylab = 'Avg. normalized expression at baseline',
     main = 'DE genes - SVS')
abline(v=0,lty=2,col=2)


