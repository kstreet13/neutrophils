source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

genelist <- unique(c(
    read.table('data/DEresults/SC_timepoint/A_vs_B/down.csv')[,1],
    read.table('data/DEresults/SC_timepoint/A_vs_B/up.csv')[,1],
    read.table('data/DEresults/SC_timepoint/A_vs_C/down.csv')[,1],
    read.table('data/DEresults/SC_timepoint/A_vs_C/up.csv')[,1],
    read.table('data/DEresults/SC_timepoint/B_vs_C/down.csv')[,1],
    read.table('data/DEresults/SC_timepoint/B_vs_C/up.csv')[,1]
))

heatdat <- as.matrix(cpm[genelist, pheno$group == 'SC'])
# force ordering
heatdat <- heatdat[,  c(paste0('SC',1:8,'A'),paste0('SC',1:8,'B'),paste0('SC',1:8,'C'))]
subpheno <- pheno[match(colnames(heatdat), pheno$sample), ]

require(NMF)
anncolors = list(
    timepoint = c(A = "green", B = "red", C = "purple"))
aheatmap(log1p(heatdat), annCol = subpheno[,"timepoint",drop=FALSE],
         Colv = NA, scale = "row",
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         filename = '~/Desktop/heat.png')
# Normalized expression values are calculated by log(CPM+1). For the heatmap, each gene is then normalized to a mean of 0 and standard deviation of 1 to show differences between timepoints.



# PCA
require(Matrix); require(matrixStats); require(sparseMatrixStats)

rv <- rowVars(cpm[, pheno$group == 'SC'])
ind <- order(rv, decreasing = TRUE)[1:5000]


pca <- BiocSingular::runPCA(t(log1p(cpm[ind, pheno$group == 'SC'])), rank = sum(pheno$group == 'SC'))

pctvar <- pca$sdev^2 / sum(pca$sdev^2)

plot(pca$x[,1:2], col=pheno$color[pheno$group=='SC'], asp=1,
     xlab = paste0('PC1 (',format(100*pctvar[1], digits = 3),'%)'),
     ylab = paste0('PC2 (',format(100*pctvar[2], digits = 3),'%)'))
# %var explained



# DE using DESeq
{
    source('setup.R')
    require(DESeq2)
    pheno$timepoint <- factor(pheno$timepoint)
    pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
    
    ind <- which(pheno$group=='SC')
    ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                                 colData = pheno[ind,],
                                 design= ~timepoint)
    res <- DESeq(ds)
    res1 <- results(res, alpha = .05, contrast = c('timepoint', 'B', 'A'))
    res2 <- results(res, alpha = .05, contrast = c('timepoint', 'C', 'A'))
    res3 <- results(res, alpha = .05, contrast = c('timepoint', 'C', 'B'))
    names(res1) <- paste0('AB_', names(res1))
    names(res2) <- paste0('AC_', names(res2))
    names(res3) <- paste0('BC_', names(res3))
    
    deseq <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseq



zscores <- (heatdat - rowMeans(heatdat)) / rowSds(heatdat)
means <- t(apply(zscores,1,function(cts){
    c(mean(cts[1:8]), mean(cts[9:16]), mean(cts[17:24]))
}))

plot(c(1,3), range(means), col='white')
for(i in 1:nrow(means)){
    lines(1:3, means[i,], col=rgb(0,0,0,.3))
}



# x: maximum (absolute) change from baseline for each gene
# y: expression at timepoint A
# 1) all genes, 2) DE genes

x <- cbind(deseq$AB_log2FoldChange, deseq$AC_log2FoldChange)
ind <- !is.na(x[,1])
x <- x[ind, ]
x <- apply(x,1,function(x){
    x[which.max(abs(x))]
})
means <- t(apply(log1p(cpm),1,function(cts){
    c(mean(cts[1:8]), mean(cts[9:16]), mean(cts[17:24]))
}))
y <- means[ind,1]


plot(x,y)




# define recovery:
# DE in A->B, not DE in A->C and B->C shift has opposite sign
#  DE/not DE vs. A-B, A-C, B-C

table(deseq$AB_padj < .05, deseq$AC_padj < .05)


# down and stay down:
# DE in A->B and in A->C, with same sign for both


