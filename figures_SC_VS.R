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



# DE using DESeq
{
    source('setup.R')
    require(DESeq2)
    pheno$timepoint <- factor(pheno$timepoint)
    pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
    
    ind <- which(pheno$group=='MVS')
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
    
    deseqMVS <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseqMVS
{
    source('setup.R')
    require(DESeq2)
    pheno$timepoint <- factor(pheno$timepoint)
    pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
    
    ind <- which(pheno$group=='SVS')
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
    
    deseqSVS <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseqSVS


############
# PATTERNS #
############
# define recovery:
# DE in A->B, not DE in A->C and B->C shift has opposite sign
recovMVS <- which(deseqMVS$AB_padj < .05 & deseqMVS$BC_padj < .05 & deseqMVS$AC_padj > .05)
recovSVS <- which(deseqSVS$AB_padj < .05 & deseqSVS$BC_padj < .05 & deseqSVS$AC_padj > .05)

# early change/no recovery:
# DE in A->B, not in B->C
earlyMVS <- which(deseqMVS$AB_padj < .05 & deseqMVS$BC_padj > .05 & deseqMVS$AC_padj < .05)
earlySVS <- which(deseqSVS$AB_padj < .05 & deseqSVS$BC_padj > .05 & deseqSVS$AC_padj < .05)

# late induction:
# not DE in A->B, DE in B->C
lateMVS <- which(deseqMVS$AB_padj > .05 & deseqMVS$BC_padj < .05 & deseqMVS$AC_padj < .05)
lateSVS <- which(deseqSVS$AB_padj > .05 & deseqSVS$BC_padj < .05 & deseqSVS$AC_padj < .05)

# continuous change:
# DE in A->B, DE in B->C, same sign
contMVS <- which(deseqMVS$AB_padj < .05 & deseqMVS$BC_padj < .05 & (sign(deseqMVS$AB_log2FoldChange)==sign(deseqMVS$BC_log2FoldChange)))
contSVS <- which(deseqSVS$AB_padj < .05 & deseqSVS$BC_padj < .05 & (sign(deseqSVS$AB_log2FoldChange)==sign(deseqSVS$BC_log2FoldChange)))


catcounts <- matrix(c(
    length(earlyMVS),length(earlySVS),
    length(recovMVS),length(recovSVS),
    length(lateMVS), length(lateSVS),
    length(contMVS), length(contSVS)
), nrow=2)
rownames(catcounts) <- c('MVS','SVS')
colnames(catcounts) <- c('Early','Recovery','Late','Progressive')

# barplot number of genes in each category
cc <- brewer.pal(5,'Set1')[c(1,3,2,5)]


barplot(catcounts, beside=TRUE, col = rep(cc,each=2),
        ylab = '# of Genes', xlab = 'Pattern')
barplot(catcounts, beside=TRUE, add = TRUE, 
        col = c(rgb(1,1,1,.3), rgb(0,0,0,.3)))

barplot(c(length(intersect(earlyMVS,earlySVS)),
          length(intersect(recovMVS,recovSVS)),
          length(intersect(lateMVS,lateSVS)),
          length(intersect(contMVS,contSVS))),
        add=TRUE, space=.5, width=2, col=rgb(0,0,0,0))

legend('topright', legend = c('MVS','SVS'), bty='n', fill=c('grey80','grey20'))




#############
# # LISTS # #
#############
earlyMVSgenes <- rownames(deseqMVS)[earlyMVS]
recovMVSgenes <- rownames(deseqMVS)[recovMVS]
lateMVSgenes <- rownames(deseqMVS)[lateMVS]
contMVSgenes <- rownames(deseqMVS)[contMVS]

earlySVSgenes <- rownames(deseqSVS)[earlySVS]
recovSVSgenes <- rownames(deseqSVS)[recovSVS]
lateSVSgenes <- rownames(deseqSVS)[lateSVS]
contSVSgenes <- rownames(deseqSVS)[contSVS]


write.table(earlyMVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/earlyMVS.txt')
write.table(recovMVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/recovMVS.txt')
write.table(lateMVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/lateMVS.txt')
write.table(contMVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/contMVS.txt')


write.table(earlySVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/earlySVS.txt')
write.table(recovSVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/recovSVS.txt')
write.table(lateSVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/lateSVS.txt')
write.table(contSVSgenes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/contSVS.txt')


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


# Figure 2E: WOULD LOVE TO HAVE A LINE GRAPH WITH 4 DIFFERENT COLORS SHOWING THE BEHAVIOR OF ALL GENES WITHIN EACH DYNAMIC CATEGORY (AS DEFINED IN THE SLIDE FOLLOWING FIGURE 2).
# Dataset: all SVS samples

# messy line plot (with colors!)
genelist <- c(earlySVSgenes, recovSVSgenes, lateSVSgenes, contSVSgenes)
cc <- rep(brewer.pal(5,'Set1')[c(1,3,2,5)], 
          times = c(length(earlySVSgenes), length(recovSVSgenes), length(lateSVSgenes), length(contSVSgenes)))
linedat <- as.matrix(cpm[genelist, pheno$group == 'SVS'])
linedat <- linedat[,  c(paste0('SVS',1:4,'A'),paste0('SVS',1:4,'B'),paste0('SVS',1:4,'C'))]

zscores <- (linedat - rowMeans(linedat)) / rowSds(linedat)
means <- cbind(
    rowMeans(zscores[,paste0('SVS',1:4,'A')]),
    rowMeans(zscores[,paste0('SVS',1:4,'B')]),
    rowMeans(zscores[,paste0('SVS',1:4,'C')])
)
ind <- sample(nrow(linedat))

plot(c(1,3), range(means), col='white')
for(i in 1:nrow(means)){
    lines(1:3, means[ind[i],], col=alpha(cc[ind[i]], alpha=.3))
}


# separated by colors
cc <- brewer.pal(5,'Set1')[c(1,3,2,5)]
layout(matrix(1:16,nrow=4,byrow = TRUE))
par(mar=c(.2,.2,.2,.2))
# early
replicate(5, plot.new())
plot(c(1,3), range(means), col='white', axes=FALSE); box(); axis(2)
for(g in earlySVSgenes){
    lines(1:3, means[g,], col=alpha(cc[1], alpha=.2))
}
# recov
plot(c(1,3), range(means), col='white', axes=FALSE); box()
for(g in recovSVSgenes){
    lines(1:3, means[g,], col=alpha(cc[2], alpha=.2))
}
plot.new()
legend('bottomleft',bty='n', col=cc, lwd=2, legend = c('Early','Recovery','Late','Continuous'), title = 'Category')
plot.new()
# late
plot(c(1,3), range(means), col='white', axes=FALSE); box(); axis(1, at=1:3, labels=LETTERS[1:3]); axis(2)
for(g in lateSVSgenes){
    lines(1:3, means[g,], col=alpha(cc[3], alpha=.2))
}
# cont
plot(c(1,3), range(means), col='white', axes=FALSE); box(); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in contSVSgenes){
    lines(1:3, means[g,], col=alpha(cc[4], alpha=.2))
}
replicate(5, plot.new())
layout(1)
par(mar=c(5,4,4,2)+.1)

