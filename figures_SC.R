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
aheatmap(log1p(heatdat), 
         annCol = subpheno[,"timepoint",drop=FALSE],
         Colv = NA, scale = "row",
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         filename = '~/Desktop/heat_SC.png')
# Expression values are normalized by dividing each count by the total for that sample and multiplying by 1 million (Counts Per Million, CPM), then adding a pseudocount of 1 and taking the natural log, to control outliers. For the heatmap, each gene's normalized expression is then standardized to a mean of 0 and standard deviation of 1 to show differences between timepoints.



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
ind <- ifelse(abs(x[,2]) > abs(x[,1]), 2, 1)
ind[is.na(ind)] <- 1
x <- sapply(1:nrow(x),function(i){
    x[i, ind[i]]
})

means <- t(apply(log1p(cpm),1,function(cts){
    c(mean(cts[1:8]), mean(cts[9:16]), mean(cts[17:24]))
}))
y <- means[,1]

# all genes
plot(x,y, col=rgb(0,0,0,.5),
     xlab = 'Max. (abs) fold change from baseline',
     ylab = 'Avg. normalized expression at baseline',
     main = 'All genes')
abline(v=0,lty=2,col=2)

# de genes
DEind <- which(deseq$AB_padj < .05 | deseq$AC_padj < .05 | deseq$BC_padj < .05)
plot(x[DEind], y[DEind], col=rgb(0,0,0,.5),
     xlab = 'Max. (abs) fold change from baseline',
     ylab = 'Avg. normalized expression at baseline',
     main = 'DE genes')
abline(v=0,lty=2,col=2)

hist(y[x < 0])
hist(y[x > 0])

############
# PATTERNS #
############
# define recovery:
# DE in A->B, not DE in A->C and B->C shift has opposite sign
#recov.ind <- which(deseq$AB_padj < .05 & deseq$AC_padj > .05 & (sign(deseq$AB_log2FoldChange) == -sign(deseq$BC_log2FoldChange)))
recov.ind <- which(deseq$AB_padj < .05 & deseq$BC_padj < .05 & deseq$AC_padj > .05)

# down early/no recovery:
# DE in A->B and in A->C, with same sign for both
#down.ind <- which(deseq$AB_padj < .05 & deseq$AC_padj < .05 & deseq$AB_log2FoldChange < 0  & deseq$AC_log2FoldChange < 0)
#norecov.ind <- which(deseq$AB_padj < .05 & deseq$AC_padj < .05 & (sign(deseq$AB_log2FoldChange) == sign(deseq$AC_log2FoldChange)))
early.ind <- which(deseq$AB_padj < .05 & deseq$BC_padj > .05)

# late induction:
# not DE in A->B, DE in B->C and A->C
#late.ind <- which(deseq$AB_padj > .05 & deseq$AC_padj < .05 & deseq$BC_padj < .05)
#late.ind <- which(deseq$AC_padj < .05 & deseq$BC_padj < .05 & (sign(deseq$AC_log2FoldChange) == sign(deseq$BC_log2FoldChange)))
late.ind <- which(deseq$AB_padj > .05 & deseq$BC_padj < .05)

# continuous change:
# DE in A->B, DE in B->C, same sign
cont.ind <- which(deseq$AB_padj < .05 & deseq$BC_padj < .05 & (sign(deseq$AB_log2FoldChange)==sign(deseq$BC_log2FoldChange)))


cc <- brewer.pal(5,'Set1')[c(1,3,2,5)]
catcounts <- c(Early = length(early.ind),
               Recovery = length(recov.ind),
               Late = length(late.ind),
               Progressive = length(cont.ind))

barplot(catcounts,
        col = cc, ylab = '# of Genes', xlab = 'Pattern')


sum(deseq$AB_padj < .05, na.rm=TRUE)
sum(deseq$AB_padj < .05 & deseq$AC_padj > .05 & (sign(deseq$AB_log2FoldChange) == -sign(deseq$BC_log2FoldChange)), na.rm = TRUE)
barplot(c(Recovery = 1924, NoRecovery = 4225-1924))

# explainer figure
plot(c(0,4), c(-2,2), col='white')
abline(h=0, lty=2)
lines(c(0,1,0), col=cc[2], lwd=3)
lines(c(0,-1,0), col=cc[2], lwd=3)
lines((1:3)-.05, c(0,1,1), col=cc[1], lwd=3)
lines((1:3)-.05, c(0,-1,-1), col=cc[1], lwd=3)
lines((1:3), c(0,0,1), col=cc[3], lwd=3)
lines((1:3), c(0,0,-1), col=cc[3], lwd=3)
lines((1:3)-.1, c(0,1,2), col=cc[4], lwd=3)
lines((1:3)-.1, c(0,-1,-2), col=cc[4], lwd=3)



pattern <- rep(NA, nrow(heatdat))
pattern[rownames(heatdat) %in% rownames(deseq)[recov.ind]] <- 'Recovery'
pattern[rownames(heatdat) %in% rownames(deseq)[early.ind]] <- 'Early'
pattern[rownames(heatdat) %in% rownames(deseq)[late.ind]] <- 'Late'
pattern[rownames(heatdat) %in% rownames(deseq)[cont.ind]] <- 'Progressive'
genemeta <- data.frame(Pattern = factor(pattern))

anncolors = list(
    timepoint = c(A = "green", B = "red", C = "purple"),
    Pattern = c(Early = cc[1], Recovery = cc[2], Late = cc[3], Progressive = cc[4]))
aheatmap(log1p(heatdat), 
         annCol = subpheno[,"timepoint",drop=FALSE],
         annRow = genemeta,
         Colv = NA, scale = "row",
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         filename = '~/Desktop/heat_SC.png')




# the B->C change tends to reverse the A->B change
cc <- rep('grey50', nrow(deseq))
cc[recov.ind] <- 3
cc[early.ind] <- 2
cc[late.ind] <- 4

# add counts of different groups
plot(range(-log10(deseq$AB_padj)*sign(deseq$AB_log2FoldChange), na.rm = TRUE),
     range(-log10(deseq$BC_padj)*sign(deseq$BC_log2FoldChange), na.rm = TRUE),
     col = 'white', asp = 1,
     xlab = 'Signed -log10 p-value for A vs. B',
     ylab = 'Signed -log10 p-value for B vs. C')
abline(h=0,lty=2);abline(v=0,lty=2)
points(-log10(deseq$AB_padj)*sign(deseq$AB_log2FoldChange), -log10(deseq$BC_padj)*sign(deseq$BC_log2FoldChange), col = alpha(cc, alpha=.5), pch=16)
legend('topright', col = c(2,4,3), pch = 16, bty = 'n', cex = 1,
       legend = c(paste0('Early change (n=',length(early.ind),')'),
                  paste0('Late change (n=',length(late.ind),')'),
                  paste0('Recovery (n=',length(recov.ind),')')))


plot(deseq$AB_log2FoldChange, deseq$BC_log2FoldChange, col = alpha(cc, alpha=.5), asp = 1)
abline(h=0,col=2);abline(v=0,col=2)






#############
# # LISTS # #
#############
early.genes <- rownames(deseq)[early.ind]
recov.genes <- rownames(deseq)[recov.ind]
late.genes <- rownames(deseq)[late.ind]
cont.genes <- rownames(deseq)[cont.ind]



write.table(early.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/earlySC.txt')
write.table(recov.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/recovSC.txt')
write.table(late.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/lateSC.txt')
write.table(cont.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/contSC.txt')
