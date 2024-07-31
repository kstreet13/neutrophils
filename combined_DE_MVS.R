# combine DE methods (1) DESeq2 and (2) at least 1.5 (not log) fold change in all samples ("consistent fold change" or cfc)

# MVS 

# for comparison
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
    
    deseq <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseq

source('setup.R')

counts <- counts[, pheno$group == 'MVS']

cpm <- t(t(1e6*counts)/colSums(counts))

pheno <- pheno[pheno$group == 'MVS', ]
# check order
all(pheno$sample == c("MVS1A","MVS1B","MVS1C","MVS2A","MVS2B","MVS2C","MVS3A","MVS3B","MVS3C"))
all(colnames(counts) == c("MVS1A","MVS1B","MVS1C","MVS2A","MVS2B","MVS2C","MVS3A","MVS3B","MVS3C"))

A <- cpm[,grep('A$', pheno$sample)]
B <- cpm[,grep('B$', pheno$sample)]
C <- cpm[,grep('C$', pheno$sample)]

AB <- B/A; AB[is.nan(AB)] <- 1 # 1/0 = Inf, 0/0 = NaN
AC <- C/A; AC[is.nan(AC)] <- 1
BC <- C/B; BC[is.nan(BC)] <- 1


require(matrixStats)
cfc <- data.frame(
    ABup = rowMins(as.matrix(AB)) > 1.5,
    ABdn = rowMaxs(as.matrix(AB)) < 2/3,
    ACup = rowMins(as.matrix(AC)) > 1.5,
    ACdn = rowMaxs(as.matrix(AC)) < 2/3,
    BCup = rowMins(as.matrix(BC)) > 1.5,
    BCdn = rowMaxs(as.matrix(BC)) < 2/3
)
rownames(cfc) <- rownames(cpm)

rm(A,B,C,AB,BC,AC)


table(deseq$AB_padj < .05 & abs(deseq$AB_log2FoldChange) > 1.5, cfc$ABup | cfc$ABdn, useNA = 'ifany')

table(deseq$AC_padj < .05 & abs(deseq$AC_log2FoldChange) > 1.5, cfc$ACup | cfc$ACdn, useNA = 'ifany')

table(deseq$BC_padj < .05 & abs(deseq$BC_log2FoldChange) > 1.5, cfc$BCup | cfc$BCdn, useNA = 'ifany')




plot(deseq$AB_log2FoldChange, -log10(deseq$AB_padj), col='grey', main = 'Timepoint A vs. B', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$AB_padj < .05 & abs(deseq$AB_log2FoldChange) > 1.5)
points(deseq$AB_log2FoldChange[ind], -log10(deseq$AB_padj)[ind], col=4, pch=16)
ind <- which(cfc$ABup)
points(deseq$AB_log2FoldChange[ind], -log10(deseq$AB_padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$ABdn)
points(deseq$AB_log2FoldChange[ind], -log10(deseq$AB_padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2', 'all FC >1.5', 'all FC <-1.5'))


plot(deseq$AC_log2FoldChange, -log10(deseq$AC_padj), col='grey', main = 'Timepoint A vs. C', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$AC_padj < .05 & abs(deseq$AC_log2FoldChange) > 1.5)
points(deseq$AC_log2FoldChange[ind], -log10(deseq$AC_padj)[ind], col=4, pch=16)
ind <- which(cfc$ACup)
points(deseq$AC_log2FoldChange[ind], -log10(deseq$AC_padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$ACdn)
points(deseq$AC_log2FoldChange[ind], -log10(deseq$AC_padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2', 'all FC >1.5', 'all FC <-1.5'))


plot(deseq$BC_log2FoldChange, -log10(deseq$BC_padj), col='grey', main = 'Timepoint B vs. C', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$BC_padj < .05 & abs(deseq$BC_log2FoldChange) > 1.5)
points(deseq$BC_log2FoldChange[ind], -log10(deseq$BC_padj)[ind], col=4, pch=16)
ind <- which(cfc$BCup)
points(deseq$BC_log2FoldChange[ind], -log10(deseq$BC_padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$BCdn)
points(deseq$BC_log2FoldChange[ind], -log10(deseq$BC_padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2', 'all FC >1.5', 'all FC <-1.5'))






# A vs B
up <- rownames(deseq)[which(deseq$AB_log2FoldChange > 1.5 & deseq$AB_padj < .05 & cfc$ABup)]
write.table(up, file='~/Desktop/DEresults/MVS_timepoint/A_vs_B/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)
dn <- rownames(deseq)[which(deseq$AB_log2FoldChange < 1.5 & deseq$AB_padj < .05 & cfc$ABdn)]
write.table(dn, file='~/Desktop/DEresults/MVS_timepoint/A_vs_B/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

# A vs C
up <- rownames(deseq)[which(deseq$AC_log2FoldChange > 1.5 & deseq$AC_padj < .05 & cfc$ACup)]
write.table(up, file='~/Desktop/DEresults/MVS_timepoint/A_vs_C/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)
dn <- rownames(deseq)[which(deseq$AC_log2FoldChange < 1.5 & deseq$AC_padj < .05 & cfc$ACdn)]
write.table(dn, file='~/Desktop/DEresults/MVS_timepoint/A_vs_C/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

# B vs C
up <- rownames(deseq)[which(deseq$BC_log2FoldChange > 1.5 & deseq$BC_padj < .05 & cfc$BCup)]
write.table(up, file='~/Desktop/DEresults/MVS_timepoint/B_vs_C/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)
dn <- rownames(deseq)[which(deseq$BC_log2FoldChange < 1.5 & deseq$BC_padj < .05 & cfc$BCdn)]
write.table(dn, file='~/Desktop/DEresults/MVS_timepoint/B_vs_C/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)


