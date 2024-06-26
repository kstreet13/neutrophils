# MMP8- Septic Shock signature: 
# Compare SS7-12 to the average of EC1-EC8

# for comparison
# DE using DESeq
{
    source('setup.R')
    require(DESeq2)
    pheno$timepoint <- factor(pheno$timepoint)
    pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
    
    ind <- which(pheno$group == 'EC' | pheno$mmp8 == 'MMP8-')
    ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                                 colData = pheno[ind,],
                                 design= ~group)
    res <- DESeq(ds)
    
    deseq <- results(res, alpha = .05, name = 'group_SS_vs_EC')
} # deseq

source('setup.R')

counts <- counts[, which(pheno$group == 'EC' | pheno$mmp8 == 'MMP8-')]

cpm <- t(t(1e6*counts)/colSums(counts))

pheno <- pheno[which(pheno$group == 'EC' | pheno$mmp8 == 'MMP8-'), ]
# check order
all(pheno$sample == c("EC1","EC2","EC3","EC4","EC5","EC6","EC7","EC8","SS1","SS2","SS3","SS4","SS5","SS6"))
all(colnames(counts) == c("EC1","EC2","EC3","EC4","EC5","EC6","EC7","EC8","SS1","SS2","SS3","SS4","SS5","SS6"))

EC <- cpm[,1:8]
ECbase <- rowMeans(EC)
SS <- cpm[,9:14]

FC <- SS/ECbase; FC[is.nan(FC)] <- 1 # 1/0 = Inf, 0/0 = NaN


require(matrixStats)
clfc <- data.frame(
    FCup = rowMins(as.matrix(FC)) > 1.5,
    FCdn = rowMaxs(as.matrix(FC)) < 2/3
)
rownames(clfc) <- rownames(cpm)

rm(FC)


table(deseq$padj < .05 & abs(deseq$log2FoldChange) > 1.5, clfc$FCup | clfc$FCdn, useNA = 'ifany')


plot(deseq$log2FoldChange, -log10(deseq$padj), col='grey', main = 'SS (MMP8-) vs. EC', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$padj < .05 & abs(deseq$log2FoldChange) > 1.5)
points(deseq$log2FoldChange[ind], -log10(deseq$padj)[ind], col=4, pch=16)
ind <- which(clfc$FCup)
points(deseq$log2FoldChange[ind], -log10(deseq$padj)[ind], col=3, cex = .75, pch=16)
ind <- which(clfc$FCdn)
points(deseq$log2FoldChange[ind], -log10(deseq$padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2 deseq', 'all FC >1.5', 'all FC <-1.5'))



# example problematic gene: AC020909.3
boxplot(cpm['AC020909.3',] ~ pheno$group); abline(h = ECbase['AC020909.3'], lty=2)
boxplot(log1p(cpm['AC020909.3',]) ~ pheno$group)

boxplot(cpm['MMP8',] ~ pheno$group); abline(h = ECbase['MMP8'], lty=2)
boxplot(log1p(cpm['MMP8',]) ~ pheno$group)



up <- rownames(deseq)[which(deseq$log2FoldChange > 1.5 & deseq$padj < .05 & clfc$FCup)]
write.table(up, file='~/Desktop/DEresults/MMP8n_vs_EC/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

dn <- rownames(deseq)[which(deseq$log2FoldChange < 1.5 & deseq$padj < .05 & clfc$FCdn)]
write.table(dn, file='~/Desktop/DEresults/MMP8n_vs_EC/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

