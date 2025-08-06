# Septic Shock signature: 
# Compare SS1-12 to the average of EC1-EC8 for CFC

# for comparison
# DE using DESeq
{
    source('setup.R')
    require(DESeq2)
    pheno$timepoint <- factor(pheno$timepoint)
    pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
    
    ind <- which(pheno$group %in% c('SS','EC'))
    ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                                 colData = pheno[ind,],
                                 design= ~group)
    res <- DESeq(ds)
    
    deseq <- results(res, alpha = .05, name = 'group_SS_vs_EC')
} # deseq

source('setup.R')

counts <- counts[, pheno$group %in% c('SS','EC')]

cpm <- t(t(1e6*counts)/colSums(counts))

pheno <- pheno[pheno$group %in% c('SS','EC'), ]
# check order
all(pheno$sample == c("EC1","EC2","EC3","EC5","EC6","EC7","EC8","SS1","SS10","SS11","SS12","SS2","SS3","SS4","SS5","SS6","SS7","SS8","SS9"))
all(colnames(counts) == c("EC1","EC2","EC3","EC5","EC6","EC7","EC8","SS1","SS10","SS11","SS12","SS2","SS3","SS4","SS5","SS6","SS7","SS8","SS9"))

EC <- cpm[,grep('EC', colnames(cpm))]
ECbase <- rowMeans(EC)
SS <- cpm[,grep('SS', colnames(cpm))]

FC <- SS/ECbase; FC[is.nan(FC)] <- 1 # 1/0 = Inf, 0/0 = NaN


require(matrixStats)
cfc <- data.frame(
    FCup = rowMins(as.matrix(FC)) > 1.5,
    FCdn = rowMaxs(as.matrix(FC)) < 2/3
)
rownames(cfc) <- rownames(cpm)

rm(FC)


table(deseq$padj < .05 & abs(deseq$log2FoldChange) > 1.5, cfc$FCup | cfc$FCdn, useNA = 'ifany')


plot(deseq$log2FoldChange, -log10(deseq$padj), col='grey', main = 'SS vs. EC', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$padj < .05 & abs(deseq$log2FoldChange) > 1.5)
points(deseq$log2FoldChange[ind], -log10(deseq$padj)[ind], col=4, pch=16)
ind <- which(cfc$FCup)
points(deseq$log2FoldChange[ind], -log10(deseq$padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$FCdn)
points(deseq$log2FoldChange[ind], -log10(deseq$padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2 deseq', 'all FC >1.5', 'all FC <-1.5'))



# example problematic gene: A2M has positive log2FC (by DESeq2), but is "consistently down"
# it's because DESeq2 handles outliers, but simple average doesn't
boxplot(cpm['A2M',] ~ pheno$group); abline(h = ECbase['A2M'], lty=2)
boxplot(log1p(cpm['A2M',]) ~ pheno$group); abline(h = log1p(ECbase['A2M']), lty=2)



up <- rownames(deseq)[which(deseq$log2FoldChange > 1.5 & deseq$padj < .05 & cfc$FCup)]
write.table(up, file='data/DEresults/SS_vs_EC/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

dn <- rownames(deseq)[which(deseq$log2FoldChange < -1.5 & deseq$padj < .05 & cfc$FCdn)]
write.table(dn, file='data/DEresults/SS_vs_EC/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)



