source('cluster.R')

require(DESeq2)

pheno$timepoint <- factor(pheno$timepoint)
pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))

ind <- which(pheno$group %in% c('SC','MVS','SVS'))

ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                             colData = pheno[ind,],
                             design= ~ group + timepoint)
res <- DESeq(ds)


tab <- results(res, name=resultsNames(res)[3])
rownames(tab)[which(tab$padj < .01)]





####

# find genes DE between any pair of labels

pheno$color <- factor(pheno$color)

de <- NULL
for(c in unique(pheno$color)){
    pheno$color <- relevel(pheno$color, ref = c)
    
    ds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = pheno,
                                 design= ~ color)
    res <- DESeq(ds)
    
    for(n in resultsNames(res)[-1]){
        tab <- results(res, name=n)
        de <- unique(c(de, rownames(tab)[which(tab$padj < .01)]))
    }
}



zscore <- cpm[de,]
zscore <- (zscore - rowMeans(zscore)) / rowSds(zscore)
png(filename = '~/Desktop/heat.png', width = 1500, height = 2000)
heatmap(as.matrix(zscore), ColSideColors = as.character(pheno$color))
dev.off()
# not much structure




