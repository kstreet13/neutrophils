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
