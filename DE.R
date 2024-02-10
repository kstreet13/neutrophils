source('cluster.R')

require(DESeq2)

pheno$timepoint <- factor(pheno$timepoint)
pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
pheno$mmp8[pheno$mmp8 == "MMP8+"] <- "MMP8p"
pheno$mmp8[pheno$mmp8 == "MMP8-"] <- "MMP8n"
pheno$mmp8 <- factor(pheno$mmp8)


# MMP8+ vs. MMP8-
ind <- which(!is.na(pheno$mmp8))

ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                             colData = pheno[ind,],
                             design= ~ mmp8)
res <- DESeq(ds)

# resultsNames(res)
tab <- results(res, name=resultsNames(res)[2])
rownames(tab)[which(tab$padj < .01)]

boxplot(log1p(counts['CHDH',]) ~ pheno$mmp8)



####

# DE one-vs-all clusters

pheno$cluster <- factor(pheno$cluster)

de <- NULL
for(c in unique(pheno$cluster)){
    clusFac <- as.character(pheno$cluster)
    clusFac[clusFac != c] <- 'other'
    clusFac <- factor(clusFac, levels = c('other',c))
    pheno$clusFac <- clusFac

    ds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = pheno,
                                 design= ~ clusFac)
    res <- DESeq(ds)
    
    for(n in resultsNames(res)[-1]){
        tab <- results(res, name=n)
        tab$cluster <- c
        #de <- unique(c(de, rownames(tab)[which(tab$padj < .01)]))
        de <- rbind(de, tab[which(tab$padj < .01), ])
    }
}
pheno$clusFac <- NULL


boxplot(log1p(counts['ZYX',]) ~ pheno$cluster)


zscore <- cpm[de,]
zscore <- (zscore - rowMeans(zscore)) / rowSds(zscore)
png(filename = '~/Desktop/heat.png', width = 1500, height = 2000)
heatmap(as.matrix(zscore), ColSideColors = as.character(pheno$cluster))
dev.off()
# not much structure




