source('cluster.R')

require(DESeq2)

pheno$timepoint <- factor(pheno$timepoint)
pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
pheno$mmp8[pheno$mmp8 == "MMP8+"] <- "MMP8p"
pheno$mmp8[pheno$mmp8 == "MMP8-"] <- "MMP8n"
pheno$mmp8 <- factor(pheno$mmp8)



# Surgical controls
###################
# all three pairwise comparisons between timepoints within the SC group
ind <- which(pheno$group=='SC')
ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                             colData = pheno[ind,],
                             design= ~timepoint)
res <- DESeq(ds)
res1 <- results(res, alpha = .05, name = 'timepoint_B_vs_A')
names(res1) <- paste0('AB_', names(res1))
res2 <- results(res, alpha = .05, name = 'timepoint_C_vs_A')
names(res2) <- paste0('AC_', names(res2))

pheno$timepoint <- factor(pheno$timepoint, levels = c('B','C','A'))
ds <- DESeqDataSetFromMatrix(countData = counts[,ind],
                             colData = pheno[ind,],
                             design= ~timepoint)
res <- DESeq(ds)
res3 <- results(res, alpha = .05, name = 'timepoint_C_vs_B')
names(res3) <- paste0('BC_', names(res3))


hits <- cbind(res1,res2,res3)
hits <- hits[which((hits$AB_padj < .05 & abs(hits$AB_log2FoldChange) > 1.5) | 
                   (hits$AC_padj < .05 & abs(hits$AC_log2FoldChange) > 1.5) | 
                   (hits$BC_padj < .05 & abs(hits$BC_log2FoldChange) > 1.5)), ]

heatmat <- as.matrix(cpm[rownames(hits), pheno$group=='SC'])
heatmat <- (heatmat-rowMeans(heatmat)) / rowSds(heatmat)

require(pheatmap)
rownames(pheno) <- pheno$sample
pheno$timepoint <- factor(pheno$timepoint, levels = c('A','B','C'))

annGene <- data.frame(
    AB = rep(NA, nrow(heatmat)),
    AC = rep(NA, nrow(heatmat)),
    BC = rep(NA, nrow(heatmat))
)
annGene$AB[hits$AB_padj < .05 & hits$AB_log2FoldChange < 0] <- 'green'
annGene$AB[hits$AB_padj < .05 & hits$AB_log2FoldChange > 0] <- 'red'
annGene$AC[hits$AC_padj < .05 & hits$AC_log2FoldChange < 0] <- 'green'
annGene$AC[hits$AC_padj < .05 & hits$AC_log2FoldChange > 0] <- 'purple'
annGene$BC[hits$BC_padj < .05 & hits$BC_log2FoldChange < 0] <- 'red'
annGene$BC[hits$BC_padj < .05 & hits$BC_log2FoldChange > 0] <- 'purple'
rownames(annGene) <- rownames(hits)

ann_colors = list(
    timepoint = c(A = "green", B = "red", C = "purple"),
    AB = c(green = "green", red = "red"),
    AC = c(green = "green", purple = "purple"),
    BC = c(red = "red", purple = "purple")
)
pheatmap(heatmat, annotation_col = pheno[pheno$group=='SC','timepoint',drop=FALSE], annotation_colors = ann_colors, labels_row = rep('',nrow(heatmat)), annotation_row = annGene)







# MMP8+ vs. MMP8-
#################
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
########################

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
        
        # fix issue with rownames
        
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




