source('cluster.R')

require(DESeq2)

pheno$timepoint <- factor(pheno$timepoint)
pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
pheno$mmp8[pheno$mmp8 == "MMP8+"] <- "MMP8p"
pheno$mmp8[pheno$mmp8 == "MMP8-"] <- "MMP8n"
pheno$mmp8 <- factor(pheno$mmp8)



# Full model??
############
p <- pheno
p$timepoint[p$timepoint==0] <- 'A'
p$timepoint <- factor(p$timepoint)

ds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = p,
                             design= ~ group + timepoint)
res <- DESeq(ds)





# Surgical controls (by timepoint)
###################
# all three pairwise comparisons between timepoints within the SC group
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

hits <- cbind(res1,res2,res3)
hits <- hits[which((hits$AB_padj < .05 & abs(hits$AB_log2FoldChange) > 1.5) | 
                   (hits$AC_padj < .05 & abs(hits$AC_log2FoldChange) > 1.5) | 
                   (hits$BC_padj < .05 & abs(hits$BC_log2FoldChange) > 1.5)), ]

heatmat <- log1p(as.matrix(cpm[rownames(hits), pheno$group=='SC']))
heatmat <- (heatmat-rowMeans(heatmat)) / rowSds(heatmat)

require(pheatmap)
rownames(pheno) <- pheno$sample
pheno$timepoint <- factor(pheno$timepoint, levels = c('A','B','C'))

annGene <- data.frame(
    AvsB = rep(NA, nrow(heatmat)),
    AvsC = rep(NA, nrow(heatmat)),
    BvsC = rep(NA, nrow(heatmat))
)
annGene$AvsB[hits$AB_padj < .05 & hits$AB_log2FoldChange < 0] <- 'A'
annGene$AvsB[hits$AB_padj < .05 & hits$AB_log2FoldChange > 0] <- 'B'
annGene$AvsC[hits$AC_padj < .05 & hits$AC_log2FoldChange < 0] <- 'A'
annGene$AvsC[hits$AC_padj < .05 & hits$AC_log2FoldChange > 0] <- 'C'
annGene$BvsC[hits$BC_padj < .05 & hits$BC_log2FoldChange < 0] <- 'B'
annGene$BvsC[hits$BC_padj < .05 & hits$BC_log2FoldChange > 0] <- 'C'
rownames(annGene) <- rownames(hits)

ann_colors = list(
    timepoint = c(A = "green", B = "red", C = "purple"),
    AvsB = c(A = "green", B = "red"),
    AvsC = c(A = "green", C = "purple"),
    BvsC = c(B = "red", C = "purple")
)
pheatmap(heatmat, annotation_col = pheno[pheno$group=='SC','timepoint',drop=FALSE], annotation_colors = ann_colors, labels_row = rep('',nrow(heatmat)), annotation_row = annGene, cluster_cols = FALSE, treeheight_row = 0)





# MVS/SVS (by timepoint)
###################
# all three pairwise comparisons between timepoints within the MVS/SVS groups
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

hits <- cbind(res1,res2,res3)
hits <- hits[which((hits$AB_padj < .05 & abs(hits$AB_log2FoldChange) > 1.5) | 
                       (hits$AC_padj < .05 & abs(hits$AC_log2FoldChange) > 1.5) | 
                       (hits$BC_padj < .05 & abs(hits$BC_log2FoldChange) > 1.5)), ]

heatmat <- log1p(as.matrix(cpm[rownames(hits), pheno$group=='SC']))
heatmat <- (heatmat-rowMeans(heatmat)) / rowSds(heatmat)

require(pheatmap)
rownames(pheno) <- pheno$sample
pheno$timepoint <- factor(pheno$timepoint, levels = c('A','B','C'))

annGene <- data.frame(
    AvsB = rep(NA, nrow(heatmat)),
    AvsC = rep(NA, nrow(heatmat)),
    BvsC = rep(NA, nrow(heatmat))
)
annGene$AvsB[hits$AB_padj < .05 & hits$AB_log2FoldChange < 0] <- 'A'
annGene$AvsB[hits$AB_padj < .05 & hits$AB_log2FoldChange > 0] <- 'B'
annGene$AvsC[hits$AC_padj < .05 & hits$AC_log2FoldChange < 0] <- 'A'
annGene$AvsC[hits$AC_padj < .05 & hits$AC_log2FoldChange > 0] <- 'C'
annGene$BvsC[hits$BC_padj < .05 & hits$BC_log2FoldChange < 0] <- 'B'
annGene$BvsC[hits$BC_padj < .05 & hits$BC_log2FoldChange > 0] <- 'C'
rownames(annGene) <- rownames(hits)

ann_colors = list(
    timepoint = c(A = "green", B = "red", C = "purple"),
    AvsB = c(A = "green", B = "red"),
    AvsC = c(A = "green", C = "purple"),
    BvsC = c(B = "red", C = "purple")
)
pheatmap(heatmat, annotation_col = pheno[pheno$group=='SC','timepoint',drop=FALSE], annotation_colors = ann_colors, labels_row = rep('',nrow(heatmat)), annotation_row = annGene, cluster_cols = FALSE, treeheight_row = 0)




# SS vs EC 
###########
ds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = pheno,
                             design= ~group)
res <- DESeq(ds)
res1 <- results(res, alpha = .05, contrast = c('group', 'EC', 'SS'))


hits <- res1
hits <- hits[which(hits$padj < .05 & abs(hits$log2FoldChange) > 1.5), ]




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




