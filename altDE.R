# alternate surgical control (SC) DE list
# at least 1.5x (not log) difference in every pair of samples between two timepoints

# for comparison
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
    
    hits <- cbind(res1,res2,res3)
    hits <- hits[which((hits$AB_padj < .05 & abs(hits$AB_log2FoldChange) > 1.5) | 
                           (hits$AC_padj < .05 & abs(hits$AC_log2FoldChange) > 1.5) | 
                           (hits$BC_padj < .05 & abs(hits$BC_log2FoldChange) > 1.5)), ]
    
    rm(res,res1,res2,res3,ind,ds)
} # hits

source('setup.R')

counts <- counts[, pheno$group == 'SC']
pheno <- pheno[pheno$group == 'SC', ]
# check order
all(pheno$sample == c("SC1A","SC2A","SC3A","SC4A","SC5A","SC6A","SC7A","SC8A","SC1B","SC2B","SC3B","SC4B","SC5B","SC6B","SC7B","SC8B","SC1C","SC2C","SC3C","SC4C","SC5C","SC6C","SC7C","SC8C"))
all(colnames(counts) == c("SC1A","SC2A","SC3A","SC4A","SC5A","SC6A","SC7A","SC8A","SC1B","SC2B","SC3B","SC4B","SC5B","SC6B","SC7B","SC8B","SC1C","SC2C","SC3C","SC4C","SC5C","SC6C","SC7C","SC8C"))

A <- counts[,1:8]
B <- counts[,9:16]
C <- counts[,17:24]

AB <- A/B; AB[is.nan(AB)] <- 1 # 1/0 = Inf, 0/0 = NaN
AC <- A/C; AC[is.nan(AC)] <- 1
BC <- B/C; BC[is.nan(BC)] <- 1


require(matrixStats)
res <- data.frame(
    ABup = rowMins(as.matrix(AB)) > 1.5,
    ABdn = rowMaxs(as.matrix(AB)) < 2/3,
    ACup = rowMins(as.matrix(AC)) > 1.5,
    ACdn = rowMaxs(as.matrix(AC)) < 2/3,
    BCup = rowMins(as.matrix(BC)) > 1.5,
    BCdn = rowMaxs(as.matrix(BC)) < 2/3
)



ABhits.DESeq <- rownames(hits)[hits$AB_padj < .05]
ABhits.FC <- rownames(counts)[which(res$ABup | res$ABdn)]
x <- unique(c(ABhits.DESeq, ABhits.FC))
table(x %in% ABhits.DESeq, x %in% ABhits.FC)


plot(hits$AB_log2FoldChange, -log10(hits$AB_padj))
abline(v = c(-1.5,1.5)); abline(h = -log10(.05))
ind <- which(rownames(hits) %in% ABhits.FC)
points(hits$AB_log2FoldChange[ind], -log10(hits$AB_padj)[ind], col=2)






