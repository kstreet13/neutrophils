source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

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
    
    deseqMVS <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseqMVS
{
    source('setup.R')
    require(DESeq2)
    pheno$timepoint <- factor(pheno$timepoint)
    pheno$group <- factor(pheno$group, levels = c('SC','MVS','SVS','EC','SS'))
    
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
    
    deseqSVS <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseqSVS


############
# PATTERNS #
############
# define recovery:
# DE in A->B, not DE in A->C and B->C shift has opposite sign
recovMVSup <- which(deseqMVS$AB_padj < .05 & deseqMVS$AB_log2FoldChange > 0 & deseqMVS$BC_padj < .05 & deseqMVS$AC_padj > .05)
recovMVSdn <- which(deseqMVS$AB_padj < .05 & deseqMVS$AB_log2FoldChange < 0 & deseqMVS$BC_padj < .05 & deseqMVS$AC_padj > .05)
recovMVS <- c(recovMVSup, recovMVSdn)

recovSVSup <- which(deseqSVS$AB_padj < .05 & deseqSVS$AB_log2FoldChange > 0 & deseqSVS$BC_padj < .05 & deseqSVS$AC_padj > .05)
recovSVSdn <- which(deseqSVS$AB_padj < .05 & deseqSVS$AB_log2FoldChange < 0 & deseqSVS$BC_padj < .05 & deseqSVS$AC_padj > .05)
recovSVS <- c(recovSVSup, recovSVSdn)


# early change/no recovery:
# DE in A->B, not in B->C
earlyMVSup <- which(deseqMVS$AB_padj < .05 & deseqMVS$AB_log2FoldChange > 0 & deseqMVS$BC_padj > .05 & deseqMVS$AC_padj < .05)
earlyMVSdn <- which(deseqMVS$AB_padj < .05 & deseqMVS$AB_log2FoldChange < 0 & deseqMVS$BC_padj > .05 & deseqMVS$AC_padj < .05)
earlyMVS <- c(earlyMVSup, earlyMVSdn)

earlySVSup <- which(deseqSVS$AB_padj < .05 & deseqSVS$AB_log2FoldChange > 0 & deseqSVS$BC_padj > .05 & deseqSVS$AC_padj < .05)
earlySVSdn <- which(deseqSVS$AB_padj < .05 & deseqSVS$AB_log2FoldChange < 0 & deseqSVS$BC_padj > .05 & deseqSVS$AC_padj < .05)
earlySVS <- c(earlySVSup, earlySVSdn)

# late induction:
# not DE in A->B, DE in B->C
lateMVSup <- which(deseqMVS$AB_padj > .05 & deseqMVS$BC_padj < .05 & deseqMVS$BC_log2FoldChange > 0 & deseqMVS$AC_padj < .05)
lateMVSdn <- which(deseqMVS$AB_padj > .05 & deseqMVS$BC_padj < .05 & deseqMVS$BC_log2FoldChange < 0 & deseqMVS$AC_padj < .05)
lateMVS <- c(lateMVSup, lateMVSdn)

lateSVSup <- which(deseqSVS$AB_padj > .05 & deseqSVS$BC_padj < .05 & deseqSVS$BC_log2FoldChange > 0 & deseqSVS$AC_padj < .05)
lateSVSdn <- which(deseqSVS$AB_padj > .05 & deseqSVS$BC_padj < .05 & deseqSVS$BC_log2FoldChange < 0 & deseqSVS$AC_padj < .05)
lateSVS <- c(lateSVSup, lateSVSdn)

# continuous change:
# DE in A->B, DE in B->C, same sign
contMVSup <- which(deseqMVS$AB_padj < .05 & deseqMVS$BC_padj < .05 & (sign(deseqMVS$AB_log2FoldChange)==sign(deseqMVS$BC_log2FoldChange)) & deseqMVS$AB_log2FoldChange > 0)
contMVSdn <- which(deseqMVS$AB_padj < .05 & deseqMVS$BC_padj < .05 & (sign(deseqMVS$AB_log2FoldChange)==sign(deseqMVS$BC_log2FoldChange)) & deseqMVS$AB_log2FoldChange < 0)
contMVS <- c(contMVSup, contMVSdn)

contSVSup <- which(deseqSVS$AB_padj < .05 & deseqSVS$BC_padj < .05 & (sign(deseqSVS$AB_log2FoldChange)==sign(deseqSVS$BC_log2FoldChange)) & deseqSVS$AB_log2FoldChange > 0)
contSVSdn <- which(deseqSVS$AB_padj < .05 & deseqSVS$BC_padj < .05 & (sign(deseqSVS$AB_log2FoldChange)==sign(deseqSVS$BC_log2FoldChange)) & deseqSVS$AB_log2FoldChange < 0)
contSVS <- c(contSVSup, contSVSdn)










#############
# # LISTS # #
#############
earlyMVSup.genes <- rownames(deseqMVS)[earlyMVSup]
earlyMVSdn.genes <- rownames(deseqMVS)[earlyMVSdn]
recovMVSup.genes <- rownames(deseqMVS)[recovMVSup]
recovMVSdn.genes <- rownames(deseqMVS)[recovMVSdn]
lateMVSup.genes <- rownames(deseqMVS)[lateMVSup]
lateMVSdn.genes <- rownames(deseqMVS)[lateMVSdn]
contMVSup.genes <- rownames(deseqMVS)[contMVSup]
contMVSdn.genes <- rownames(deseqMVS)[contMVSdn]

earlySVSup.genes <- rownames(deseqSVS)[earlySVSup]
earlySVSdn.genes <- rownames(deseqSVS)[earlySVSdn]
recovSVSup.genes <- rownames(deseqSVS)[recovSVSup]
recovSVSdn.genes <- rownames(deseqSVS)[recovSVSdn]
lateSVSup.genes <- rownames(deseqSVS)[lateSVSup]
lateSVSdn.genes <- rownames(deseqSVS)[lateSVSdn]
contSVSup.genes <- rownames(deseqSVS)[contSVSup]
contSVSdn.genes <- rownames(deseqSVS)[contSVSdn]

earlyMVSgenes <- c(earlyMVSup.genes, earlyMVSdn.genes)
recovMVSgenes <- c(recovMVSup.genes, recovMVSdn.genes)
lateMVSgenes <- c(lateMVSup.genes, lateMVSdn.genes)
contMVSgenes <- c(contMVSup.genes, contMVSdn.genes)
earlySVSgenes <- c(earlySVSup.genes, earlySVSdn.genes)
recovSVSgenes <- c(recovSVSup.genes, recovSVSdn.genes)
lateSVSgenes <- c(lateSVSup.genes, lateSVSdn.genes)
contSVSgenes <- c(contSVSup.genes, contSVSdn.genes)



write.table(earlyMVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/earlyMVSup.txt')
write.table(recovMVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/recovMVSup.txt')
write.table(lateMVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/lateMVSup.txt')
write.table(contMVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/contMVSup.txt')
write.table(earlyMVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/earlyMVSdn.txt')
write.table(recovMVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/recovMVSdn.txt')
write.table(lateMVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/lateMVSdn.txt')
write.table(contMVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/MVS/contMVSdn.txt')


write.table(earlySVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/earlySVSup.txt')
write.table(recovSVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/recovSVSup.txt')
write.table(lateSVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/lateSVSup.txt')
write.table(contSVSup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/contSVSup.txt')
write.table(earlySVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/earlySVSdn.txt')
write.table(recovSVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/recovSVSdn.txt')
write.table(lateSVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/lateSVSdn.txt')
write.table(contSVSdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SVS/contSVSdn.txt')


############
# BARPLOTS #
############
catcounts <- matrix(c(
    length(earlyMVS),length(earlySVS),
    length(recovMVS),length(recovSVS),
    length(lateMVS), length(lateSVS),
    length(contMVS), length(contSVS)
), nrow=2)
rownames(catcounts) <- c('MVS','SVS')
colnames(catcounts) <- c('Early','Recovery','Late','Progressive')

# barplot number of genes in each category (SVS only)
barplot(catcounts[2,], ylab = 'Count', main = 'SVS DE Gene Patterns')




################################
# Figure 2E: A LINE GRAPH WITH 4 DIFFERENT COLORS SHOWING THE BEHAVIOR OF ALL GENES WITHIN EACH DYNAMIC CATEGORY (AS DEFINED IN THE SLIDE FOLLOWING FIGURE 2).
# Dataset: all SVS samples
################################

# messy line plot 
genelist <- c(earlySVSgenes, recovSVSgenes, lateSVSgenes, contSVSgenes)
linedat <- as.matrix(cpm[genelist, pheno$group == 'SVS'])
linedat <- linedat[,  c(paste0('SVS',1:4,'A'),paste0('SVS',1:4,'B'),paste0('SVS',1:4,'C'))]
zscores <- (linedat - rowMeans(linedat)) / rowSds(linedat)
means <- cbind(
    rowMeans(zscores[,paste0('SVS',1:4,'A')]),
    rowMeans(zscores[,paste0('SVS',1:4,'B')]),
    rowMeans(zscores[,paste0('SVS',1:4,'C')])
)
ind <- sample(nrow(linedat))

# separated, without colors
cc <- rep('black',4)
layout(matrix(1:4,nrow=2,byrow = TRUE))
# early
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Early SVS', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in earlySVSgenes){
    lines(1:3, means[g,], col=alpha(cc[1], alpha=.2))
}
# recov
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Recovery SVS', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in recovSVSgenes){
    lines(1:3, means[g,], col=alpha(cc[2], alpha=.2))
}
# late
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Late SVS', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in lateSVSgenes){
    lines(1:3, means[g,], col=alpha(cc[3], alpha=.2))
}
# cont
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Continuous SVS', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in contSVSgenes){
    lines(1:3, means[g,], col=alpha(cc[4], alpha=.2))
}
layout(1)
par(mar=c(5,4,4,2)+.1)

