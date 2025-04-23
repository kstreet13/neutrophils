source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

genelist <- unique(c(
    read.table('data/DEresults/SC_timepoint/A_vs_B/down.csv')[,1],
    read.table('data/DEresults/SC_timepoint/A_vs_B/up.csv')[,1],
    read.table('data/DEresults/SC_timepoint/A_vs_C/down.csv')[,1],
    read.table('data/DEresults/SC_timepoint/A_vs_C/up.csv')[,1],
    read.table('data/DEresults/SC_timepoint/B_vs_C/down.csv')[,1],
    read.table('data/DEresults/SC_timepoint/B_vs_C/up.csv')[,1]
))

heatdat <- as.matrix(cpm[genelist, pheno$group == 'SC'])
# force ordering
heatdat <- heatdat[,  c(paste0('SC',1:8,'A'),paste0('SC',1:8,'B'),paste0('SC',1:8,'C'))]
subpheno <- pheno[match(colnames(heatdat), pheno$sample), ]

require(NMF)
anncolors = list(
    timepoint = c(A = "green", B = "red", C = "purple"))
aheatmap(log1p(heatdat), 
         annCol = subpheno[,"timepoint",drop=FALSE],
         Colv = NA, scale = "row",
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         filename = '~/Desktop/heat_SC.png')
# Expression values are normalized by dividing each count by the total for that sample and multiplying by 1 million (Counts Per Million, CPM), then adding a pseudocount of 1 and taking the natural log, to control outliers. For the heatmap, each gene's normalized expression is then standardized to a mean of 0 and standard deviation of 1 to show differences between timepoints.



# PCA
require(Matrix); require(matrixStats); require(sparseMatrixStats)

rv <- rowVars(cpm[, pheno$group == 'SC'])
ind <- order(rv, decreasing = TRUE)[1:5000]


pca <- BiocSingular::runPCA(t(log1p(cpm[ind, pheno$group == 'SC'])), rank = sum(pheno$group == 'SC'))

pctvar <- pca$sdev^2 / sum(pca$sdev^2)

plot(pca$x[,1:2], col=pheno$color[pheno$group=='SC'], asp=1,
     xlab = paste0('PC1 (',format(100*pctvar[1], digits = 3),'%)'),
     ylab = paste0('PC2 (',format(100*pctvar[2], digits = 3),'%)'))
# %var explained



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
    
    deseq <- cbind(res1,res2,res3)
    
    rm(res,res1,res2,res3,ind,ds)
} # deseq


# messy black lines showing average gene expression at each timepoint. Kinda looks like a bat
zscores <- (heatdat - rowMeans(heatdat)) / rowSds(heatdat)
means <- t(apply(zscores,1,function(cts){
    c(mean(cts[1:8]), mean(cts[9:16]), mean(cts[17:24]))
}))

plot(c(1,3), range(means), col='white')
for(i in 1:nrow(means)){
    lines(1:3, means[i,], col=rgb(0,0,0,.3))
}



# x: maximum (absolute) change from baseline for each gene
# y: expression at timepoint A
# 1) all genes, 2) DE genes

x <- cbind(deseq$AB_log2FoldChange, deseq$AC_log2FoldChange)
ind <- ifelse(abs(x[,2]) > abs(x[,1]), 2, 1)
ind[is.na(ind)] <- 1
x <- sapply(1:nrow(x),function(i){
    x[i, ind[i]]
})

# stats for paper:
absx <- abs(x)
length(which(absx > 1.5 & absx <= 2.5))
length(which(absx > 2.5 & absx <= 5))
length(which(absx > 5 & absx <= 10))
length(which(absx > 10))



means <- t(apply(log1p(cpm),1,function(cts){
    c(mean(cts[1:8]), mean(cts[9:16]), mean(cts[17:24]))
}))
y <- means[,1]

# all genes
plot(x,y, col=rgb(0,0,0,.5),
     xlab = 'Max. (abs) fold change from baseline',
     ylab = 'Avg. normalized expression at baseline',
     main = 'All genes')
abline(v=0,lty=2,col=2)

# de genes
DEind <- which(deseq$AB_padj < .05 | deseq$AC_padj < .05 | deseq$BC_padj < .05)
plot(x[DEind], y[DEind], col=rgb(0,0,0,.5),
     xlab = 'Max. (abs) fold change from baseline',
     ylab = 'Avg. normalized expression at baseline',
     main = 'DE genes')
abline(v=0,lty=2,col=2)

hist(y[x < 0])
hist(y[x > 0])

############
# PATTERNS #
############
# define recovery:
# DE in A->B, not DE in A->C and B->C shift has opposite sign
recovSCup <- which(deseq$AB_padj < .05 & deseq$AB_log2FoldChange > 0 & deseq$BC_padj < .05 & deseq$AC_padj > .05)
recovSCdn <- which(deseq$AB_padj < .05 & deseq$AB_log2FoldChange < 0 & deseq$BC_padj < .05 & deseq$AC_padj > .05)
recovSC <- c(recovSCup, recovSCdn)

# early change/no recovery:
# DE in A->B, not in B->C
earlySCup <- which(deseq$AB_padj < .05 & deseq$AB_log2FoldChange > 0 & deseq$BC_padj > .05 & deseq$AC_padj < .05)
earlySCdn <- which(deseq$AB_padj < .05 & deseq$AB_log2FoldChange < 0 & deseq$BC_padj > .05 & deseq$AC_padj < .05)
earlySC <- c(earlySCup, earlySCdn)

# late induction:
# not DE in A->B, DE in B->C
lateSCup <- which(deseq$AB_padj > .05 & deseq$BC_padj < .05 & deseq$BC_log2FoldChange > 0 & deseq$AC_padj < .05)
lateSCdn <- which(deseq$AB_padj > .05 & deseq$BC_padj < .05 & deseq$BC_log2FoldChange < 0 & deseq$AC_padj < .05)
lateSC <- c(lateSCup, lateSCdn)

# continuous change:
# DE in A->B, DE in B->C, same sign
contSCup <- which(deseq$AB_padj < .05 & deseq$BC_padj < .05 & (sign(deseq$AB_log2FoldChange)==sign(deseq$BC_log2FoldChange)) & deseq$AB_log2FoldChange > 0)
contSCdn <- which(deseq$AB_padj < .05 & deseq$BC_padj < .05 & (sign(deseq$AB_log2FoldChange)==sign(deseq$BC_log2FoldChange)) & deseq$AB_log2FoldChange < 0)
contSC <- c(contSCup, contSCdn)



# BARPLOT
cc <- brewer.pal(5,'Set1')[c(1,3,2,5)]
catcounts <- c(Early = length(earlySC),
               Recovery = length(recovSC),
               Late = length(lateSC),
               Progressive = length(contSC))

barplot(catcounts, xlab = 'Pattern', ylab = 'Count', main = 'SC DE Gene Patterns')


sum(deseq$AB_padj < .05, na.rm=TRUE)
sum(deseq$AB_padj < .05 & deseq$AC_padj > .05 & (sign(deseq$AB_log2FoldChange) == -sign(deseq$BC_log2FoldChange)), na.rm = TRUE)
barplot(c(Recovery = 1924, NoRecovery = 4225-1924))

# explainer figure
plot(c(0,4), c(-2,2), col='white')
abline(h=0, lty=2)
lines(c(0,1,0), col=cc[2], lwd=3)
lines(c(0,-1,0), col=cc[2], lwd=3)
lines((1:3)-.05, c(0,1,1), col=cc[1], lwd=3)
lines((1:3)-.05, c(0,-1,-1), col=cc[1], lwd=3)
lines((1:3), c(0,0,1), col=cc[3], lwd=3)
lines((1:3), c(0,0,-1), col=cc[3], lwd=3)
lines((1:3)-.1, c(0,1,2), col=cc[4], lwd=3)
lines((1:3)-.1, c(0,-1,-2), col=cc[4], lwd=3)



pattern <- rep(NA, nrow(heatdat))
pattern[rownames(heatdat) %in% rownames(deseq)[recov.ind]] <- 'Recovery'
pattern[rownames(heatdat) %in% rownames(deseq)[early.ind]] <- 'Early'
pattern[rownames(heatdat) %in% rownames(deseq)[late.ind]] <- 'Late'
pattern[rownames(heatdat) %in% rownames(deseq)[cont.ind]] <- 'Progressive'
genemeta <- data.frame(Pattern = factor(pattern))

anncolors = list(
    timepoint = c(A = "green", B = "red", C = "purple"),
    Pattern = c(Early = cc[1], Recovery = cc[2], Late = cc[3], Progressive = cc[4]))
aheatmap(log1p(heatdat), 
         annCol = subpheno[,"timepoint",drop=FALSE],
         annRow = genemeta,
         Colv = NA, scale = "row",
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         filename = '~/Desktop/heat_SC.png')




# the B->C change tends to reverse the A->B change
cc <- rep('grey50', nrow(deseq))
cc[recov.ind] <- 3
cc[early.ind] <- 2
cc[late.ind] <- 4

# add counts of different groups
plot(range(-log10(deseq$AB_padj)*sign(deseq$AB_log2FoldChange), na.rm = TRUE),
     range(-log10(deseq$BC_padj)*sign(deseq$BC_log2FoldChange), na.rm = TRUE),
     col = 'white', asp = 1,
     xlab = 'Signed -log10 p-value for A vs. B',
     ylab = 'Signed -log10 p-value for B vs. C')
abline(h=0,lty=2);abline(v=0,lty=2)
points(-log10(deseq$AB_padj)*sign(deseq$AB_log2FoldChange), -log10(deseq$BC_padj)*sign(deseq$BC_log2FoldChange), col = alpha(cc, alpha=.5), pch=16)
legend('topright', col = c(2,4,3), pch = 16, bty = 'n', cex = 1,
       legend = c(paste0('Early change (n=',length(early.ind),')'),
                  paste0('Late change (n=',length(late.ind),')'),
                  paste0('Recovery (n=',length(recov.ind),')')))


plot(deseq$AB_log2FoldChange, deseq$BC_log2FoldChange, col = alpha(cc, alpha=.5), asp = 1)
abline(h=0,col=2);abline(v=0,col=2)






#############
# # LISTS # #
#############
earlySCup.genes <- rownames(deseq)[earlySCup]
earlySCdn.genes <- rownames(deseq)[earlySCdn]
recovSCup.genes <- rownames(deseq)[recovSCup]
recovSCdn.genes <- rownames(deseq)[recovSCdn]
lateSCup.genes <- rownames(deseq)[lateSCup]
lateSCdn.genes <- rownames(deseq)[lateSCdn]
contSCup.genes <- rownames(deseq)[contSCup]
contSCdn.genes <- rownames(deseq)[contSCdn]


write.table(earlySCup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/earlySCup.txt')
write.table(recovSCup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/recovSCup.txt')
write.table(lateSCup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/lateSCup.txt')
write.table(contSCup.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/contSCup.txt')
write.table(earlySCdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/earlySCdn.txt')
write.table(recovSCdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/recovSCdn.txt')
write.table(lateSCdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/lateSCdn.txt')
write.table(contSCdn.genes, quote = FALSE, col.names = FALSE, row.names = FALSE,
            file='~/Desktop/gene_lists/SC/contSCdn.txt')




################################
# Figure ??: A LINE GRAPH WITH 4 DIFFERENT COLORS SHOWING THE BEHAVIOR OF ALL GENES WITHIN EACH DYNAMIC CATEGORY (AS DEFINED IN THE SLIDE FOLLOWING FIGURE 2).
# Dataset: all SC samples
################################
earlySCgenes <- unique(c(earlySCup.genes, earlySCdn.genes))
recovSCgenes <- unique(c(recovSCup.genes, recovSCdn.genes))
lateSCgenes <- unique(c(lateSCup.genes, lateSCdn.genes))
contSCgenes <- unique(c(contSCup.genes, contSCdn.genes))

# messy line plot 
genelist <- c(earlySCgenes, recovSCgenes, lateSCgenes, contSCgenes)
linedat <- as.matrix(cpm[genelist, pheno$group == 'SC'])
linedat <- linedat[,  c(paste0('SC',1:4,'A'),paste0('SC',1:4,'B'),paste0('SC',1:4,'C'))]
zscores <- (linedat - rowMeans(linedat)) / rowSds(linedat)
means <- cbind(
    rowMeans(zscores[,paste0('SC',1:4,'A')]),
    rowMeans(zscores[,paste0('SC',1:4,'B')]),
    rowMeans(zscores[,paste0('SC',1:4,'C')])
)
ind <- sample(nrow(linedat))

# separated, without colors
cc <- rep('black',4)
layout(matrix(1:4,nrow=2,byrow = TRUE))
# early
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Early SC', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in earlySCgenes){
    lines(1:3, means[g,], col=alpha(cc[1], alpha=.2))
}
# recov
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Recovery SC', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in recovSCgenes){
    lines(1:3, means[g,], col=alpha(cc[2], alpha=.2))
}
# late
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Late SC', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in lateSCgenes){
    lines(1:3, means[g,], col=alpha(cc[3], alpha=.2))
}
# cont
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'Continuous SC', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in contSCgenes){
    lines(1:3, means[g,], col=alpha(cc[4], alpha=.2))
}
layout(1)
par(mar=c(5,4,4,2)+.1)

