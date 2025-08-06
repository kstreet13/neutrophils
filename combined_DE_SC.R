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
    
    deseq <- cbind(res1,res2,res3)

    rm(res,res1,res2,res3,ind,ds)
} # deseq

source('setup.R')

counts <- counts[, pheno$group == 'SC']

cpm <- t(t(1e6*counts)/colSums(counts))

pheno <- pheno[pheno$group == 'SC', ]
# check order
all(pheno$sample == c("SC1A","SC1B","SC1C","SC2A","SC2B","SC2C","SC3A","SC3B","SC3C","SC4A","SC4B","SC4C","SC5A","SC5B","SC5C","SC6A","SC6B","SC6C","SC7A","SC7B","SC7C","SC8A","SC8B","SC8C"))
all(colnames(counts) == c("SC1A","SC1B","SC1C","SC2A","SC2B","SC2C","SC3A","SC3B","SC3C","SC4A","SC4B","SC4C","SC5A","SC5B","SC5C","SC6A","SC6B","SC6C","SC7A","SC7B","SC7C","SC8A","SC8B","SC8C"))

A <- cpm[,grep('A$', pheno$sample)]
B <- cpm[,grep('B$', pheno$sample)]
C <- cpm[,grep('C$', pheno$sample)]

AB <- B/A; AB[is.nan(AB)] <- 1 # 1/0 = Inf, 0/0 = NaN
AC <- C/A; AC[is.nan(AC)] <- 1
BC <- C/B; BC[is.nan(BC)] <- 1


require(matrixStats)
cfc <- data.frame(
    ABup = rowMins(as.matrix(AB)) > 1.5,
    ABdn = rowMaxs(as.matrix(AB)) < 2/3,
    ACup = rowMins(as.matrix(AC)) > 1.5,
    ACdn = rowMaxs(as.matrix(AC)) < 2/3,
    BCup = rowMins(as.matrix(BC)) > 1.5,
    BCdn = rowMaxs(as.matrix(BC)) < 2/3
)
rownames(cfc) <- rownames(cpm)

rm(A,B,C,AB,BC,AC)


table(deseq$AB_padj < .05 & abs(deseq$AB_log2FoldChange) > 1.5, cfc$ABup | cfc$ABdn, useNA = 'ifany')

table(deseq$AC_padj < .05 & abs(deseq$AC_log2FoldChange) > 1.5, cfc$ACup | cfc$ACdn, useNA = 'ifany')

table(deseq$BC_padj < .05 & abs(deseq$BC_log2FoldChange) > 1.5, cfc$BCup | cfc$BCdn, useNA = 'ifany')




plot(deseq$AB_log2FoldChange, -log10(deseq$AB_padj), col='grey', main = 'Timepoint A vs. B', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$AB_padj < .05 & abs(deseq$AB_log2FoldChange) > 1.5)
points(deseq$AB_log2FoldChange[ind], -log10(deseq$AB_padj)[ind], col=4, pch=16)
ind <- which(cfc$ABup)
points(deseq$AB_log2FoldChange[ind], -log10(deseq$AB_padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$ABdn)
points(deseq$AB_log2FoldChange[ind], -log10(deseq$AB_padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2 deseq', 'all FC >1.5', 'all FC <-1.5'))


plot(deseq$AC_log2FoldChange, -log10(deseq$AC_padj), col='grey', main = 'Timepoint A vs. C', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$AC_padj < .05 & abs(deseq$AC_log2FoldChange) > 1.5)
points(deseq$AC_log2FoldChange[ind], -log10(deseq$AC_padj)[ind], col=4, pch=16)
ind <- which(cfc$ACup)
points(deseq$AC_log2FoldChange[ind], -log10(deseq$AC_padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$ACdn)
points(deseq$AC_log2FoldChange[ind], -log10(deseq$AC_padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2 deseq', 'all FC >1.5', 'all FC <-1.5'))


plot(deseq$BC_log2FoldChange, -log10(deseq$BC_padj), col='grey', main = 'Timepoint B vs. C', xlab='Avg Log2 FC', ylab='-log10 Padj', pch=16)
ind <- which(deseq$BC_padj < .05 & abs(deseq$BC_log2FoldChange) > 1.5)
points(deseq$BC_log2FoldChange[ind], -log10(deseq$BC_padj)[ind], col=4, pch=16)
ind <- which(cfc$BCup)
points(deseq$BC_log2FoldChange[ind], -log10(deseq$BC_padj)[ind], col=3, cex = .75, pch=16)
ind <- which(cfc$BCdn)
points(deseq$BC_log2FoldChange[ind], -log10(deseq$BC_padj)[ind], col=2, cex = .75, pch=16)
abline(v = c(-1.5,1.5), col=4); abline(h = -log10(.05), col=4)
legend('topleft',bty='n',cex=.75,pch=16, col = c(4,3,2), legend = c('DESeq2 deseq', 'all FC >1.5', 'all FC <-1.5'))


# A vs B
up <- rownames(deseq)[which(deseq$AB_log2FoldChange > 1.5 & deseq$AB_padj < .05 & cfc$ABup)]
write.table(up, file='~/Desktop/DEresults/SC_timepoint/A_vs_B/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)
dn <- rownames(deseq)[which(deseq$AB_log2FoldChange < -1.5 & deseq$AB_padj < .05 & cfc$ABdn)]
write.table(dn, file='~/Desktop/DEresults/SC_timepoint/A_vs_B/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

# A vs C
up <- rownames(deseq)[which(deseq$AC_log2FoldChange > 1.5 & deseq$AC_padj < .05 & cfc$ACup)]
write.table(up, file='~/Desktop/DEresults/SC_timepoint/A_vs_C/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)
dn <- rownames(deseq)[which(deseq$AC_log2FoldChange < -1.5 & deseq$AC_padj < .05 & cfc$ACdn)]
write.table(dn, file='~/Desktop/DEresults/SC_timepoint/A_vs_C/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)

# B vs C
up <- rownames(deseq)[which(deseq$BC_log2FoldChange > 1.5 & deseq$BC_padj < .05 & cfc$BCup)]
write.table(up, file='~/Desktop/DEresults/SC_timepoint/B_vs_C/up.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)
dn <- rownames(deseq)[which(deseq$BC_log2FoldChange < -1.5 & deseq$BC_padj < .05 & cfc$BCdn)]
write.table(dn, file='~/Desktop/DEresults/SC_timepoint/B_vs_C/down.csv', row.names = FALSE, quote = FALSE, col.names = FALSE)




# The number of consistently, substantially and significantly changed genes with
# absolute log fold changes 1.5-2.5; 2.5-5, 5-10, >10 for both the surgical
# control signature and severe vasoplegic syndrome signature (to be added to the
# text). For comparisons A-B and A-C
# A vs B
up <- deseq[which(deseq$AB_log2FoldChange > 1.5 & deseq$AB_padj < .05 & cfc$ABup), ]
dn <- deseq[which(deseq$AB_log2FoldChange < 1.5 & deseq$AB_padj < .05 & cfc$ABdn), ]
SC_AB_up <- table(cut(abs(up$AB_log2FoldChange), breaks = c(1.5,2.5,5,10,Inf)))
SC_AB_dn <- table(cut(abs(dn$AB_log2FoldChange), breaks = c(1.5,2.5,5,10,Inf)))

# A vs C
up <- deseq[which(deseq$AC_log2FoldChange > 1.5 & deseq$AC_padj < .05 & cfc$ACup), ]
dn <- deseq[which(deseq$AC_log2FoldChange < 1.5 & deseq$AC_padj < .05 & cfc$ACdn), ]
SC_AC_up <- table(cut(abs(up$AC_log2FoldChange), breaks = c(1.5,2.5,5,10,Inf)))
SC_AC_dn <- table(cut(abs(dn$AC_log2FoldChange), breaks = c(1.5,2.5,5,10,Inf)))

cbind(SC_AB_up, SC_AB_dn, SC_AC_up, SC_AC_dn)


###############
# Max FC bins #
###############
# Max fold change 1.5-2.5,2.5-5,5-10 at AB, BC, or AC
maxFC <- apply(cbind(deseq$AB_log2FoldChange, deseq$BC_log2FoldChange, deseq$AC_log2FoldChange), 1, function(x){
    if(is.na(x[1])){
        return(NA)
    }else{
        return(x[which.max(abs(x))])
    }
})
whichMaxFC <- apply(cbind(deseq$AB_log2FoldChange, deseq$BC_log2FoldChange, deseq$AC_log2FoldChange), 1, function(x){
    if(is.na(x[1])){
        return(0)
    }else{
        return(c("AB","BC","AC")[which.max(abs(x))])
    }
})

# only use significant genes
genelist <- unique(c(
    read.table('~/Desktop/DEresults/SC_timepoint/A_vs_B/down.csv')[,1],
    read.table('~/Desktop/DEresults/SC_timepoint/A_vs_B/up.csv')[,1],
    read.table('~/Desktop/DEresults/SC_timepoint/A_vs_C/down.csv')[,1],
    read.table('~/Desktop/DEresults/SC_timepoint/A_vs_C/up.csv')[,1],
    read.table('~/Desktop/DEresults/SC_timepoint/B_vs_C/down.csv')[,1],
    read.table('~/Desktop/DEresults/SC_timepoint/B_vs_C/up.csv')[,1]
))
idx <- which(rownames(deseq) %in% genelist)

tab <- table(cut(maxFC[idx], breaks=c(-Inf,-10,-5,-2.5,-1.5,1.5,2.5,5,10,Inf)),
             whichMaxFC[idx])

write.csv(tab, file='~/Desktop/DEresults/SC_timepoint/maxAbsFCbins.csv')


############
# PATTERNS #
############
ABdn <- read.table('~/Desktop/DEresults/SC_timepoint/A_vs_B/down.csv')[,1]
ABup <- read.table('~/Desktop/DEresults/SC_timepoint/A_vs_B/up.csv')[,1]
ACdn <- read.table('~/Desktop/DEresults/SC_timepoint/A_vs_C/down.csv')[,1]
ACup <- read.table('~/Desktop/DEresults/SC_timepoint/A_vs_C/up.csv')[,1]
BCdn <- read.table('~/Desktop/DEresults/SC_timepoint/B_vs_C/down.csv')[,1]
BCup <- read.table('~/Desktop/DEresults/SC_timepoint/B_vs_C/up.csv')[,1]

# "DE" means by all 3 criteria
# recovery:
# DE in A->B, DE in B->C (opposite direction), not DE in A->C
recovSCup <- which(rownames(deseq) %in% ABup & rownames(deseq) %in% BCdn & !(rownames(deseq) %in% c(ACup,ACdn)))
recovSCdn <- which(rownames(deseq) %in% ABdn & rownames(deseq) %in% BCup & !(rownames(deseq) %in% c(ACup,ACdn)))
recovSC <- c(recovSCup, recovSCdn)

# early change/no recovery:
# DE in A->B, A->C, not in B->C
earlySCup <- which(rownames(deseq) %in% ABup & !(rownames(deseq) %in% c(BCup,BCdn)) & rownames(deseq) %in% ACup)
earlySCdn <- which(rownames(deseq) %in% ABdn & !(rownames(deseq) %in% c(BCup,BCdn)) & rownames(deseq) %in% ACdn)
earlySC <- c(earlySCup, earlySCdn)

# late induction:
# not DE in A->B, DE in B->C, A->C
lateSCup <- which(!(rownames(deseq) %in% c(ABup,ABdn)) & rownames(deseq) %in% BCup & rownames(deseq) %in% ACup)
lateSCdn <- which(!(rownames(deseq) %in% c(ABup,ABdn)) & rownames(deseq) %in% BCdn & rownames(deseq) %in% ACdn)
lateSC <- c(lateSCup, lateSCdn)

# continuous change:
# DE in A->B, DE in B->C, same sign
contSCup <- which(rownames(deseq) %in% ABup & rownames(deseq) %in% BCup)
contSCdn <- which(rownames(deseq) %in% ABdn & rownames(deseq) %in% BCdn)
contSC <- c(contSCup, contSCdn)




# get gene names, not just indices
earlySCup.genes <- rownames(deseq)[earlySCup]
earlySCdn.genes <- rownames(deseq)[earlySCdn]
recovSCup.genes <- rownames(deseq)[recovSCup]
recovSCdn.genes <- rownames(deseq)[recovSCdn]
lateSCup.genes <- rownames(deseq)[lateSCup]
lateSCdn.genes <- rownames(deseq)[lateSCdn]
contSCup.genes <- rownames(deseq)[contSCup]
contSCdn.genes <- rownames(deseq)[contSCdn]



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

