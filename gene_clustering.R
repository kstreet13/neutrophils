
require(readxl)
counts <- read_excel('data/Trauma neutrophil counts.xlsx', col_names = FALSE)

cd <- t(counts[1:3,])
colnames(cd) <- cd[1,]
cd <- cd[-1,]
cd <- as.data.frame(cd)
cd$time <- as.numeric(gsub('h','',cd$Time))
cd$time[is.na(cd$time)] <- -4
rownames(cd) <- cd$Sample.name

colnames(counts) <- counts[1,]
counts <- counts[-(1:3),]
counts <- as.matrix(counts)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
mode(counts) <- 'numeric'





require(DESeq2)
cd$timepoint <- factor(cd$time)

ds <- DESeqDataSetFromMatrix(countData = round(counts),
                             colData = cd,
                             design= ~timepoint)

res <- DESeq(ds)
fc <- sapply(resultsNames(res)[2:5], function(n){
    results(res, name = n)$log2FoldChange
})
rownames(fc) <- names(res)
# 1 problematic gene
# fc <- fc[rowVars(fc) > 0, ]

# only cluster genes with at least 1 significant difference from baseline
keep <- sapply(resultsNames(res)[2:5], function(n){
    results(res, name = n)$pvalue < .05
})
keep <- which(rowSums(keep) > 0)

# add 0s before normalizing so that starting point is considered in clustering
fc0 <- cbind(0, fc[keep, ])
fc.norm <- (fc0 - rowMeans(fc0)) / rowSds(fc0)

###################
# gene clustering #
###################
require(mclust)
set.seed(1)
mc <- Mclust(fc.norm, G=2:6)
table(mc$classification)

# or just k-means
# set.seed(1)
# km <- kmeans(fc.norm, centers = 6)

times <- sort(unique(cd$time))

layout(matrix(1:6, ncol=2, byrow=TRUE))
for(k in 1:6){
    plot(range(times), range(fc.norm[mc$classification==k,]), col='white', axes = FALSE,
         main = paste0('Cluster ',k, ' (', sum(mc$classification==k), ' genes)'),
         xlab = 'Time (hrs)', ylab = 'Normalized log2 FC')
    box(); axis(2); axis(1, at = times[-1])
    for(ii in which(mc$classification==k)){
        lines(times, fc.norm[ii,], col=rgb(0,0,0,.1))
    }
    for(gene in c('MMP9', 'MMP8', 'CEACAM1', 'BMX'))
    if(gene %in% rownames(fc.norm)[mc$classification==k]){
        lines(times, fc.norm[gene,], col=2)
        legend('topright',bty='n',legend = gene, text.col = 2)
    }
}





layout(matrix(1:6, ncol=2, byrow=TRUE))
for(k in 1:6){
    plot(range(times), range(fc0[mc$classification==k,]), col='white', axes = FALSE,
         main = paste0('Cluster ',k, ' (', sum(mc$classification==k), ' genes)'),
         xlab = 'Time (hrs)', ylab = 'log2 FC')
    box(); axis(2); axis(1, at = times[-1])
    abline(h = 0, lty=2)
    for(ii in which(mc$classification==k)){
        lines(times, fc0[ii,], col=rgb(0,0,0,.1))
    }
    for(gene in c('MMP9', 'MMP8', 'CEACAM1', 'BMX'))
        if(gene %in% rownames(fc0)[mc$classification==k]){
            lines(times, fc0[gene,], col=2)
            legend('topright',bty='n',legend = gene, text.col = 2)
        }
}


layout(1)
plot(range(times), range(fc0[c('MMP9', 'MMP8', 'CEACAM1', 'BMX'),]), col='white', axes = FALSE,
     main = 'Selected Genes',
     xlab = 'Time (hrs)', ylab = 'log2 FC')
box(); axis(2); axis(1, at = times[-1])
abline(h = 0, lty=2)
abline(v = 0, lty=2)
lines(times, fc0['MMP9',], col=1)
lines(times, fc0['MMP8',], col=2)
lines(times, fc0['CEACAM1',], col=3)
lines(times, fc0['BMX',], col=4)
legend('topright', bty='n', legend =c('MMP9', 'MMP8', 'CEACAM1', 'BMX'), col = 1:4, lty=1)

