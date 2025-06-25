source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

#######
# PCA #
#######
require(matrixStats); require(sparseMatrixStats)
rm <- rowMeans(cpm)
rv <- rowVars(cpm)
filt <- which(rm > 1 & rv > 1)
#filt <- which(rowVars(cpm) >= sort(rowVars(cpm), decreasing = TRUE)[10000])

require(BiocSingular)
pca <- BiocSingular::runPCA(t(log1p(cpm[filt, ])), rank=64)


### OLD VERSIONS
# PC1 & PC2
# plot(pca$x, asp=1, col=pheno$color, cex=3)
# 3D
# require(plotly)
# plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3, color = pheno$color, colors = sort(unique(pheno$color)))
###


# open circles for controls, filled circles for mild vasoplegic syndrome, filled
# triangles or squares (or whatever other shape) for severe vasoplegic
# syndrome. Instead of having dark and light red/purple indicate mild vs.
# severe, the shapes will indicate this difference and there will be only one
# shade of green, red, and purple.

# set up new plotting variables
pheno$shape <- 1
pheno$shape[pheno$group=='MVS'] <- 16
pheno$shape[pheno$group=='SVS'] <- 15
pheno$color2 <- pheno$color
pheno$color2[pheno$color=='firebrick'] <- 'red'
pheno$color2[pheno$color=='forestgreen'] <- 'green'
pheno$color2[pheno$color=='slateblue4'] <- 'purple'
pheno$color2[pheno$mmp8=='MMP8+'] <- 'brown4'
#pheno$color2[pheno$color=='purple'] <- 'blue'
#pheno$color2[pheno$group=='EC'] <- 'purple'
pheno$shape_plotly <- as.character(pheno$shape)
pheno$shape_plotly[pheno$shape == '1'] <- 'circle-open'
pheno$shape_plotly[pheno$shape == '15'] <- 'square'
pheno$shape_plotly[pheno$shape == '16'] <- 'circle'


plot(pca$x, asp=1, col=pheno$color2, pch=pheno$shape, cex=2)


# legend
png("~/Desktop/legend.png", width = 600, height = 875, res = 200)
plot.new()
legend('topleft', legend = c('A','B','C'), title='Timepoint', col=c('green','red','blue'), pch=16, bty='n')
legend('bottomleft', legend = c('SC','MVS','SVS','EC','SS MMP8+', 'SS MMP8-'), title='Condition', col=c(1,1,1,'purple','brown4','orange'), pch=c(1,16,15,1,1), bty='n')
dev.off()


# 3D
require(plotly)
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3,
        color = pheno$color2, colors = sort(unique(pheno$color2)),
        symbol = ~pheno$shape_plotly, symbols = sort(unique(pheno$shape_plotly)))




### the one with surgical controls only
rvSC <- rowVars(cpm[, pheno$group == 'SC'])
rmSC <- rowMeans(cpm[, pheno$group == 'SC'])
filtSC <- which(rmSC > 1 & rvSC > 1)
phenoSC <- pheno[pheno$group == 'SC', ]
pcaSC <- BiocSingular::runPCA(t(log1p(cpm[filt, pheno$group == 'SC'])), rank=sum(pheno$group == 'SC'))

png("~/Desktop/pca1.png", width = 800, height = 700, res = 200)
plot(pcaSC$x, asp=1, col=phenoSC$color2, pch=phenoSC$shape, cex=1, las=1)
dev.off()

# loadings
write.csv(pcaSC$x[,1:2], file = '~/Desktop/loadings1.csv', quote = FALSE)



### the one with surgical controls and all vasoplegic patients
rvVS <- rowVars(cpm[, pheno$group %in% c('SC','MVS','SVS')])
rmVS <- rowMeans(cpm[, pheno$group %in% c('SC','MVS','SVS')])
filtVS <- which(rmVS > 1 & rvVS > 1)
phenoVS <- pheno[pheno$group %in% c('SC','MVS','SVS'), ]
pcaVS <- BiocSingular::runPCA(t(log1p(cpm[filt, pheno$group %in% c('SC','MVS','SVS')])),
                              rank=sum(pheno$group %in% c('SC','MVS','SVS')))

png("~/Desktop/pca2.png", width = 800, height = 700, res = 200)
plot(pcaVS$x, asp=1, col=phenoVS$color2, pch=phenoVS$shape, cex=1)
dev.off()

# loadings
write.csv(pcaVS$x[,1:2], file = '~/Desktop/loadings2.csv', quote = FALSE)



### the 3D one with all the patients
require(plotly)
plot_ly(data.frame(pca$x), x=~PC1, y=~PC2, z=~PC3,
        color = pheno$color2, colors = sort(unique(pheno$color2)),
        symbol = ~pheno$shape_plotly, symbols = sort(unique(pheno$shape_plotly)))




