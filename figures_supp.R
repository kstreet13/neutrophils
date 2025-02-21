source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

heatdat <- as.matrix(cpm)
# force ordering
heatdat <- heatdat[,  c(grep('EC', colnames(cpm)),
                        grep('SC', colnames(cpm)),
                        grep('SS', colnames(cpm)),
                        grep('MVS', colnames(cpm)),
                        grep('SVS', colnames(cpm)))]
subpheno <- pheno[match(colnames(heatdat), pheno$sample), ]
subpheno$annTime <- subpheno$timepoint
subpheno$annTime[!subpheno$group %in% c('EC','SC','SS')] <- paste0(subpheno$annTime[!subpheno$group %in% c('EC','SC','SS')],'_VS')
subpheno$annTime[subpheno$group == 'SS'] <- paste0(subpheno$annTime[subpheno$group == 'SS'],'_SS')

require(NMF)
anncolors = list(
    annTime = c(`0` = "blue", `0_SS` = "orange",
                A = "green", B = "red", C = "purple",
                A_VS = "forestgreen", B_VS = "firebrick", C_VS = "slateblue4"))


aheatmap(log1p(heatdat), 
         annCol = subpheno[,c("group", "annTime"),drop=FALSE],
         Colv = NA, scale = "row",
         Rowv = FALSE,
         labRow = NA, labCol = NA,
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         width = 8.5,
         height = 20,
         filename = '~/Desktop/heat_all.png')
# Expression values are normalized by dividing each count by the total for that sample and multiplying by 1 million (Counts Per Million, CPM), then adding a pseudocount of 1 and taking the natural log, to control outliers. For the heatmap, each gene's normalized expression is then standardized to a mean of 0 and standard deviation of 1 to show differences between timepoints.



###################
# similar heatmap, but using restricted gene list (provided by Anaar) and including SS patients (split by MMP8)
###################
require(readxl)
genelist <- read_excel('data/SVS MMP8 Final Model Lists.xlsx', col_names = FALSE)
genelist <- genelist[,1, drop=TRUE]

heatdat <- as.matrix(cpm[genelist, ])
# force ordering
heatdat <- heatdat[,  c(grep('SC', colnames(cpm)),
                        grep('MVS', colnames(cpm)),
                        grep('SVS', colnames(cpm)),
                        grep('SS', colnames(cpm)))]
subpheno <- pheno[match(colnames(heatdat), pheno$sample), ]
subpheno$annTime <- paste0(subpheno$group, '_', subpheno$timepoint)
subpheno$annTime[subpheno$group == 'SS'] <- paste0(subpheno$group, '_', subpheno$mmp8)[subpheno$group == 'SS']
subpheno$annTime <- gsub('SVS_', 'VS_', subpheno$annTime)
subpheno$annTime <- gsub('MVS_', 'VS_', subpheno$annTime)

# reorder
for(g in c('SC','MVS','SVS')){
    ord <- order(subpheno$timepoint[subpheno$group==g])
    heatdat[, subpheno$group == g] <- heatdat[, subpheno$group == g][, ord]
    subpheno[subpheno$group == g, ] <- subpheno[subpheno$group == g, ][ord, ]
}
ord <- order(subpheno$mmp8[subpheno$group=='SS'])
heatdat[, subpheno$group == 'SS'] <- heatdat[, subpheno$group == 'SS'][, ord]
subpheno[subpheno$group == 'SS', ] <- subpheno[subpheno$group == 'SS', ][ord, ]


require(NMF)
anncolors = list(
    annTime = c(`SS_MMP8+` = "brown4", `SS_MMP8-`= "orange",
                SC_A = "green", SC_B = "red", SC_C = "purple",
                VS_A = "forestgreen", VS_B = "firebrick", VS_C = "slateblue4"))


aheatmap(log1p(heatdat), 
         annCol = subpheno[,c("group", "annTime"),drop=FALSE],
         Colv = NA, scale = "row",
         labRow = NA, labCol = colnames(heatdat),
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         width = 8.5,
         height = 20,
         filename = '~/Desktop/heat_new.png')



###################
# smaller heatmap with averages
###################
rm(list=ls())
source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

require(readxl)
genelist <- read_excel('data/MMP8+ SS.xlsx', col_names = FALSE)
genelist <- genelist[,1, drop=TRUE]

heatdat <- log1p(as.matrix(cpm[genelist, ]))
heatavg <- rowMeans(heatdat[, pheno$group == 'EC'])

for(g in c('SC','MVS','SVS')){
    for(t in c('A','B','C')){
        heatavg <- cbind(heatavg, rowMeans(heatdat[, pheno$group == g & pheno$timepoint == t]))
        colnames(heatavg)[ncol(heatavg)] <- paste0(g,'_',t)
    }
}
colnames(heatavg)[1] <- 'EC'

heatavg <- cbind(heatavg, rowMeans(heatdat[, which(pheno$mmp8 == "MMP8-")]))
colnames(heatavg)[ncol(heatavg)] <- 'SS_MMP8-'
heatavg <- cbind(heatavg, rowMeans(heatdat[, which(pheno$mmp8 == "MMP8+")]))
colnames(heatavg)[ncol(heatavg)] <- 'SS_MMP8+'

heatavg <- heatavg[, c(2:ncol(heatavg),1)]


subpheno <- data.frame(annTime = colnames(heatavg))
subpheno$annTime <- gsub('SVS_', 'VS_', subpheno$annTime)
subpheno$annTime <- gsub('MVS_', 'VS_', subpheno$annTime)
subpheno$annTime <- factor(subpheno$annTime) # has to be a factor, for some reason (hours!)

require(NMF)
anncolors = list(
    annTime = c(SC_A = "green", SC_B = "red", SC_C = "purple",
                VS_A = "forestgreen", VS_B = "firebrick", VS_C = "slateblue4",
                `SS_MMP8+` = "brown4", `SS_MMP8-` = "orange",
                EC = "blue"))
    

aheatmap(heatavg, 
         annCol = subpheno,
         Colv = NA, scale = "row",
         labRow = NA, labCol = colnames(heatavg),
         color = "RdBu:100",
         distfun = "pearson",
         annColors = anncolors,
         width = 8.5,
         height = 20,
         filename = '~/Desktop/heat_avg.png')




###################
# line drawing of specific genes in SVS patients
###################
rm(list=ls())
source('setup.R')
require(readxl)
genelist <- read_excel('data/SVS MMP8 Final Model Lists.xlsx', col_names = FALSE)
genelist <- genelist[,1, drop=TRUE]

cpm <- t(t(1e6*counts)/colSums(counts))

linedat <- as.matrix(cpm[genelist, pheno$group == 'SVS'])
linedat <- linedat[,  c(paste0('SVS',1:4,'A'),paste0('SVS',1:4,'B'),paste0('SVS',1:4,'C'))]
zscores <- (linedat - rowMeans(linedat)) / rowSds(linedat)
means <- cbind(
    rowMeans(zscores[,paste0('SVS',1:4,'A')]),
    rowMeans(zscores[,paste0('SVS',1:4,'B')]),
    rowMeans(zscores[,paste0('SVS',1:4,'C')])
)
# separated, without colors
cc <- rep('black',4)
layout(1)
plot(c(1,3), range(means), col='white', axes=FALSE,
     main = 'SVS samples; SVS MMPM8 genes', xlab = 'Timepoint', ylab = 'Avg. Expression (z-score)')
box(); axis(2); axis(1, at=1:3, labels=LETTERS[1:3])
for(g in genelist){
    lines(1:3, means[g,], col=alpha(cc[1], alpha=.2))
}





