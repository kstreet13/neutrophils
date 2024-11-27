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

