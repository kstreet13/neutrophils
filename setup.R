set.seed(1)
require(readxl)
require(Matrix)
counts <- as.matrix(read_excel('data/Raw Counts.xlsx', sheet=1))
rownames(counts) <- counts[,1]
counts <- counts[,-1]

# exclude patient 32 (for now)
counts <- counts[, -which(colnames(counts) %in% c('32A','32B','32C'))]

# correct mislabeling
# patient 6 should be EC3, is labelled as SC3, which isn't a sample
counts[1,][counts[1,] == 'SC3'] <- 'EC3'

ptID <- colnames(counts)
samp <- counts[1,]
counts <- counts[-1,]
mode(counts) <- 'numeric'
counts <- Matrix(counts, sparse = TRUE)
colnames(counts) <- samp

# remove rows with 0 counts
counts <- counts[rowSums(counts) > 0, ]

# pheno data
# Clinical Group:	
# EC	Emergency Department Control (Healthy Control, no inflammation)
# SS	Septic Shock (SS1-6 are MMP8+ and SS7-12 are MMP8-)
# SC	Surgical Control (no evidence of vasoplegic syndrome)
# SVS	Severe vasoplegic syndrome
# MVS	Mild vasoplegic syndrome
# Timepoint:	
# A	Baseline
# B	0 hrs post-cardiopulmonary bypass (CPB)
# C	24 hrs post-cardiopulmonary bypass (CPB)
group <- gsub("([[:digit:]]+)([[:alnum:]]*)$", '',samp)
timepoint <- gsub("([[:alpha:]]+)([[:digit:]]+)", '',samp)
timepoint[timepoint==''] <- "0"
patient <- gsub("([[:alpha:]]*)$", '',samp)
mmp8 <- rep(NA, ncol(counts))
mmp8[patient %in% paste0('SS',1:6)] <- 'MMP8+'
mmp8[patient %in% paste0('SS',7:12)] <- 'MMP8-'

pheno <- data.frame(
    sample = samp,
    group = group,
    timepoint = timepoint,
    patient = patient,
    patientID = ptID,
    mmp8 = mmp8
)
rownames(pheno) <- samp

rm(group, timepoint, patient, mmp8, samp, ptID)

# color
# Blue: ED (ambulatory, healthy) controls
# Orange: septic shock
# Green: surgical control, Preop time point ("A" time point)
# Forestgreen: Vasoplegic syndrome (mild and severe), Preop timepoint ("A" time point)
# Red: surgical control, 0hr Post-CPB time point ("B" timepoint)
# Firebrick: Vasoplegic syndrome (mild and severe), 0hr Post-CPB timepoint ("B" time point)
# Purple: surgical control, 24hr Post-CPB time point ("C" timepoint)
# Slateblue4:  Vasoplegic syndrome (mild and severe), 24hr Post-CPB timepoint ("C" time point)
require(dplyr)

pheno$color <- case_when(
    pheno$group == 'EC' ~ 'blue',
    pheno$group == 'SS' ~ 'orange',
    pheno$group == 'SC' & pheno$timepoint =='A' ~ 'green',
    pheno$group %in% c('MVS','SVS') & pheno$timepoint =='A' ~ 'forestgreen',
    pheno$group == 'SC' & pheno$timepoint =='B' ~ 'red',
    pheno$group %in% c('MVS','SVS') & pheno$timepoint =='B' ~ 'firebrick',
    pheno$group == 'SC' & pheno$timepoint =='C' ~ 'purple',
    pheno$group %in% c('MVS','SVS') & pheno$timepoint =='C' ~ 'slateblue4'
)

pheno <- pheno[order(pheno$sample), ]
counts <- counts[ ,order(colnames(counts))]

# check
all(pheno$sample == colnames(counts))
