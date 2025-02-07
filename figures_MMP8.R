source('setup.R')
cpm <- t(t(1e6*counts)/colSums(counts))

# mmp8 expression levels by timepoint
# split SS into mmp8+/-
mmp8 <- log1p(cpm["MMP8",])

mns <- sapply(unique(pheno$group), function(gp){
    sapply(unique(pheno$timepoint), function(tp){
        ind <- which(pheno$group == gp & pheno$timepoint == tp)
        if(length(ind) > 0){
            return(mean(mmp8[ind]))
        }else{
            return(NA)
        }
    })
})
sds <- sapply(unique(pheno$group), function(gp){
    sapply(unique(pheno$timepoint), function(tp){
        ind <- which(pheno$group == gp & pheno$timepoint == tp)
        if(length(ind) > 0){
            return(sd(mmp8[ind]))
        }else{
            return(NA)
        }
    })
})
mns <- mns[,-which(colnames(mns)=='SS')]
mns <- cbind(mns, mmp8p = c(mean(mmp8[which(pheno$mmp8=='MMP8+')]),NA,NA,NA))
mns <- cbind(mns, mmp8n = c(mean(mmp8[which(pheno$mmp8=='MMP8-')]),NA,NA,NA))
sds <- sds[,-which(colnames(sds)=='SS')]
sds <- cbind(sds, mmp8p = c(sd(mmp8[which(pheno$mmp8=='MMP8+')]),NA,NA,NA))
sds <- cbind(sds, mmp8n = c(sd(mmp8[which(pheno$mmp8=='MMP8-')]),NA,NA,NA))
lwr <- mns - 1.96*sds
lwr[lwr < 0] <- 0
upr <- mns + 1.96*sds


png(filename = '~/Desktop/mmp8expr.png', width = 1000, height = 600, res = 150)
layout(matrix(1:2, nrow=1))

plot(c(1,3), c(0,max(upr,na.rm=TRUE)), col='white',
     xlab='Timepoint', ylab='Avg log(CPM+1)', main='MMP8 Expression')
# EC (HC)
rect(1,lwr[1,'EC'],3,upr[1,'EC'], col=rgb(0,0,0,.1), border = NA)
lines(c(1,3),mns[c(1,1),'EC'], lwd=2)
# SS - MMP8+
rect(1,lwr[1,'mmp8p'],3,upr[1,'mmp8p'], col=rgb(.1,.1,1,.1), border = NA)
lines(c(1,3),mns[c(1,1),'mmp8p'], lwd=2, col=4)
# SS - MMP8-
rect(1,lwr[1,'mmp8n'],3,upr[1,'mmp8n'], col=rgb(1,.1,.1,.1), border = NA)
lines(c(1,3),mns[c(1,1),'mmp8n'], lwd=2, col=2)
# SVS
lines(1:3, mns[2:4,'SVS'], lwd=2, lty=2, col='grey50')
arrows(1:3, lwr[2:4,'SVS'], 1:3, upr[2:4,'SVS'], col='grey50',
       lty = 1, length = .05, angle = 90, code = 3)
# MVS
lines(1:3, mns[2:4,'MVS'], lwd=2, lty=3, col='grey75')
arrows(1:3, lwr[2:4,'MVS'], 1:3, upr[2:4,'MVS'], col='grey75',
       lty = 1, length = .05, angle = 90, code = 3)
# SC
lines(1:3, mns[2:4,'SC'], lwd=2)
arrows(1:3, lwr[2:4,'SC'], 1:3, upr[2:4,'SC'],
       lty = 1, length = .05, angle = 90, code = 3)


par(mar=c(0,0,0,0))
plot.new()

legend('left', legend = c('MMP8+ SS','MMP8- SS','HC','SVS','MVS','SC'),
       lwd = 10, col = c(rgb(.1,.1,1,.1),rgb(1,.1,.1,.1),rgb(0,0,0,.1),0,0,0))
legend('left', legend = c('MMP8+ SS','MMP8- SS','HC','SVS','MVS','SC'),
       lty = c(1,1,1,2,3,1), lwd=2, col = c(4,2,1,'grey50','grey75',1), bty='n')

dev.off()
par(mar=c(5,4,4,2)+.1)






# Gene lists for All sepsis, MMP8+ sepsis, and MMP8- sepsis with increased
# versus decreased genes relative to healthy controls, for each clinical
# category (6 lists total:  All sepsis up, all sepsis down, MMP8+ up, MMP8+
# down, MMP8-up, MMP8-down)




