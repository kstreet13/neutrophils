library(readxl)
data <- read_excel("~/Downloads/data.xlsx")
#View(data)
data$ACEI_bin <- as.numeric(grepl('ACE/ARB use within', data$`Preoperative ACE inhibitor/ARB  use`))

sofa <- read_excel("~/Downloads/sofa for surgery pts.xlsx")
data <- merge(data, sofa, all.x = TRUE)
rm(sofa)

complications <- read_excel("~/Downloads/complications.xlsx")
complications <- complications[-13,]
data <- merge(data, complications, by = 'ï»¿\"Record ID\"', all.x = TRUE)
rm(complications)

###########
# Table 1 #
###########

data$tab1group <- NA
data$tab1group[data$Key == 'SC'] <- 'SC'
data$tab1group[data$Key %in% c('SVS','MVS')] <- 'VS'

# continuous variables
contvars <- c('Age','BMI','Cardiopulmonary bypass duration (minutes)',
              'Aortic crossclamp time (minutes)','sofa_baseline',
              'sofa_postcpb','post-CPB CI (L/min/m2)',
              'post-CPB SVRI  (dyne*sec)/cm5',
              'post-CPB Vasoactive-Inotropic Score','sofa_day1',
              '24h  CI (L/min/m2)','24h  SVRI  (dyne*sec)/cm5',
              '24h  Vasoactive-Inotropic Score')

all(contvars %in% names(data))

for(var in contvars){
    print(var)
    print(summary(data[[var]][which(data$tab1group=='SC')]))
    print(summary(data[[var]][which(data$tab1group=='VS')]))
    print(wilcox.test(data[[var]][which(data$tab1group=='SC')],
                      data[[var]][which(data$tab1group=='VS')]))
}

# binary variable
tab <- table(data$Sex, data$tab1group)
fisher.test(tab)



###########
# Table 2 #
###########

data$tab2group <- NA
data$tab2group[data$Key == 'HC'] <- 'HC'
data$tab2group[data$Key %in% c('MMP8+','MMP8-')] <- 'SS'

# ASSUMPTION: all missing SOFA scores are 0
data$sofa_baseline[is.na(data$sofa_baseline)] <- 0
data$sofa_postcpb[is.na(data$sofa_postcpb)] <- 0
data$sofa_day1[is.na(data$sofa_day1)] <- 0
data$`SOFA score on enrollment`[is.na(data$`SOFA score on enrollment`)] <- 0

# continuous variables
# VIS = Post-CPB
contvars <- c('Age','BMI','sofa_baseline',
              'SOFA score on enrollment',
              'post-CPB Vasoactive-Inotropic Score')

all(contvars %in% names(data))

for(var in contvars){
    print(var)
    print(summary(data[[var]][which(data$tab2group=='HC')]))
    print(summary(data[[var]][which(data$tab2group=='SS')]))
    print(wilcox.test(data[[var]][which(data$tab2group=='HC')],
                      data[[var]][which(data$tab2group=='SS')]))
}

# binary variable
tab <- table(data$Sex, data$tab2group)
fisher.test(tab)


# Questions

# Baseline SOFA IQRs don't match








###########
# Table 3 #
###########

contvars <- c('Age','BMI','Cardiopulmonary bypass duration (minutes)',
              'Aortic crossclamp time (minutes)','sofa_baseline',
              'sofa_postcpb','post-CPB CI (L/min/m2)',
              'post-CPB SVRI  (dyne*sec)/cm5',
              'post-CPB Vasoactive-Inotropic Score',
              'sofa_day1','24h  CI (L/min/m2)',
              '24h  SVRI  (dyne*sec)/cm5','24h  Vasoactive-Inotropic Score',
              'Total number of complications','ICU free days',
              'Hospital fee days')
all(contvars %in% names(data))

# binary variables
binvars <- c('Sex','ACEI_bin',
             'Periop BetaBlocker','Systemic steroid  administration periop',
             'Modified ultrafiltration','Major Complications (1=yes, 0=no)',
             '60-day in-hospital mortality.y')
all(binvars %in% names(data))

myWilcoxTest <- function(x,y){
    r1 <- format(summary(x)[c(3,2,5)], digits=4)
    r2 <- format(summary(y)[c(3,2,5)], digits=4)
    pval <- wilcox.test(x,y)$p.value
    paste0(r1[1],' (',r1[2],',',r1[3],')  ,  ',
           r2[1],' (',r2[2],',',r2[3],')  ,  ', 
           format(pval, digits=5, scientific=FALSE))
}


# Column 1: SC vs MVS
#####################
data$t3g1 <- NA
data$t3g1[data$Key == 'SC'] <- 'SC'
data$t3g1[data$Key == 'MVS'] <- 'MVS'

for(var in contvars){
    print(var)
    print(summary(data[[var]][which(data$t3g1=='SC')]))
    print(summary(data[[var]][which(data$t3g1=='MVS')]))
    print(wilcox.test(data[[var]][which(data$t3g1=='SC')],
                      data[[var]][which(data$t3g1=='MVS')]))
}
sapply(contvars, function(var){
    myWilcoxTest(data[[var]][which(data$t3g1=='SC')],
                 data[[var]][which(data$t3g1=='MVS')])
})

for(var in binvars){
    print(var)
    tab <- table(data[[var]], data$t3g1)
    print(tab)
    print(fisher.test(tab))
}


# Column 2: SC vs SVS
#####################
data$t3g2 <- NA
data$t3g2[data$Key == 'SC'] <- 'SC'
data$t3g2[data$Key == 'SVS'] <- 'SVS'

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][which(data$t3g2=='SC')],
                      data[[var]][which(data$t3g2=='SVS')])$p.value)
}
sapply(contvars, function(var){
    myWilcoxTest(data[[var]][which(data$t3g2=='SC')],
                 data[[var]][which(data$t3g2=='SVS')])
})

for(var in binvars){
    print(var)
    tab <- table(data[[var]], data$t3g2)
    print(tab)
    print(fisher.test(tab)$p.value)
}


# Column 3: MVS vs SVS
######################
data$t3g3 <- NA
data$t3g3[data$Key == 'MVS'] <- 'MVS'
data$t3g3[data$Key == 'SVS'] <- 'SVS'

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][which(data$t3g3=='MVS')],
                      data[[var]][which(data$t3g3=='SVS')])$p.value)
}
sapply(contvars, function(var){
    myWilcoxTest(data[[var]][which(data$t3g3=='MVS')],
                 data[[var]][which(data$t3g2=='SVS')])
})

for(var in binvars){
    print(var)
    tab <- table(data[[var]], data$t3g3)
    print(tab)
    print(fisher.test(tab)$p.value)
}


