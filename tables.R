library(readxl)
data <- read_excel("~/Downloads/data.xlsx")
View(data)

###########
# Table 1 #
###########

data$tab1group <- NA
data$tab1group[data$Key == 'SC'] <- 'SC'
data$tab1group[data$Key %in% c('SVS','MVS')] <- 'VS'

# continuous variables
contvars <- c('Age','BMI','Cardiopulmonary bypass duration (minutes)',
              'Aortic crossclamp time (minutes)','post-CPB CI (L/min/m2)',
              'post-CPB SVRI  (dyne*sec)/cm5',
              'post-CPB Vasoactive-Inotropic Score','24h  CI (L/min/m2)',
              '24h  SVRI  (dyne*sec)/cm5','24h  Vasoactive-Inotropic Score')

all(contvars %in% names(data))

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][data$tab1group=='SC'],
                      data[[var]][data$tab1group=='VS']))
}

# binary variable
tab <- table(data$Sex, data$tab1group)
fisher.test(tab)


# Questions
# Preop SOFA (guessing = "SOFA score on enrollment"?)
wilcox.test(data$`SOFA score on enrollment`[data$tab1group=='SC'],
            data$`SOFA score on enrollment`[data$tab1group=='VS'])
# NOT ENOUGH DATA

# 0hr Post-CPB
## SOFA ???

# 24hr Post-CPB
## SOFA ???



###########
# Table 2 #
###########

data$tab2group <- NA
data$tab2group[data$Key == 'HC'] <- 'HC'
data$tab2group[data$Key %in% c('MMP8+','MMP8-')] <- 'SS'

# continuous variables
# VIS = Post-CPB
contvars <- c('Age','BMI','post-CPB Vasoactive-Inotropic Score')

all(contvars %in% names(data))

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][data$tab2group=='HC'],
                      data[[var]][data$tab2group=='SS']))
}

# binary variable
tab <- table(data$Sex, data$tab2group)
fisher.test(tab)


# Questions

# Baseline SOFA Score (if known) (0 if previously healthy):
# SOFA score on enrollment
# NOT ENOUGH DATA








###########
# Table 3 #
###########

contvars <- c('Age','BMI','Cardiopulmonary bypass duration (minutes)',
              'Aortic crossclamp time (minutes)','post-CPB CI (L/min/m2)',
              'post-CPB SVRI  (dyne*sec)/cm5',
              'post-CPB Vasoactive-Inotropic Score','24h  CI (L/min/m2)',
              '24h  SVRI  (dyne*sec)/cm5','24h  Vasoactive-Inotropic Score')
all(contvars %in% names(data))

# binary variables
binvars <- c('Sex','Preoperative ACE inhibitor/ARB  use',
             'Periop BetaBlocker','Systemic steroid  administration periop',
             'Modified ultrafiltration','Major Complications (1=yes, 0=no)')
all(binvars %in% names(data))

# Questions
# SOFA
# "Major complications per patient" change to "Major complications"
# ICU-free days
# Hospital-free days
# 60-day in-hospital mortality

# Column 1: SC vs MVS
#####################
data$t3g1 <- NA
data$t3g1[data$Key == 'SC'] <- 'SC'
data$t3g1[data$Key == 'MVS'] <- 'MVS'

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][data$t3g1=='SC'],
                      data[[var]][data$t3g1=='MVS']))
}
for(var in binvars){
    print(var)
    tab <- table(data[[var]], data$t3g1)
    print(fisher.test(tab))
}


# Column 2: SC vs SVS
#####################
data$t3g2 <- NA
data$t3g2[data$Key == 'SC'] <- 'SC'
data$t3g2[data$Key == 'SVS'] <- 'SVS'

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][data$t3g2=='SC'],
                      data[[var]][data$t3g2=='SVS'])$p.value)
}
for(var in binvars){
    print(var)
    tab <- table(data[[var]], data$t3g2)
    print(fisher.test(tab)$p.value)
}


# Column 3: MVS vs SVS
######################
data$t3g3 <- NA
data$t3g3[data$Key == 'MVS'] <- 'MVS'
data$t3g3[data$Key == 'SVS'] <- 'SVS'

for(var in contvars){
    print(var)
    print(wilcox.test(data[[var]][data$t3g3=='MVS'],
                      data[[var]][data$t3g3=='SVS'])$p.value)
}
for(var in binvars){
    print(var)
    tab <- table(data[[var]], data$t3g3)
    print(fisher.test(tab)$p.value)
}


