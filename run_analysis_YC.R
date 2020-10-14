
# SET/CHANGE your working directory to you current working directory
setwd('/data/NHLBI_BCB/Arai_Lab/PGRScorer')

# read input data
########################################################################
#
# Create one input file (space delimited) with the following columns/variables:
# 
# 1. Individual ID
# 2. GRS (calculated using calculate_PRS.sh from https://github.com/seifudd/PGRScorer)
# 3. age
# 4. sex (Male = 0, Female = 1)
# 5. CAC
#
# replace FILENAME with name of input file (space delimited)
#
########################################################################
pheno = read.table(file="FILENAME", header=T, row.names=1)

########################################################################
#
# Creating dummy data here to run/test the models - START
#
# COMMENT OUT or DELETE when reading real input data above
#
########################################################################
n.sample <- 500
pheno <- data.frame(
  GRS = rnorm(n.sample, 0, 50),
  age = sample(20:70, n.sample, replace = T),
  sex = sample(0:1, n.sample, replace = T)
)

pheno$CAC <- pheno$GRS + exp(rnorm(n.sample, 5, 1))
pheno$CAC[which(pheno$CAC < 0)] <- 0
########################################################################
#
# Creating dummy data here to run/test the models - END
#
# COMMENT OUT or DELETE when reading real input data above
#
########################################################################

# logtransform the CAC variable
pheno$CAClog = log(pheno$CAC + 1)

# stratify CAC in categories based on clinical recommendations [ CACS=0, 0<CACS<100, 100<CACS<400, CACS >400 (ordered 0,1,2,3) ]
pheno$CACcateg <- 0
pheno$CACcateg[which(pheno$CAC > 0 & pheno$CAC <= 100)] <- 1
pheno$CACcateg[which(pheno$CAC > 100 & pheno$CAC <= 400)] <- 2
pheno$CACcateg[which(pheno$CAC > 400)] <- 3

# stratify the GRS into a quintile variable
pheno$GRSquint <- cut(pheno$GRS, quantile(pheno$GRS, seq(0, 1, 0.2)), labels = F, include.lowest = T)

########################################################################
#
# Creating the models for All, Males and Females
#
#
########################################################################
for(type in c("All","Male","Female")){
  
  temp <- pheno
  if(type == 'Male'){
    temp <- temp[which(temp$sex == 0),]
  } else if(type == 'Female'){
    temp <- temp[which(temp$sex == 1),]
  }

# fit models
  # [CAC continous (after log transformation of CAC score +1)]~[GRS as continuous]
  f_GcCc <- if(type != 'All') 'CAClog ~ GRS + age' else 'CAClog ~ GRS + age + factor(sex)'

  # [CAC continous (after log transformation of CAC score +1)]~[GRS as quintile]
  f_GqCc <- if(type != 'All') 'CAClog ~ GRSquint + age' else 'CAClog ~ GRSquint + age + factor(sex)'

  # [CAC continous (after log transformation of CAC score +1)]~[GRS as quintile - dichotomized]
  f_GqDichCc <- if(type != 'All') 'CAClog ~ factor(GRSquint) + age' else 'CAClog ~ factor(GRSquint) + age + factor(sex)'

  # [CAC in categories based on clinical recommendations]~[GRS as continuous]
  f_GcCq <- if(type != 'All') 'CACcateg ~ GRS + age' else 'CACcateg ~ GRS + age + factor(sex)'

  # [CAC in categories based on clinical recommendations]~[GRS as quiltile]
  f_GqCq <- if(type != 'All') 'CACcateg ~ GRSquint + age' else 'CACcateg ~ GRSquint + age + factor(sex)'

  # [CAC in categories based on clinical recommendations]~[GRS as quiltile - dichotomized]
  f_GqDichCq <- if(type != 'All') 'CACcateg ~ factor(GRSquint) + age' else 'CACcateg ~ factor(GRSquint) + age + factor(sex)'
  

# run models above
  m_GcCc <- lm(formula(f_GcCc), temp)
  m_GqCc <- lm(formula(f_GqCc), temp)
  m_GqDichCc <- lm(formula(f_GqDichCc), temp[which(temp$GRSquint %in% c(1,5)),])
  
  m_GcCq <- lm(formula(f_GcCq), temp)
  m_GqCq <- lm(formula(f_GqCq), temp)
  m_GqDichCq <- lm(formula(f_GqDichCq), temp[which(temp$GRSquint %in% c(1,5)),])
  
# extract model results
  s_GcCc <- as.data.frame(coef(summary(m_GcCc)))
  s_GqCc <- as.data.frame(coef(summary(m_GqCc)))
  s_GqDichCc <- as.data.frame(coef(summary(m_GqDichCc)))
  s_GcCq <- as.data.frame(coef(summary(m_GcCq)))
  s_GqCq <- as.data.frame(coef(summary(m_GqCq)))
  s_GqDichCq <- as.data.frame(coef(summary(m_GqDichCq)))
  
# save model results
  write.csv(s_GcCc, paste0('GRScont_CACcont_',type,'.csv'), quote = F, row.names = T)
  write.csv(s_GqCc, paste0('GRSquint_CACcont_',type,'.csv'), quote = F, row.names = T)
  write.csv(s_GqDichCc, paste0('GRSdich_CACcont_',type,'.csv'), quote = F, row.names = T)
  write.csv(s_GcCq, paste0('GRScont_CACcateg_',type,'.csv'), quote = F, row.names = T)
  write.csv(s_GqCq, paste0('GRSquint_CACcateg_',type,'.csv'), quote = F, row.names = T)
  write.csv(s_GqDichCq, paste0('GRSdich_CACcateg_',type,'.csv'), quote = F, row.names = T)

}

##############################################################################################################
#
# Notes from Pat, Larry and MJ
#
##############################################################################################################

# 1. GRS as continuous
# 2. GRS in quintiles (to be consistent with the previous paper), 
# 3. and also run a model comparing the highest quintile to the lowest quintile for GRS - divide the GRS to high low - dichotomize
# 
# Then for CAC: 
# 4. CAC continuous (after log transformation of CAC score +1) and 
# 5. CAC in categories based on clinical recommendations - CACS=0, 0<CACS<100, 100<CACS<400, CACS >400 (ordered 0,1,2,3) - data$CAC
# 
# That would give us 6 models in total for everyone, 
# 6 models for women alone and 
# 6 models for men alone.  
# That would be 18 models in total.