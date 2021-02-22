########################################################
###                                                  ###
###                  Piecewise SEM                   ###
###                                                  ###
###               CONVENTIONAL SEM's                 ###
###                                                  ###
########################################################

# Script by Koenraad Van Meerbeek and  Elisa Van Cleemput, 2019

#--------------------------------------------------------------------------------
# Clean workspace
rm(list=ls())

#######################################################################
################IMPORT AND EXPLORE#####################################
#######################################################################
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD")
dataset<-read.csv("Data_CWM_FD_PFTraits_Funct.csv", header=TRUE,sep=";")


#Check structure dataset
str(dataset)
head(dataset)

#Check normality and linear relationships
plot(dataset$SLA,dataset$logDB)
plot(dataset$LDMC,dataset$logDB)
plot(dataset$LNC_mass,dataset$logDB)
plot(dataset$Height,dataset$logDB)
plot(dataset$Invasion,dataset$logDB)
plot(dataset$FDis_1,dataset$logDB)

plot(dataset$SLA,dataset$logP)
plot(dataset$LDMC,dataset$logP)
plot(dataset$LNC_mass,dataset$logP)
plot(dataset$LCaMgC_mass_log,dataset$logP)
plot(dataset$Invasion,dataset$logP)
plot(dataset$FDis_2,dataset$logP)

plot(dataset$SLA,dataset$tea_S)
plot(dataset$LDMC,dataset$tea_S)
plot(dataset$LNC_mass,dataset$tea_S)
plot(dataset$LCaMgC_mass_log,dataset$tea_S)
plot(dataset$Invasion,dataset$tea_S)
plot(dataset$FDis_2,dataset$tea_S)

hist(dataset$logDB)
hist(dataset$logP)
hist(dataset$tea_S)

#######################################################################
###Construct SEM dataset
#Scale dataset
colnames(dataset)
scaleddata<-scale(dataset[,6:ncol(dataset)]) #standardize continuous variables (except for "Invasion")
semdata<-cbind(dataset[,c(1:5)],scaleddata)
semdata$Invasion <- semdata$Invasion/100 # set Invasion between 0 and 1 (std. estimates don't change because of this)


#Piecewise cannot yet handle interactions with categorical data
#Make dummy variables from variable 'species'
library(fastDummies)
semdata <- fastDummies::dummy_cols(semdata, select_columns = "Species")
semdata$Invsol<-semdata$Invasion*semdata$Species_Solidago
semdata$Invimp<-semdata$Invasion*semdata$Species_Impatiens

#######################################################################
################PIECEWISE SEM #########################################
#######################################################################
# Use the development version of piecewiseSEM instead of cran version because the latter contains bugs in Fisher's C test
# library(installr)
# uninstall.packages("piecewiseSEM")
# library(devtools)
# devtools::install_github("jslefche/piecewiseSEM@devel", build_vignette = TRUE) , force=T
library(piecewiseSEM) # version 2.1.0
library(nlme)
# vignette('piecewiseSEM')

#SEM logDB with species
psem2<- psem(
  lme(SLA ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LDMC ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LNC_mass ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(Height ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(FDis_1 ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(logDB ~ Species_Solidago + Invsol + Invimp + FDis_1 + SLA + LDMC + LNC_mass + Height, random = ~ 1 |Site, na.action = na.omit, data=semdata),
 
  LDMC %~~% SLA,
  LNC_mass %~~% SLA,
  Height %~~% SLA,
  FDis_1 %~~% SLA,
  LNC_mass %~~% LDMC,
  Height %~~% LDMC,
  FDis_1 %~~% LDMC,
  Height %~~% LNC_mass,
  FDis_1 %~~% Height
)
summary(psem2, intercepts = T,.progressBar = F)


#SEM tea_S with species
plot(semdata$LDMC,semdata$FDis_2)
cor(semdata$SLA,semdata$Height)
plot(semdata$SLA,semdata$tea_S)
cor(semdata$SLA,semdata$LNC_mass)

psem4<- psem(
  lme(SLA ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LDMC ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LNC_mass ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LCaMgC_mass_log ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(FDis_2 ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(tea_S ~ Species_Solidago + Invsol + Invimp + FDis_2 + SLA + LDMC + LNC_mass + LCaMgC_mass_log, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  
  LDMC %~~% SLA,
  LNC_mass %~~% SLA,
  LCaMgC_mass_log %~~% SLA,
  LNC_mass %~~% LDMC,
  LCaMgC_mass_log %~~% LDMC,
  FDis_2 %~~% LDMC,
  LCaMgC_mass_log %~~% LNC_mass,
  # # FDis_2 %~~% LNC_mass, # not significant in functional model, but significant in optical model
  FDis_2 %~~% LCaMgC_mass_log
)
summary(psem4, intercepts = F,.progressBar = F)

#SEM Polsen with species
plot(semdata$LDMC,semdata$logP)
cor(semdata$SLA,semdata$LNC_mass)

psem6<- psem(
  lme(SLA ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LDMC ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LNC_mass ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(LCaMgC_mass_log ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(FDis_2 ~ Species_Solidago + Invsol + Invimp, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  lme(logP ~ Species_Solidago + Invsol + Invimp + FDis_2 + SLA + LDMC + LNC_mass + LCaMgC_mass_log, random = ~ 1 |Site, na.action = na.omit, data=semdata),
  
  LDMC %~~% SLA,
  LNC_mass %~~% SLA,
  LCaMgC_mass_log %~~% SLA,
  LNC_mass %~~% LDMC,
  LCaMgC_mass_log %~~% LDMC,
  FDis_2 %~~% LDMC,
  LCaMgC_mass_log %~~% LNC_mass,
  # FDis_2 %~~% LNC_mass, # not significant in functional model, but significant in optical model
  FDis_2 %~~% LCaMgC_mass_log
)
summary(psem6, intercepts = F,.progressBar = F)


setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/SEM")
write.csv(summary(psem2, intercepts = F,.progressBar = F)$coefficients, file = "SEM.PFT.logDB.withSpecies.csv")
write.csv(summary(psem4, intercepts = F,.progressBar = F)$coefficients, file = "SEM.PFT.teaS.withSpecies.csv")
write.csv(summary(psem6, intercepts = F,.progressBar = F)$coefficients, file = "SEM.PFT.logP.withSpecies.csv")
