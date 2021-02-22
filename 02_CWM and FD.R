########################################################
###                                                  ###
###         Calculate CWM and FD indices             ###
###                                                  ###
########################################################

# Script by Elisa Van Cleemput, 2017-2019

#--------------------------------------------------------------------------------
# Clean workspace
rm(list=ls())

#--------------------------------------------------------------------------------
# Load libraries
library(xlsx)
library(FD)
library(dplyr)
library(plyr) # ldply
library(qdapTools) # list_vect2df

#--------------------------------------------------------------------------------
# Basic functions

# Load conventioanlly and optically measured functional traits
data_PFTraits <- function(wd1,wd2,traitlist) {
  setwd(wd1)
  PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
  PFT$LNC_mass <- 10 * PFT$LNC # % dry mass to mg/g
  PFT$LCaC_mass <- 10 * PFT$Ca # % dry mass to mg/g
  PFT$LMgC_mass <- 10 * PFT$Mg # % dry mass to mg/g
  PFT$LCaMgC_mass <- PFT$LCaC_mass + PFT$LMgC_mass
  PFT$LCaC_mass_log <- log10(PFT$LCaC_mass)
  PFT$LMgC_mass_log <- log10(PFT$LMgC_mass)
  PFT$LCaMgC_mass_log <- PFT$LCaC_mass_log + PFT$LMgC_mass_log

  PFT <- PFT[,c("Species.code","Site",traitlist)]
  
  # Calculate average trait values per species
  PFT_average <- aggregate(PFT[,3:ncol(PFT)], list(PFT$Species.code), mean)
  rownames(PFT_average) <- paste(PFT_average$Group.1,"av",sep="_")
  PFT_average<-PFT_average[,-1]
  
  # Merge average trait data with site specific trait data
  rownames(PFT) <- paste(PFT$Species.code,PFT$Site,sep="_")
  PFT_specific <- PFT[, 3:ncol(PFT)]
  PFT_all<-rbind(PFT_average,PFT_specific)
  
  # Store IDs of trait records so that we can afterwards filter the vegetation survey matrix to available records
  setwd(wd2)
  write.table(rownames(PFT_all),file="traitrecords_conv.csv",sep=";",row.names = F)
  
  return(PFT_all)  
}
data_PFTraits_NAT <- function(wd1,wd2,traitlist) {
  setwd(wd1)
  PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
  PFT$LNC_mass <- 10 * PFT$LNC # % dry mass to mg/g
  PFT$LCaC_mass <- 10 * PFT$Ca # % dry mass to mg/g
  PFT$LMgC_mass <- 10 * PFT$Mg # % dry mass to mg/g
  PFT$LCaMgC_mass <- PFT$LCaC_mass + PFT$LMgC_mass
  PFT$LCaC_mass_log <- log10(PFT$LCaC_mass)
  PFT$LMgC_mass_log <- log10(PFT$LMgC_mass)
  PFT$LCaMgC_mass_log <- PFT$LCaC_mass_log + PFT$LMgC_mass_log

  PFT <- PFT[,c("Species.code","Site",traitlist)]
  
  # REMOVE ALL ROWS REFERRING TO S. GIGANTEA OR I. GLANDULIFERA
  PFT <- PFT[!(PFT$Species.code=="SOLGIG" | PFT$Species.code=="IMPGLA"),] # 6 times SOLGIG and 4 times IMPGLA
  
  # PFT_average <- aggregate(PFT[,which(names(PFT) == 'LCC'):ncol(PFT)], list(PFT$Species.code), mean)
  PFT_average <- aggregate(PFT[,3:ncol(PFT)], list(PFT$Species.code), mean)
  rownames(PFT_average) <- paste(PFT_average$Group.1,"av",sep="_")
  PFT_average<-PFT_average[,-1]
  
  # Merge average trait data with site specific trait data
  rownames(PFT) <- paste(PFT$Species.code,PFT$Site,sep="_")
  PFT_specific <- PFT[, 3:ncol(PFT)]
  PFT_all<-rbind(PFT_average,PFT_specific)
  
  # Store IDs of trait records so that we can afterwards filter the vegetation survey matrix to avaible records
  setwd(wd2)
  write.table(rownames(PFT_all),file="traitrecords_conv_NAT.csv",sep=";",row.names = F)
  
  return(PFT_all)  
}
data_POTraits <- function(wd1,wd2,traitlist,type) {
    setwd(wd1)
    PFT_pred_plsr <- read.csv(paste("PLSR_best_predictions_",".csv",sep=type),sep=";",header=TRUE, row.names=1)
    PFT_pred_plsr$LCaMgC_mass_log <- PFT_pred_plsr$LCaC_mass_log + PFT_pred_plsr$LMgC_mass_log
    
    PFT_pred_plsr_acc <- read.csv("PLSR_best_predictions_accuracy.csv",sep=";",header=TRUE, row.names=1)
    
    # Calculate average trait values per species
    PFT_pred_plsr$Species.code <- substr(rownames(PFT_pred_plsr), 1, 6)
    PFT_pred_plsr_average <- aggregate(PFT_pred_plsr[,1:which(names(PFT_pred_plsr) == 'Species.code')-1], list(PFT_pred_plsr$Species.code), mean)
    rownames(PFT_pred_plsr_average) <- paste(PFT_pred_plsr_average$Group.1,"av",sep="_")
    PFT_pred_plsr_average<-PFT_pred_plsr_average[,-c(1)]
    
    # Merge average trait data with site specific trait data
    PFT_pred_plsr_specific <- PFT_pred_plsr[, 1:which(names(PFT_pred_plsr) == 'Species.code')-1]
    rownames(PFT_pred_plsr_specific) <- substr(rownames(PFT_pred_plsr_specific), 1, nchar(rownames(PFT_pred_plsr_specific))-9)
    
    PFT_pred_plsr_all<-rbind(PFT_pred_plsr_average,PFT_pred_plsr_specific)
    PFT_pred_plsr_all <- PFT_pred_plsr_all[,traitlist]
    
    # Store IDs of trait records so that we can afterwards filter the vegetation survey matrix to avaible records
    setwd(wd2)
    write.table(rownames(PFT_pred_plsr_all),file="traitrecords_opt.csv",sep=";",row.names = F)

    output <- list(PFT_pred = PFT_pred_plsr_all,
                   PFT_pred_acc = PFT_pred_plsr_acc)

  return(output)
} # extend with CN, carot, Ca estimations etc.
data_POTraits_NAT <- function(wd1,wd2,traitlist,type) {
    setwd(wd1)
    PFT_pred_plsr <- read.csv(paste("PLSR_best_predictions_",".csv",sep=type),sep=";",header=TRUE, row.names=1)
    PFT_pred_plsr$LCaMgC_mass_log <- PFT_pred_plsr$LCaC_mass_log + PFT_pred_plsr$LMgC_mass_log
    
    PFT_pred_plsr_acc <- read.csv("PLSR_best_predictions_accuracy.csv",sep=";",header=TRUE, row.names=1)
    
    PFT_pred_plsr$Species.code <- substr(rownames(PFT_pred_plsr), 1, 6)
    
    # REMOVE ALL ROWS REFERRING TO S. GIGANTEA OR I. GLANDULIFERA
    PFT_pred_plsr <- PFT_pred_plsr[!(PFT_pred_plsr$Species.code=="SOLGIG" | PFT_pred_plsr$Species.code=="IMPGLA"),] # 6 times SOLGIG and 4 times IMPGLA
    
    # Calculate average trait values per species
    PFT_pred_plsr_average <- aggregate(PFT_pred_plsr[,1:which(names(PFT_pred_plsr) == 'Species.code')-1], list(PFT_pred_plsr$Species.code), mean)
    rownames(PFT_pred_plsr_average) <- paste(PFT_pred_plsr_average$Group.1,"av",sep="_")
    PFT_pred_plsr_average<-PFT_pred_plsr_average[,-c(1)]
    
    # Merge average trait data with site specific trait data
    PFT_pred_plsr_specific <- PFT_pred_plsr[, 1:which(names(PFT_pred_plsr) == 'Species.code')-1]
    rownames(PFT_pred_plsr_specific) <- substr(rownames(PFT_pred_plsr_specific), 1, nchar(rownames(PFT_pred_plsr_specific))-9)
    
    PFT_pred_plsr_all<-rbind(PFT_pred_plsr_average,PFT_pred_plsr_specific)
    PFT_pred_plsr_all <- PFT_pred_plsr_all[,traitlist]
    
    # Store IDs of trait records so that we can afterwards filter the vegetation survey matrix to avaible records
    setwd(wd2)
    write.table(rownames(PFT_pred_plsr_all),file="traitrecords_opt_NAT.csv",sep=";",row.names = F)
    
    output <- list(PFT_pred = PFT_pred_plsr_all,
                   PFT_pred_acc = PFT_pred_plsr_acc)
    
  return(output)
}

# Load vegetation survey (= species x plots)
data_plots <- function(wd3) {
  setwd(wd3)
  Plots<-read.xlsx("Species_x_plots.xlsx",4)
  rownames(Plots)<-Plots[,'Species_Code']
  Plots_meta <- Plots[, 1:which(names(Plots)=='Species_Code')]
  Plots <- Plots[, which(names(Plots)=='S1.1'):ncol(Plots)]
  return(Plots)
}
data_plots_NAT <- function(wd3) {
  setwd(wd3)
  Plots<-read.xlsx("Species_x_plots.xlsx",4)
  rownames(Plots)<-Plots[,'Species_Code']
  
  # REMOVE ALL ROWS REFERRING TO S. GIGANTEA OR I. GLANDULIFERA
  Plots <- Plots[!(Plots$Species_Code=="SOLGIG" | Plots$Species_Code=="IMPGLA"),] # 6 times SOLGIG and 4 times IMPGLA
  
  Plots_meta <- Plots[, 1:which(names(Plots)=='Species_Code')]
  Plots <- Plots[, which(names(Plots)=='S1.1'):ncol(Plots)]
  return(Plots)
}

# Create site-specific trait composition per plot
Create_compos <- function(Plots,site_name,plotsbegin, plotsend,specieslist) {
    Plots_site <- Plots %>% dplyr::select(plotsbegin:plotsend)
    Plots_site<-Plots_site[apply(Plots_site[], 1, function(x) !all(x==0)),]
    rownames(Plots_site)<-paste(rownames(Plots_site),site_name, sep = "_")
    
    for (i in 1:length(specieslist)){
      rownames(Plots_site) <- replace(rownames(Plots_site),
                                      rownames(Plots_site) == paste(specieslist[i],site_name,sep="_"),
                                      paste(substr(specieslist[i],1,6),"av",sep="_"))
    }
  return(Plots_site)
}

# Calculation of CWM (Lavorel et al. 2008) and other functional diversty indices of Villéger et al. (2008) and Laliberté & Legendre (2010)
# First, check which species in the vegetation survey have a trait measurement, preferably site-specific but average if possible
# Then, Prepare data for CWM calculation
# Lastly perform calculations
Community_values <- function(Plots_site, sp_data) {
  # the dbFD funtion has following requirements:
  # 'x' = data frame of functional traits (or reflectance) = sp_data 
  # 'a' = matrix with abundances of species: rows are sites and species are columns. --> we have to transpose Plots_site
  # The number and order of species in 'a' must match the number and order of species in 'x'. 
  Plots_site<-Plots_site[rownames(Plots_site) %in% rownames(sp_data),]
  sp_data_site<-sp_data[rownames(sp_data) %in% rownames(Plots_site),]
  
  # order alphabetically:
  Plots_site<-Plots_site[order(rownames(Plots_site)),]
  sp_data_site<-sp_data_site[order(rownames(sp_data_site)),]
  
  # There is 1 plot that has no trait data (S4.1), because too low veg cover 
  # dbFD won't work in that case -> delete this plots
  if (length(which(apply(Plots_site, 2, sum) == 0)) > 0 ) {
    Plots_site <- Plots_site[, - which(apply(Plots_site, 2, sum) == 0)]
  }
  
  Plots_site.t<-t(Plots_site)
  
  # We should restrict on interesting traits when calculating diversity. Maybe they doesn't make sense at all because we only measured most abundant spp.
  # FRic = functional richness = convex hull volume (Villéger et al. 2008)
  # FEve = functional eveness
  # Fdiv = functional divergence
  # FDis = functional dispersion: weighted average distance to centroid (Laliberté and Legendre 2010). For communities composed of only one species, dbFD returns a FDis value of 0.
  Calc<-dbFD(sp_data_site, Plots_site.t, w.abun=T, stand.x=T, calc.FRic =T, m="max", calc.CWM=T, calc.FDiv=T)
  
  return(Calc)
}

#--------------------------------------------------------------------------------
# Functions thath apply the basic functions to different datasets and save the results

Plots_CWM_FD <- function (wd1,wd2,wd3,approach,traitlist,type) {
  # load data
  if (approach == "PFTraits") {
    traits <- data_PFTraits(wd1,wd2,traitlist)
    Plots <- data_plots(wd3)
  } else if (approach == "PFTraits_NAT"){
    traits <- data_PFTraits_NAT(wd1,wd2,traitlist)
    Plots <- data_plots_NAT(wd3)
    # There is 1 plot containing only IAS and no other species (S8.2), dbFD won't work because no cover -> delete this plots
    Plots <- Plots[, - which(apply(Plots, 2, sum) == 0)] 
  } else if (approach == "POTraits"){
    traits <- data_POTraits(wd1,wd2,traitlist,type)$PFT_pred
    traits_acc <- data_POTraits(wd1,wd2,traitlist,type)$PFT_pred_acc
    Plots<- data_plots(wd3)
  } else if (approach == "POTraits_NAT"){
    traits <- data_POTraits_NAT(wd1,wd2,traitlist,type)$PFT_pred
    traits_acc <- data_POTraits_NAT(wd1,wd2,traitlist,type)$PFT_pred_acc
    Plots<- data_plots_NAT(wd3)
    Plots <- Plots[, - which(apply(Plots, 2, sum) == 0)] 
  }
  
  # Apply "Create_compos" to each site and replace records by average or another site-specific equivalent if necessary and available
  Plots_BattSD <- Create_compos(Plots,"BattSD","I5.1","I6.4",c("CALSEP","PHRAUS"))
  rownames(Plots_BattSD) <- replace(rownames(Plots_BattSD),rownames(Plots_BattSD) == "GLEHED_BattSD","GLEHED_BattSP")
  Plots_BattSP <- Create_compos(Plots,"BattSP","I8.1","I8.6",c("EPITET"))
  rownames(Plots_BattSP) <- replace(rownames(Plots_BattSP),rownames(Plots_BattSP) == "CIRARV_BattSP","CIRARV_BattSD")
  rownames(Plots_BattSP) <- replace(rownames(Plots_BattSP),rownames(Plots_BattSP) == "IMPGLA_BattSP","IMPGLA_BattSD")
  rownames(Plots_BattSP) <- replace(rownames(Plots_BattSP),rownames(Plots_BattSP) == "URTDIO_BattSP","URTDIO_BattSD")
  Plots_BattN <- Create_compos(Plots,"BattN","I3.1","I4.4",c())
  Plots_BattM <- Create_compos(Plots,"BattN","I7.1","I7.4",c("EPITET","FILULM"))
  Plots_DoBeI <- Create_compos(Plots,"DoBeI","I1.1","I2.4",c())
  Plots_DoBeS <- Create_compos(Plots,"DoBeS","S10.1","S10.4",c("CIRARV","CIROLE","CALSEP","EQUPAL","GLEHED","PHAARU","URTDIO"))
  Plots_DeMatCem <- Create_compos(Plots,"DeMatCem","S8.1","S8.4",c("EPITET","EQUARV","HYPPER","RANREP","URTDIO"))
  rownames(Plots_DeMatCem) <- replace(rownames(Plots_DeMatCem),rownames(Plots_DeMatCem) == "DACGLO_DeMatCem","DACGLO_DeMatSpo")
  rownames(Plots_DeMatCem) <- replace(rownames(Plots_DeMatCem),rownames(Plots_DeMatCem) == "EUPCAN_DeMatCem","EUPCAN_DeMatSpo")
  Plots_DeMatSpo <- Create_compos(Plots,"DeMatSpo","S6.1","S7.4",c("CIRARV","EPITET","HYPPER","RANREP","URTDIO"))
  rownames(Plots_DeMatSpo) <- replace(rownames(Plots_DeMatSpo),rownames(Plots_DeMatSpo) == "HOLLAN_DeMatSpo","HOLLAN_DeMatCem")
  Plots_HaachtBev <- Create_compos(Plots,"HaachtBev","S11.1","S11.4",c("ARRELA","CIRARV","EPIHIR","EPITET","HYPPER","RANREP"))
  Plots_HaachtHol <- Create_compos(Plots,"HaachtHol","S1.1","S2.4",c("ANGSYL","CIRARV","EPITET","EQUPAL","URTDIO"))
  Plots_HeverleeKapel <- Create_compos(Plots,"HeverleeKapel","S9.1","S9.4",c("CIRARV","CIROLE","CALSEP","EPIHIR","EPITET","EQUARV","EUPCAN","GLEHED","RUBFRU","TANVUL","URTDIO"))
  Plots_PapenIn <- Create_compos(Plots,"PapenIn","I9.1","I9.4",c("ANGSYL","CIRARV","EQUPAL","FILULM","GLEHED","HOLLAN","RUBFRU","TANVUL","URTDIO"))
  Plots_PapenOut <- Create_compos(Plots,"PapenIn","I10.1","I10.4",c("ANGSYL","CIRARV","EQUPAL","FILULM","GLEHED","HOLLAN","RUBFRU","TANVUL","URTDIO"))
  Plots_HaachtSchor <- Create_compos(Plots,"av","S3.1","S4.4",c())
  Plots_Diepenbeek <- Create_compos(Plots,"av","S5.1","S5.4",c())
  
  # Calculate CWM and functional diversty indices 
  Plots_list <- list(BattSD=Plots_BattSD,BattSP=Plots_BattSP,BattN=Plots_BattN,BattM=Plots_BattM,
                     DoBeI=Plots_DoBeI,DoBeS=Plots_DoBeS,
                     DeMatCem=Plots_DeMatCem,DeMatSpo=Plots_DeMatSpo,
                     HaachtBev=Plots_HaachtBev,HaachtHol=Plots_HaachtHol,
                     HeverleeKapel=Plots_HeverleeKapel,PapenIn=Plots_PapenIn,PapenOut=Plots_PapenOut,
                     HaachtSchor=Plots_HaachtSchor,Diepenbeek=Plots_Diepenbeek)
  
  Calc_Comm = function(Data,traits){
    Community_values(Data,traits)
  }
  
  Plots_Comm <- lapply(Plots_list, Calc_Comm,traits=traits)
  return(Plots_Comm)
}
Plots_CWM_FD_save <- function (wd2,Plots_Comm,approach,save) {
  # Paste all info of different sites in 1 data frame for CWM and 1 data frame for FD metrics
  
  listOfDataFrames_CWM <- lapply(Plots_Comm, function(x) x[["CWM"]])
  df_CWM<-do.call(rbind.data.frame, listOfDataFrames_CWM)  # least time-efficient but controls for matching column names

  # get rownames from all data.frames in a list
  Site_Plot <- ldply(lapply(listOfDataFrames_CWM, rownames), data.frame)
  colnames(Site_Plot) <- c("Site", "Plot")
 
  setwd(wd2)
  
  if (save == "CWM"){
    df_CWM <- cbind(Site_Plot,df_CWM)
    write.table(df_CWM,file=paste("CWM_",".csv",sep=approach),sep=";",row.names = F)
  } else if (save == "FD"){
    FDlist <- c("FRic","FEve","FDiv","FDis","RaoQ")
    FDmatrix <- matrix(data=NA, ncol=5,nrow=nrow(Site_Plot))
    colnames(FDmatrix) <- FDlist
    for (i in 1: length(FDlist)){
      listOfDataFrames_FD <- lapply(Plots_Comm, function(x) x[[FDlist[i]]])
      df_FD <- list_vect2df(listOfDataFrames_FD)
      colnames(df_FD) <- c("Site","Plot",FDlist[i])
      df_FD <- df_FD[order(match(df_FD[,"Plot"],Site_Plot[,"Plot"])),]
      df_FD_values <- df_FD[FDlist[i]]
      FDmatrix[,FDlist[i]]<-df_FD_values[,1]
    }
    df_FD <- cbind(Site_Plot,FDmatrix)
    write.table(df_FD,file=paste("FD_",".csv",sep=approach),sep=";",row.names = F)
  }
}

#--------------------------------------------------------------------------------
# APPLY
# -------------- 1) Conventionally measured functional traits
traits <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")
traits_FD1 <- c("SLA","LDMC","LNC_mass","Height")
traits_FD2 <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log")
wd1 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data"
wd2 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD"
wd3 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data"

### CWMs 
PFT_CWM <- Plots_CWM_FD(wd1, wd2, wd3, approach="PFTraits", traitlist=traits, type="")
Plots_CWM_FD_save(wd2, PFT_CWM, "PFTraits", save="CWM")
# 2 plots were deleted for the calculation of CWM_NAT (S4.1 and S8.2)
PFT_CWM_NAT <- Plots_CWM_FD(wd1, wd2, wd3, approach="PFTraits_NAT", traitlist=traits, type="")
Plots_CWM_FD_save(wd2, PFT_CWM_NAT, "PFTraits_NAT", save="CWM")

### FDs
PFT_FD1 <- Plots_CWM_FD(wd1, wd2, wd3, approach="PFTraits", traitlist=traits_FD1, type="")
PFT_FD2 <- Plots_CWM_FD(wd1, wd2, wd3, approach="PFTraits", traitlist=traits_FD2, type="")
Plots_CWM_FD_save(wd2, PFT_FD1, "PFTraits_1", save="FD")
Plots_CWM_FD_save(wd2, PFT_FD2, "PFTraits_2", save="FD")

# -------------- 1) Optically measured functional traits
traits <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")
traits_FD1 <- c("SLA","LDMC","LNC_mass","Height")
traits_FD2 <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log")
wd1 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/c) PLSR best models overview/"
wd2 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD"
wd3 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data"

### CWMs
POT_CWM <- Plots_CWM_FD(wd1, wd2, wd3, approach="POTraits", traitlist=traits, type="cal")
Plots_CWM_FD_save(wd2, POT_CWM, "POTraits", save="CWM")
POT_CWM_NAT <- Plots_CWM_FD(wd1, wd2, wd3, approach="POTraits_NAT", traitlist=traits, type="cal")
Plots_CWM_FD_save(wd2, POT_CWM_NAT, "POTraits_NAT", save="CWM")

### FDs
POT_FD1 <- Plots_CWM_FD(wd1, wd2, wd3, approach="POTraits", traitlist=traits_FD1, type="cal")
Plots_CWM_FD_save(wd2, POT_FD1, "POTraits_1", save="FD")
POT_FD2 <- Plots_CWM_FD(wd1, wd2, wd3, approach="POTraits", traitlist=traits_FD2, type="cal")
Plots_CWM_FD_save(wd2, POT_FD2, "POTraits_2", save="FD")