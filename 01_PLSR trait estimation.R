######################################################################
###                                                                ###
###     Multivariate trait estimation from hyperspectral data      ###
###                                                                ###
######################################################################

### The objective of this script is to estimate traits from hyperspectral data.
# We will use PLSR and assess how wel the estimates reflect conventionally measured traits.

### The input of these steps are the average spectra of the bouquets/patches of 
# - the unmixed spectra
# - with or without brightness normalization (_"br"): this is performed in "Load_PFT_spec"

# Script by Elisa Van Cleemput, 2017-18
#--------------------------------------------------------------------------------
# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------
# Load functional and spectral data
Load_PFT_spec <- function(proc,traitset,LUTdata, width) {
  
  #--------------------------------------------------------------------------------
  # 1) Conventionally measured functional traits
  
  PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
  PFT <- PFT[order(PFT$ID),] # order the observations according to ID
  target <- which(names(PFT) == 'ID')[1]
  PFT_meta<-PFT[,1:target]
  PFT_traits<-PFT[,-c(1:target)]

  PFT_traits$LNC_mass <- 10 * PFT_traits$LNC # % dry mass to mg/g
  PFT_traits$LCaC_mass <- 10 * PFT_traits$Ca # % dry mass to mg/g
  PFT_traits$LMgC_mass <- 10 * PFT_traits$Mg # % dry mass to mg/g

  PFT_traits$LNC_mass_log <- log10(PFT_traits$LNC_mass)
  PFT_traits$SLA_log <- log10(PFT_traits$SLA)
  PFT_traits$LDMC_log <- log10(PFT_traits$LDMC)
  PFT_traits$Height_log <- log10(PFT_traits$Height)
  PFT_traits$LCaC_mass_log <- log10(PFT_traits$LCaC_mass)
  PFT_traits$LMgC_mass_log <- log10(PFT_traits$LMgC_mass)

  PFT_traits$LNC_mass_sqrt <- sqrt(PFT_traits$LNC_mass)
  PFT_traits$SLA_sqrt <- sqrt(PFT_traits$SLA)
  PFT_traits$LDMC_sqrt <- sqrt(PFT_traits$LDMC)
  PFT_traits$Height_sqrt <- sqrt(PFT_traits$Height)
  PFT_traits$LCaC_mass_sqrt <- sqrt(PFT_traits$LCaC_mass)
  PFT_traits$LMgC_mass_sqrt <- sqrt(PFT_traits$LMgC_mass)
  
  #--------------------------------------------------------------------------------
  # 2) Reflectance
  
  spec <- read.csv(paste("Species_x_reflectance_",sep=proc,".csv"),sep=";",header=TRUE)
  target <- which(names(spec) == 'nm350')[1]
  spectra<-spec[,target:ncol(spec)]
  meta<-spec[,1:target-1]
  rownames(meta)=meta[,"ID"]
  
  # Preprocessing steps: spectral data
  # 1) smoothing
  # 2) Create library for 1st deriv or 2nd deriv spectra
  # 3a) spectral binning (optional) + removal of noise bands
  # 3b) removal of noise bands
  # 4) brightness normalization (Feilhauer et al. 2010)
  
  
  ### 1) Create speclib and Smooth
  library(hsdar)
  wl=as.numeric(gsub("nm","",names(spectra)))
  speclib <- speclib(as.matrix(spectra),wl)
  SI(speclib)<-meta
  idSpeclib(speclib)<-as.character(meta$ID)
  
  library(stringr) #str_sub
  SI(speclib)$Instrument <- str_sub(SI(speclib)[,1],-3,-1)
  libr_sv<-subset(speclib,Instrument == ".sv")
  libr_asd<-subset(speclib,Instrument == "asd")
  libr_sv<- smoothSpeclib(libr_sv,method="sgolay", n=101)
  libr_asd<- smoothSpeclib(libr_asd,method="sgolay", n=51)
  spectra_sv <- spectra(libr_sv)
  meta_sv <- SI(libr_sv)
  spectra_asd <- spectra(libr_asd)
  meta_asd <- SI(libr_asd)
  spectra_total <- rbind(spectra_asd,spectra_sv)
  wl <- seq(350,2500,1)
  colnames(spectra_total)<-wl
  meta_total <- rbind(meta_asd,meta_sv)
  Refl <- cbind(meta_total, spectra_total)
  Refl <- Refl[order(Refl$ID),] # order the observations according to ID
  meta_ordered <- Refl[,c(1:which(names(Refl)=="Instrument"))]
  speclib_smooth <- speclib(as.matrix(Refl[,-c(1:which(names(Refl)=="Instrument"))]), wl)
  SI(speclib_smooth)<-meta_ordered
  
  # speclib_smooth<- smoothSpeclib(speclib,method="sgolay", n=51)
  
  ### 2) Create library for 1st deriv or 2nd deriv spectra
  if (isTRUE(grepl("1deriv",LUTdata,fixed=T))){
    speclib_smooth <- derivative.speclib(speclib_smooth, m = 1)
  } else if (isTRUE(grepl("2deriv",LUTdata,fixed=T))){
    speclib_smooth <- derivative.speclib(speclib_smooth, m = 2)
  } else {
    speclib_smooth <- speclib_smooth
  }
  
  ### 3a) Spectral binning (optional) + removal of noise bands
  if (width == "10"){
    # function for spectral binning
    wl <- c(seq(400, 1340, 10), seq(1460,1780,10), seq (1970, 2400, 10))
    resamp <- function (spec, wl, wl.out){
      out <- approx (wl, spec, xout=wl.out)$y
    }
    
    # apply function to spectra --> result = smoothed spectra in matrix: binned + chosen spectral regions only
    spectra_smooth <- speclib_smooth@spectra@spectra_ma # smoothed spectra in matrix
    wavelength <- speclib_smooth@wavelength # wavelengths in vector
    spectra <- t (apply (spectra_smooth, 1, resamp, wl=wavelength, wl.out=wl)) 
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation
    speclib_plot<-speclib(spectra,wl) #when visualising binned spectra
    SI(speclib_plot)<-meta_ordered
    colnames(spectra) <- wl
  }  else if (width == "1"){
    ### 3b) Removal of noise bands in case of no spectral binning
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation
    spectra <- speclib_smooth@spectra@spectra_ma # smoothed spectra in matrix
    wl<-speclib_smooth@wavelength
    colnames(spectra) <- wl
  }
  
  #### 4) Brightness normalization (optional) (only for raw and unmixed spectra, not 1deriv and 2deriv)
  if (isTRUE(grepl("br",LUTdata,fixed=T))){
    # Binned or detailed spectra
    brightness<-sqrt(apply(spectra^2,1,sum))
    spectra<-spectra/brightness
    speclib_smooth<-speclib(spectra,wl)
    SI(speclib_smooth)=meta_ordered
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation, this is binned in case of binning!
    colnames(spectra) <- wl
  } 
  
  output <- list ("PFT" = PFT,
                  "PFT_traits" = PFT_traits,
                  "PFT_meta" = PFT_meta,
                  "spectra" = spectra,
                  "spec_meta" = meta_ordered,
                  "wl" = wl)
  return(output)
}

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data")
data_PFT_spec <- Load_PFT_spec(proc="unmixed", LUTdata = "_refl_sg", width = 1)
data_PFT_spec_br <- Load_PFT_spec(proc="unmixed", LUTdata = "_refl_sg_br", width = 1)

#--------------------------------------------------------------------------------
# Exploration: identify outliers in the dataset
#----------------------------
# http://r-statistics.co/Outlier-Treatment-With-R.html:
# For a given continuous variable, outliers are those observations that lie outside 1.5 * IQR,
# where IQR, the 'Inter Quartile Range' is the difference between 75th and 25th quartiles. 

traits <- c("LNC_mass","LNC_mass_log","LNC_mass_sqrt", 
            "SLA","SLA_log", "SLA_sqrt", 
            "LDMC", "LDMC_log", "LDMC_sqrt",
            "Height","Height_log","Height_sqrt",
            "LCaC_mass","LCaC_mass_log","LCaC_mass_sqrt",
            "LMgC_mass","LMgC_mass_log","LMgC_mass_sqrt")

my_outliers <- function(data, traits, name){
  PFT <- data$PFT_traits
  meta <- data$PFT_meta
  
  trait_outliers <- matrix(nrow=10, ncol=length(traits))
  colnames(trait_outliers) <- traits
  
  tiff(paste("Boxplot_traits_",".tiff",sep=name),res=300,width=15,height=8,units='in')
  # x11(width=1280,height=720)
  fig <- par(mfrow=c(3,4))
  
  for (i in 1:length(traits)){
    outlier_values <- boxplot.stats(PFT[,traits[i]])$out
    boxplot(PFT[,traits[i]], main=traits[i], boxwex=0.1)

    outlier_idx <- rownames(PFT[which(PFT[,traits[i]] %in% outlier_values),])
    outlier_ID <- meta[outlier_idx,"ID"]
    
    trait_outliers[,i]<-c(paste(as.character(outlier_ID), round(outlier_values,2), sep=": "),
                          rep(NA,10-length(outlier_ID)))
  }
  dev.off()
  
  write.table(trait_outliers, file=paste("Outliers_traits_",".csv",sep=name), sep =";", row.names=F, col.names=T) 
  
  return(trait_outliers)
}

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results")
outliers <- my_outliers(data_PFT_spec, traits, name="convent")

#--------------------------------------------------------------------------------
# PLSR phase 1: determine optimal number of components
#----------------------------
# The trait_ncomp_ini function fits a PLSR for each trait (response variable) given to the function, and
# determines and stores the ideal number of latent variables (LVs), according to three approaches.

# Input: - data = list containing "PFT_traits", "PFT_meta", "spectra", "spec_meta" and "wl"
#        - traits =  vector with trait names
#        - ncomp_ini = set a random initial number of no. of LVs to limit the extent of the plot (the same for each trait)
#        - ncomp_treshold = a series of numbers that allows you to visualise a self-chosen no. of LVs per trait
#        - name = string to store the results with a unique name

# "Cross-validation is commonly used to determine the optimal number of components to take into account"
# Based on my_plsr_ini results we define no. of components to be retained (subjectively, based on RMSEP and R²)
#   Mevik & Wehrens (2018); Salter (2018)

# https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf

# By default; pls scales the data

library(pls)
my_plsr_ini <- function(data, traits, ncomp_ini, ncomp_treshold, name){
  PFT <- data$PFT_traits
  spectra <- data$spectra
  wl <- data$wl
  
  # Store n components
  ncomp_onesigma <- vector(mode="numeric",length(traits))
  ncomp_permut <- vector(mode="numeric",length(traits))
  ncomp_press <- vector(mode="numeric",length(traits))
    
  for (i in 1 : length(traits)){
    PLSR.trait<-plsr(PFT[,traits[i]]~spectra, ncomp = ncomp_ini, validation="LOO")

    tiff(paste(paste(name,"_",sep= traits[i]),"comp.tiff",sep= toString(ncomp_ini)),res=300,width=15,height=8,units='in')
    # x11(width=1280,height=720)
    fig <- par(mfrow=c(2,4))
    plot(RMSEP(PLSR.trait), legendpos = "topright", main=paste(traits[i], "(LOO-CV)", sep=" ")) # same as: plot(PLSR.trait, plottype = "validation") OR validationplot(PLSR.trait)
    plot(x=seq(0,ncomp_ini),y=RMSEP(PLSR.trait)$val[1,,]/mean(PFT[,traits[i]]), main=paste(traits[i], "(LOO-CV)", sep=" "), ylab = "% RMSEP") 
    points(x=seq(0,ncomp_ini),y=RMSEP(PLSR.trait)$val[2,,]/mean(PFT[,traits[i]]), main=paste(traits[i], "(LOO-CV)", sep=" "), col = "red", new=F) 
    plot(R2(PLSR.trait, estimate="CV"), legendpos = "bottomright", main=paste(traits[i], "(LOO-CV)", sep=" ")) # R2 = Q² = 1 - PRESS/TSS
    
    ncomp.onesigma <- selectNcomp(PLSR.trait, method = "onesigma", plot = TRUE, main = "one-sigma heuristic") 
    ncomp.press <- which.min(PLSR.trait$validation$PRESS)
    abline(v=ncomp.press, col="red",lty=2)
    abline(v=ncomp_treshold[i], col="green",lty=2)
    
    plot(RMSEP(PLSR.trait, estimate="train"), legendpos = "topright", main=paste(traits[i], "(train)", sep=" "))
    plot(x=seq(0,ncomp_ini),y=RMSEP(PLSR.trait, estimate="train")$val[1,,]/mean(PFT[,traits[i]]), main=paste(traits[i], "(train)", sep=" "), ylab = "% RMSEP") 
    plot(R2(PLSR.trait, estimate="train"), legendpos = "bottomright", main=paste(traits[i], "(train)", sep=" ")) # R2 = R² = 1 - RSS/TSS

    ncomp.permut<- selectNcomp(PLSR.trait, method = "randomization", plot = TRUE, main = "permutation approach")
    par(fig)
    dev.off()
    # https://stats.stackexchange.com/questions/12900/when-is-r-squared-negative
    
    ncomp_onesigma[i] <- ncomp.onesigma
    ncomp_permut[i] <- ncomp.permut
    ncomp_press[i] <- ncomp.press
  }
  
  output <- list (ncomp_onesigma = ncomp_onesigma,
                  ncomp_permut = ncomp_permut,
                  ncomp_press = ncomp_press)
}

# Results (also part of the input, for visualisation):
trait_ncomp_ini <- rep(5,length(traits))

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/a) PLSR find ncomp/")
ncomp_spec <- my_plsr_ini(data_PFT_spec, traits, 20, trait_ncomp_ini, name = "PLSR_spec_")
ncomp_spec_br <- my_plsr_ini(data_PFT_spec_br, traits, 20, trait_ncomp_ini, name = "PLSR_spec_br_")

#--------------------------------------------------------------------------------
# PLSR phase 2: Fit models and save predicted trait values
#----------------------------
# The my_plsr function fits a PLSR for each trait (response variable) given to the function, with the no. of latent variables fixed
# Models are fitted with LOO-CV

# Output: - Figure per trait showing scatterplot, statistics, regression coefficients and loadings
#         - csv with model statistics
#         - csv with model predictions of calibration
#         - csv with model predictions of LOO-CV

# Input: - data = list containing "PFT_traits", "PFT_meta", "spectra", "spec_meta" and "wl"
#        - traits =  vector with trait names
#        - ncomp_traits = no. of latent variables for each trait
#        - name = string to store the results with a unique name
my_plsr <- function(data, traits, ncomp_traits, name)  {
  PFT <- data$PFT_traits
  spectra <- data$spectra
  wl <- data$wl
  
  traits_predicted_cal <- matrix(nrow = nrow(PFT), ncol = length(traits))
  colnames(traits_predicted_cal) <- traits
  rownames(traits_predicted_cal) <- data$spec_meta$ID
  traits_predicted_cv <- traits_predicted_cal
  
  traits_pred_acc <- matrix(nrow = 13, ncol = length(traits))
  rownames(traits_pred_acc) <- c("RMSEP_CV","nRMSEP_CV","RMSEP_CV_adj","nRMSEP_CV_adj","RMSE_train","nRMSE_train",
                                 "R2_CV","R2_train",
                                 "TheilU_CV_1","TheilU_CV_2", "TheilU_train_1", "TheilU_train_2",
                                 "ncomp")
  colnames(traits_pred_acc) <- traits
  trait_coeff <- matrix (nrow = length(traits), ncol =dim(spectra)[2])
  colnames(trait_coeff) <- wl
  rownames(trait_coeff) <- traits
  
  for (i in 1 : length(traits)){
    if (ncomp_traits[i] == 0){
      print(paste("bad plsr model for", traits[i]))
    } else {
      PLSR.trait<-plsr(PFT[,traits[i]]~spectra, ncomp = ncomp_traits[i], validation="LOO")
      
      tiff(paste(name,".tiff",sep= traits[i]),res=300,width=15,height=8,units='in')
      # x11(width=1280,height=720)
      fig <- par(mfrow=c(2,3))

      # 1. Cross-validated predictions vs. measured values
      RMSEP_CV      <- round(RMSEP(PLSR.trait)$val[1,,ncomp_traits[i]+1],2)
      RMSEP_CV_adj      <- round(RMSEP(PLSR.trait)$val[2,,ncomp_traits[i]+1],2)
      # nRMSEP_CV     <- round(RMSEP(PLSR.trait)$val[1,,ncomp_traits[i]+1]/PLSR.trait$Ymean,2)
      # nRMSEP_CV     <- round(RMSEP(PLSR.trait)$val[1,,ncomp_traits[i]+1]/(max(PLSR.trait$fitted.values)-min(PLSR.trait$fitted.values)),2)
      nRMSEP_CV     <- round(RMSEP(PLSR.trait)$val[1,,ncomp_traits[i]+1]/(max(PFT[,traits[i]])-min(PFT[,traits[i]])),2)
      nRMSEP_CV_adj     <- round(RMSEP(PLSR.trait)$val[2,,ncomp_traits[i]+1]/(max(PFT[,traits[i]])-min(PFT[,traits[i]])),2)
      RMSE_train   <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp_traits[i]+1],2)
      # nRMSE_train  <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp_traits[i]+1]/PLSR.trait$Ymean,2)
      # nRMSE_train  <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp_traits[i]+1]/(max(PLSR.trait$fitted.values)-min(PLSR.trait$fitted.values)),2)
      nRMSE_train  <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp_traits[i]+1]/(max(PFT[,traits[i]])-min(PFT[,traits[i]])),2)
      R2_CV         <- round(R2(PLSR.trait)$val[1,,ncomp_traits[i]+1],2)
      R2_train      <- round(R2(PLSR.trait, estimate="train")$val[1,,ncomp_traits[i]+1],2)
      traits_pred_acc["RMSEP_CV",i] <- RMSEP_CV
      traits_pred_acc["nRMSEP_CV",i] <- nRMSEP_CV
      traits_pred_acc["RMSEP_CV_adj",i] <- RMSEP_CV_adj
      traits_pred_acc["nRMSEP_CV_adj",i] <- nRMSEP_CV_adj
      traits_pred_acc["RMSE_train",i] <- RMSE_train
      traits_pred_acc["nRMSE_train",i] <- nRMSE_train
      traits_pred_acc["R2_CV",i] <- R2_CV
      traits_pred_acc["R2_train",i] <- R2_train
      traits_pred_acc["ncomp",i] <- ncomp_traits[i]
      
      par(pty="s")
      plot(PLSR.trait, ncomp = ncomp_traits[i], asp = 1, line = F, main=traits[i]) # same as plottype = "prediction" OR predplot(PLSR.trait, ncomp = ncomp_traits[i], which="validation")
      # legend("bottomright",bty="n",legend=c(RMSEP_CV, nRMSEP_CV, RMSE_train, nRMSE_train, R2_CV, R2_train))
      abline(0,1, lty=2)
      
      par(pty="s")
      plot(PLSR.trait, ncomp = ncomp_traits[i], asp = 1, line = F, main=traits[i], col="white")
      legend("center",bty="n",legend=c(paste("RMSEP_CV =",RMSEP_CV),  paste("nRMSEP_CV =",nRMSEP_CV),
                                       paste("RMSE_train =",RMSE_train), paste("nRMSE_train =", nRMSE_train),
                                       paste("R2_CV =",R2_CV), paste("R2_train =", R2_train)))
      
      # 2. Correlation plot (loadings)
      # fraction of the variance of the variables explained by the components (circles = 50% and 100%)
      if(ncomp_traits[i] == 1){
        plot(PLSR.trait, ncomp = ncomp_traits[i], asp = 1, line = F, main=traits[i], col="white")
        legend("center", bty="n", legend=c("Only one component, so impossible to create correlation loadings plot"))
      } else {
        par(pty="m")
        plot(PLSR.trait, plottype = "correlation", main="Correlation loadings plot") 
      }
      
      # 3. Regression coefficients
      # Difference between loadings and coefficients:
      #   - interpreting individual components: use loadings (or loading.weights) = regression coefficients of each latent variable
      #   - interpreting the whole model: use regression coefficients or VIP scores = combination of regression coefficients of latent variables and regression coefficients of overall model
      # Trait = a*(LV1)+b*(LV2)+c*(LV3)+d*(LV4)
      # LV1 = e*(wl1)+f*(wl2)+g*(wl3)+h*(wl4)+...
      
      # 3a. Combination of all latent variables
      plot(PLSR.trait, plottype = "coef",ylim=c(min(PLSR.trait$coefficients),max(PLSR.trait$coefficients)),
           labels = "numbers", xlab = "nm", main="Regression coeffcients")
      xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
      xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
      yy_mask1<-c(seq(min(PLSR.trait$coefficients),min(PLSR.trait$coefficients),length.out=length(xx_mask1)/2),
                  seq(max(PLSR.trait$coefficients),max(PLSR.trait$coefficients),length.out=length(xx_mask1)/2))
      yy_mask2<-c(seq(min(PLSR.trait$coefficients),min(PLSR.trait$coefficients),length.out=length(xx_mask2)/2),
                  seq(max(PLSR.trait$coefficients),max(PLSR.trait$coefficients),length.out=length(xx_mask2)/2))
      polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
      polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
      abline(h = 0, col=c("darkgray"))
      
      # 3b. Each latent variable separately
      maxvalues <- array(0,ncomp_traits[i])
      minvalues <- array(0,ncomp_traits[i])
      for (j in 1:ncomp_traits[i]){
        # in the previous graph, max and min were extracted from: coef0 <- coef(PLSR.trait, ncomp = 5, comps = NULL, intercept = F)
        # By default: intercept = F
        coef_j<- coef(PLSR.trait, ncomp = ncomp_traits[i], comps = j, intercept = F)
        maxvalues[j] <- max(coef_j)
        minvalues[j] <- min(coef_j)
      }
      plot(PLSR.trait, plottype = "coef", comps = 1:ncomp_traits[i],legendpos = "bottomright",
           ylim=c(min(minvalues),max(maxvalues)),labels = "numbers", xlab = "nm", main="Regression coeffcients")
      xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
      xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
      yy_mask1<-c(seq(min(minvalues),min(minvalues),length.out=length(xx_mask1)/2),
                  seq(max(maxvalues),max(maxvalues),length.out=length(xx_mask1)/2))
      yy_mask2<-c(seq(min(minvalues),min(minvalues),length.out=length(xx_mask2)/2),
                  seq(max(maxvalues),max(maxvalues),length.out=length(xx_mask2)/2))
      polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
      polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
      
      
      # 4. Loadings
      plot(PLSR.trait, "loadings", comps = 1:ncomp_traits[i], legendpos = "bottomright",labels = "numbers", xlab = "nm", main=paste("Loadings", traits[i],sep=" "))
      xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
      xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
      yy_mask1<-c(seq(min(PLSR.trait$loadings),min(PLSR.trait$loadings),length.out=length(xx_mask1)/2),
                  seq(max(PLSR.trait$loadings),max(PLSR.trait$loadings),length.out=length(xx_mask1)/2))
      yy_mask2<-c(seq(min(PLSR.trait$loadings),min(PLSR.trait$loadings),length.out=length(xx_mask2)/2),
                  seq(max(PLSR.trait$loadings),max(PLSR.trait$loadings),length.out=length(xx_mask2)/2))
      polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
      polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
      abline(h = 0, col=c("darkgray"))
      explvar(PLSR.trait)
      
      # 5. Score plot
      # plot(PLSR.trait, plottype = "scores", comps = 1:ncomp_traits[i], main="Score plot")
      par(fig)
      dev.off()
      
      # Save predicted values
      traits_predicted_cal[,i] <- predict(PLSR.trait, ncomp = ncomp_traits[i]) # training predictions # idem to:
      traits_predicted_cv[,i] <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp_traits[i], which="validation"))[,"predicted"] # LOO-CV predictions 
      
      # Store regression coefficients
      trait_coeff[i,] <- coef(PLSR.trait)
      
      # Extra's from Schweiger et al. 2017
      # ---- Theil's U
      obs_pred_cv <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp_traits[i], which="validation"))
      obs_pred_cal <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp_traits[i], which="train"))

      library(DescTools)
      traits_pred_acc["TheilU_CV_1",i] <- round(TheilU(obs_pred_cv$measured, obs_pred_cv$predicted, type = 1),2)
      traits_pred_acc["TheilU_CV_2",i] <- round(TheilU(obs_pred_cv$measured, obs_pred_cv$predicted, type = 2),2)
      traits_pred_acc["TheilU_train_1",i] <- round(TheilU(obs_pred_cal$measured, obs_pred_cal$predicted, type = 1),2)
      traits_pred_acc["TheilU_train_2",i] <- round(TheilU(obs_pred_cal$measured, obs_pred_cal$predicted, type = 2),2)
      
    }
  }
  write.table(traits_predicted_cal, file=paste(name,".csv",sep="predictions_cal"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_predicted_cv, file=paste(name,".csv",sep="predictions_Loo-cv"), sep =";", row.names=T, col.names=NA) 
  
  write.table(trait_coeff, file=paste(name,".csv",sep="regr_coeff"), sep =";", row.names=T, col.names=NA) 
  
  traits_pred_acc <- t(traits_pred_acc)
  write.table(traits_pred_acc, file=paste(name,".csv",sep="Accuracy"), sep =";", row.names=T, col.names=NA) 
  
  return(traits_pred_acc)
}

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/b) PLSR models/")
plsr_acc_press_PFT_spec <- my_plsr(data_PFT_spec, traits, ncomp_spec$ncomp_press, name = "PLSR_spec_")
plsr_acc_press_PFT_spec_br <- my_plsr(data_PFT_spec_br, traits,  ncomp_spec_br$ncomp_press, name = "PLSR_spec_br_")

#--------------------------------------------------------------------------------
# PLSR phase 3: Identify best model per trait
#----------------------------
# For each trait, the model_best_results function compares the model statistics obtained by different plsr modelling approaches
# The output is an overview of the best model per trait, based on the different model statistics

# Input: - traits =  vector with trait names
#        - data = list containing different plsr approaches for each trait
#        - name = string to store the results with a unique name

model_best_results <- function(traits, data, name){
    traits_pred_acc_best <- matrix(nrow = length(traits)*2, ncol = 12)
    traits_pred_acc_best <- as.data.frame(traits_pred_acc_best)
    colnames(traits_pred_acc_best) <- c("RMSEP_CV","nRMSEP_CV","RMSEP_CV_adj","nRMSEP_CV_adj","RMSE_train","nRMSE_train","R2_CV","R2_train","TheilU_CV_1","TheilU_CV_2", "TheilU_train_1", "TheilU_train_2")
  
  for (i in 1:length(traits)){
    for (j in 1:6){
      temp <- sapply(data, function(x) x[i,j])
      traits_pred_acc_best[i*2-1,j] <- min(temp)
      rownames(traits_pred_acc_best)[i*2-1] <- traits[i]
      rownames(traits_pred_acc_best)[i*2] <- paste("techn",traits[i],sep=".")
      if(length(labels(temp)[which.min(temp)]) == 0){
        traits_pred_acc_best[i*2,j] <- NA
      } else {
        traits_pred_acc_best[i*2,j] <- labels(temp)[which.min(temp)]
      }}
    for (j in 7:8){
      temp <- sapply(data, function(x) x[i,j])
      traits_pred_acc_best[i*2-1,j] <- max(temp)
      if(length(labels(temp)[which.max(temp)]) == 0){
        traits_pred_acc_best[i*2,j] <- NA
      } else {
        traits_pred_acc_best[i*2,j] <- labels(temp)[which.max(temp)]
      }}
    for (j in 9:ncol(traits_pred_acc_best)){
      temp <- sapply(data, function(x) x[i,j])
      traits_pred_acc_best[i*2-1,j] <- min(temp)
      rownames(traits_pred_acc_best)[i*2-1] <- traits[i]
      rownames(traits_pred_acc_best)[i*2] <- paste("techn",traits[i],sep=".")
      if(length(labels(temp)[which.min(temp)]) == 0){
        traits_pred_acc_best[i*2,j] <- NA
      } else {
        traits_pred_acc_best[i*2,j] <- labels(temp)[which.min(temp)]
      }}
    }
  
  # Save overview
  write.table(traits_pred_acc_best, file=paste(name,".csv",sep="best_models"), sep =";", row.names=T, col.names=NA) 
  
  return(traits_pred_acc_best)
}

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/c) PLSR best models overview/")

# For each trait we fitted a traditional and a brightness-normalized plsr
models <- list("Unmixed" = plsr_acc_press_PFT_spec, 
               "Unmixed_br" = plsr_acc_press_PFT_spec_br)
plsr_best_press_PFT_spec <- model_best_results(traits, models, name="PLSR_")

#--------------------------------------------------------------------------------
# PLSR phase 4: Store best prediction per trait
#----------------------------
# The function extract_best_results consults the previously stored model predictions and
# concatenates the results fom the best models in 1 table

# Within this function, the results generated through the my_plsr function are loaded --> manually defined)

# Input
# Inspect the table generated by the model_best_results function,
# and indicate for each trait which transformation of the trait and which plsr approach 
# you want to maintain for further ecological analyses
best_trait_transf <- c("SLA","LDMC","LNC_mass","LCaC_mass_log","LMgC_mass_log","Height")
best_spec_transf <- c("unmixed_br","unmixed","unmixed_br","unmixed_br","unmixed_br","unmixed")

extract_best_results <- function (traits, best_data){
  Pred_cal_unmixed <- read.csv("PLSR_spec_predictions_cal.csv",sep=";",header=TRUE, row.names=1)
  Pred_cal_unmixed_br <- read.csv("PLSR_spec_br_predictions_cal.csv",sep=";",header=TRUE, row.names=1)
  Pred_cv_unmixed <- read.csv("PLSR_spec_predictions_Loo-cv.csv",sep=";",header=TRUE, row.names=1)
  Pred_cv_unmixed_br <- read.csv("PLSR_spec_br_predictions_Loo-cv.csv",sep=";",header=TRUE, row.names=1)
  
  Pred_unmixed_stats <- read.csv("PLSR_spec_Accuracy.csv",sep=";",header=TRUE, row.names=1)
  Pred_unmixed_br_stats <- read.csv("PLSR_spec_br_Accuracy.csv",sep=";",header=TRUE, row.names=1)
  
  Pred_unmixed_regr_coeff <- read.csv("PLSR_spec_regr_coeff.csv",sep=";",header=TRUE, row.names=1)
  Pred_unmixed_br_regr_coeff <- read.csv("PLSR_spec_br_regr_coeff.csv",sep=";",header=TRUE, row.names=1)
  
  # Create empty matrices to fill
  best_predictions_cal <- matrix(nrow=nrow(Pred_cal_unmixed), ncol=length(traits))
  rownames(best_predictions_cal) <- rownames(Pred_cal_unmixed)
  colnames(best_predictions_cal) <- traits
  best_predictions_cv <- best_predictions_cal
  best_predictions_stats <- as.data.frame(matrix(nrow=length(traits), ncol = ncol(Pred_unmixed_stats)))
  rownames(best_predictions_stats) <- traits
  colnames(best_predictions_stats) <- colnames(Pred_unmixed_stats)
  best_predictions_regr_coeff <- data.frame()
  
  for (i in 1:length(traits)){
    if (best_data[i] == "unmixed"){
      best_predictions_cal[,i] <- Pred_cal_unmixed[,c(traits[i])]
      best_predictions_cv[,i] <- Pred_cv_unmixed[,c(traits[i])]
      best_predictions_stats[i,] <- Pred_unmixed_stats[c(traits[i]),]
      best_predictions_regr_coeff <- rbind(best_predictions_regr_coeff,Pred_unmixed_regr_coeff[c(traits[i]),])
      
    } else if (best_data[i] == "unmixed_br"){
      best_predictions_cal[,i] <- Pred_cal_unmixed_br[,c(traits[i])]
      best_predictions_cv[,i] <- Pred_cv_unmixed_br[,c(traits[i])]
      best_predictions_stats[i,] <- Pred_unmixed_br_stats[c(traits[i]),]
      best_predictions_regr_coeff <- rbind(best_predictions_regr_coeff,Pred_unmixed_br_regr_coeff[c(traits[i]),][c(traits[i]),])
    }
  }
  best_predictions_stats <- cbind(best_data, best_predictions_stats)
  best_predictions_stats.print <- best_predictions_stats[,c("R2_train","R2_CV", "nRMSE_train", "nRMSEP_CV", "TheilU_train_1", "TheilU_CV_1", "ncomp")]

  output <- list("best_predictions_cal" = best_predictions_cal,
                 "best_predictions_cv" = best_predictions_cv,
                 "best_predictions_stats.print" = best_predictions_stats.print,
                 "best_predictions_regr_coeff" = best_predictions_regr_coeff)
  
  return(output)
}

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/b) PLSR models/")
best_results <- extract_best_results(best_trait_transf, best_spec_transf)

# Stores the output
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/c) PLSR best models overview/")
write.table(best_results$best_predictions_cal, file="PLSR_best_predictions_cal.csv", sep =";", row.names=T, col.names=NA) 
write.table(best_results$best_predictions_cv, file="PLSR_best_predictions_cv.csv", sep =";", row.names=T, col.names=NA) 
write.table(best_results$best_predictions_stats.print, file="PLSR_best_predictions_accuracy.csv", sep =";", row.names=T, col.names=NA) 

#--------------------------------------------------------------------------------
# Visualise regression coefficients
#----------------------------
library(reshape2)
library(ggplot2)
wl <- as.numeric(sub("X","",names(best_results$best_predictions_regr_coeff),fixed = TRUE))
regr_coeff <- data.frame(cbind(wl, t(best_results$best_predictions_regr_coeff)))
colnames(regr_coeff) <-c("wl","SLA","LDMC","Leaf N", "Leaf Ca (log)", "Leaf Mg (log)", "Height")

stand <- function(x) {
  stand_factor <- max(abs(x))
  return(x/stand_factor)
}
regr_coeff_stand <- as.data.frame(apply(regr_coeff, 2, stand))
regr_coeff_stand$wl <- wl

stats.melted <- melt(regr_coeff_stand, id.vars = "wl")
colnames(stats.melted)[3] <- "coeff"
stats.melted$refl_br <- apply(data_PFT_spec_br$spectra, 2, mean)
stats.melted$refl <- apply(data_PFT_spec$spectra, 2, mean)

# Visualise regression coefficients
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/PLSR regression coefficients/")
# x11(width=1000,height=720)
ggplot(stats.melted) + 
  geom_point(aes(x = wl, y = coeff, group=1), size=0.5, color = "dodgerblue4") +
  geom_point(aes(x = wl, y = refl*2, group=1), size=0.5, color = "gray50") +
  scale_y_continuous(sec.axis = sec_axis(~./2*100, name = "Reflectance [%]", breaks = c(0, 25, 50))) +
  geom_point(aes(x = wl, y = coeff, group=1), size=0.5, color = "dodgerblue4") +
  ylab("standardized regression coefficients") +
  xlab("Wavelength (nm)") +
  facet_wrap(~ variable) +
  theme_bw()+
  theme(plot.margin = margin(10, 15, 10, 10))
ggsave("Regr_coefficients.jpeg",width = 25,height = 15,units=("cm"),dpi=600)
# ggsave("Regr_coefficients.jpeg",dpi=600)


