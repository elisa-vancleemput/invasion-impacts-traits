########################################################
###                                                  ###
###              Exploratory analayses               ###
###                         &                        ###
###                Linear mixed models               ###
###                                                  ###
########################################################

# Exploratory analyses:
#   1. PCA
#   2. Correlations between traits

# Linear mixed models

# In the first part (A) of this script we merge and organize the data necessary for the analyses.
# We store this data so that:
# In the second part of this script the necessary data can directly be loaded (B) and analysed (C - ...)

# Script by Elisa Van Cleemput, 2019.

#--------------------------------------------------------------------------------
# Clean workspace
rm(list=ls())

library(xlsx)
library(stringr) #str_sub, str_replace_all
library(FD)

#--------------------------------------------------------------------------------
# A) IMPORT AND ORGANIZE DATA
#----------------------------
# ---- 1) Per plot: import functioning and CWM values
data_functioning <- function (wd1,Invader) {
  setwd(wd1)
  dataset_imp<-read.csv("Plots_x_functions_Impatiens.csv",sep=";",header=TRUE)
  dataset_imp$Species <- "Impatiens"
  # Delete I8-3 en I8-4 because of different management policy (mowing) and therefore different species
  dataset_imp<-dataset_imp[!(dataset_imp$Plot=="I8-3" | dataset_imp$Plot=="I8-4"),]
  dataset_sol<-read.csv("Plots_x_functions_Solidago.csv",sep=";",header=TRUE)
  dataset_sol$Species <- "Solidago"
  
  if (Invader == "Impatiens"){
    dataset_functioning <- dataset_imp
  } else if (Invader == "Solidago") {
    dataset_functioning <- dataset_sol
  } else if (Invader == "both") {
    dataset_functioning <- rbind(dataset_imp,dataset_sol)
  }
  
  dataset_functioning$Plot <- str_replace_all(dataset_functioning[,"Plot"],"-",".")
  
  colnames(dataset_functioning)[colnames(dataset_functioning)=="C"] <- "C_soil"
  colnames(dataset_functioning)[colnames(dataset_functioning)=="N"] <- "N_soil"
  colnames(dataset_functioning)[colnames(dataset_functioning)=="CN"] <- "CN_soil"
  colnames(dataset_functioning)[colnames(dataset_functioning)=="S"] <- "tea_S"
  colnames(dataset_functioning)[colnames(dataset_functioning)=="k"] <- "tea_k"
  
  return(dataset_functioning)
}
data_CWM <- function(wd2,approach) {
  setwd(wd2)
  dataset_CWM <- read.csv(paste("CWM_",".csv",sep=approach),sep=";",header=T)
  
  return(dataset_CWM)
} 
data_FD <- function(wd2,approach) {
  setwd(wd2)
  dataset_FD1<- read.csv(paste("FD_",".csv",sep=paste(approach,"1",sep="_")),sep=";",header=T)
  dataset_FD2<- read.csv(paste("FD_",".csv",sep=paste(approach,"2",sep="_")),sep=";",header=T)

  colnames(dataset_FD1) <- paste(colnames(dataset_FD1),"1",sep="_")
  colnames(dataset_FD2) <- paste(colnames(dataset_FD2),"2",sep="_")
  
  dataset_FD <- cbind(dataset_FD1,dataset_FD2[,-c(1,2)])
  colnames(dataset_FD)[c(1,2)] <- c("Site", "Plot")
  
  return(dataset_FD)
}
# Combine the plot's data on functioning and traits
data_functioning_CWM_FD <- function(wd1,wd2,Invader, traitset, jumpcorr, approach){
  data_funct_Invader <- data_functioning(wd1, Invader)
  data_CWM_approach <- data_CWM(wd2, approach)
  data_FD_approach <- data_FD(wd2, approach)
  
  # find Plot ID's which are common in both datasets
  # retain only CWM/FD of those plots corresponding to the selected Invader
  # combine data on CWM/FD and functioning in 1 data frame --> make sure they have the same plot order first!
  intersection_plots <- intersect(data_funct_Invader[,"Plot"],data_CWM_approach[,"Plot"])
  data_funct_Invader<-data_funct_Invader[data_funct_Invader[,"Plot"] %in% intersection_plots, ]
  data_CWM_approach<-data_CWM_approach[data_CWM_approach[,"Plot"] %in% intersection_plots, ]
  data_CWM_approach <- data_CWM_approach[order(match(data_CWM_approach[,"Plot"],data_funct_Invader[,"Plot"])),]
  data_FD_approach<-data_FD_approach[data_FD_approach[,"Plot"] %in% intersection_plots, ]
  data_FD_approach <- data_FD_approach[order(match(data_FD_approach[,"Plot"],data_funct_Invader[,"Plot"])),]
  
  data_all <- cbind(data_funct_Invader[,"Species"],data_funct_Invader[,-ncol(data_funct_Invader)],
                    data_CWM_approach[,c(3:ncol(data_CWM_approach))],
                    data_FD_approach[,c(3:ncol(data_FD_approach))])

  colnames(data_all)[1] <- "Species"
  return(data_all)
}

wd1 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data"
wd2 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD"
data_both_PFT <- data_functioning_CWM_FD(wd1, wd2, Invader="both", approach="PFTraits")
# Only retain functions of interest
data_both_PFT <- data_both_PFT[,-which(names(data_both_PFT) %in% c("CN_soil","MassLoss_Green","MassLoss_Rooibos"))]
data_both_POT <- data_functioning_CWM_FD(wd1, wd2, Invader="both", approach="POTraits")
data_both_POT <- data_both_POT[,-which(names(data_both_POT) %in% c("CN_soil","MassLoss_Green","MassLoss_Rooibos"))]

# ---- 2) Remove I1.3 because missing biomass samples lead to underestimated biomass in the plot
data_both_PFT <- data_both_PFT[-which(data_both_PFT[,"Plot"] == "I1.3"),]
data_both_POT <- data_both_POT[-which(data_both_POT[,"Plot"] == "I1.3"),]

# ---- 3) Change unit of soil characteristics
data_both_PFT$C_soil <- data_both_PFT$C_soil * 10 # from % to mg/g
data_both_PFT$N_soil <- data_both_PFT$N_soil * 10 # from % to mg/g
data_both_PFT$Polsen <- data_both_PFT$Polsen * 10 # from mg/100g to mg/kg

# ---- 4)
#tea_S and tea_k have each 1 missing value
#Impute missing data with mice package
library(mice)
tempData <- mice(data_both_PFT,m=5,maxit=100,meth='pmm',seed=500,remove_collinear = FALSE)
data_both_PFT <- complete(tempData,1)
tempData <- mice(data_both_POT,m=5,maxit=100,meth='pmm',seed=500,remove_collinear = FALSE)
data_both_POT <- complete(tempData,1)

# ---- 5) Log transform ecosystem functions
data_both_PFT$logDB<-log10(data_both_PFT$DB)
data_both_PFT$logC<-log10(data_both_PFT$C_soil)
data_both_PFT$logN<-log10(data_both_PFT$N_soil)
data_both_PFT$logP<-log10(data_both_PFT$Polsen)
data_both_POT$logDB<-log10(data_both_POT$DB)
data_both_POT$logC<-log10(data_both_POT$C_soil)
data_both_POT$logN<-log10(data_both_POT$N_soil)
data_both_POT$logP<-log10(data_both_POT$Polsen)

# ---- 6) Save the data
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD")
write.table(data_both_PFT,file="Data_CWM_FD_PFTraits_Funct.csv",row.names = F)
write.table(data_both_POT,file="Data_CWM_FD_POTraits_Funct.csv",row.names = F)

# ---- 7) We will also perform analyses on CWM values of native species only
#  Repeat part of the above steps to create this data and stores it
######## We can do the same for CWM trait values of native species
wd1 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data"
wd2 <- "C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD"
data_IAScover <- data_functioning(wd1,"both")

# ---- 7a) Conventionally measured traits
data_CWM_NAT <- data_CWM(wd2, "PFTraits_NAT")
# We deleted two plots during calculation of CWM_NAT
intersection_plots <- intersect(data_IAScover[,"Plot"],data_CWM_NAT[,"Plot"])
data_IAScover<-data_IAScover[data_IAScover[,"Plot"] %in% intersection_plots, ]
data_CWM_NAT <- data_CWM_NAT[order(match(data_CWM_NAT[,"Plot"],data_IAScover[,"Plot"])),]
colnames(data_CWM_NAT) <- paste(colnames(data_CWM_NAT),"NAT",sep="_")
data_CWM_NAT <- cbind(data_IAScover[,c("Species","Invasion","Plot","PlotDesign","Site")],
                      data_CWM_NAT[,c(3:ncol(data_CWM_NAT))])
# Remove I1.3 (because also done for full analyses)
data_CWM_NAT <- data_CWM_NAT[-which(data_CWM_NAT[,"Plot"] == "I1.3"),]

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD")
write.table(data_CWM_NAT,file="Data_CWM_NAT_PFTraits.csv",row.names = F)

# ---- 7b) Optically measured traits
data_CWM_NAT <- data_CWM(wd2, "POTraits_NAT")
intersection_plots <- intersect(data_IAScover[,"Plot"],data_CWM_NAT[,"Plot"])
data_IAScover<-data_IAScover[data_IAScover[,"Plot"] %in% intersection_plots, ]
data_CWM_NAT <- data_CWM_NAT[order(match(data_CWM_NAT[,"Plot"],data_IAScover[,"Plot"])),]
colnames(data_CWM_NAT) <- paste(colnames(data_CWM_NAT),"NAT",sep="_")
data_CWM_NAT <- cbind(data_IAScover[,c("Species","Invasion","Plot","PlotDesign","Site")],
                      data_CWM_NAT[,c(3:ncol(data_CWM_NAT))])
data_CWM_NAT <- data_CWM_NAT[-which(data_CWM_NAT[,"Plot"] == "I1.3"),]
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD")
write.table(data_CWM_NAT,file="Data_CWM_NAT_POTraits.csv",row.names = F)

#--------------------------------------------------------------------------------
# B) Immediately import the data
#----------------------------
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/CWM and FD")
data_both_PFT <- read.table("Data_CWM_FD_PFTraits_Funct.csv",sep=";", header=T)
data_both_POT <- read.table("Data_CWM_FD_POTraits_Funct.csv",sep=";", header=T)

data_both_PFT_NAT <- read.table("Data_CWM_NAT_PFTraits.csv", header=T)
data_both_POT_NAT <- read.table("Data_CWM_NAT_POTraits.csv", header=T)

traits <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")
traits_NAT <- c("SLA_NAT","LDMC_NAT","LNC_mass_NAT","LCaMgC_mass_log_NAT","Height_NAT")
names(data_both_PFT_NAT)[names(data_both_PFT_NAT) == traits_NAT] <- traits
names(data_both_POT_NAT)[names(data_both_POT_NAT) == traits_NAT] <- traits
#--------------------------------------------------------------------------------
# C) Compare control plots of sites invaded by the two IAS
#----------------------------

library(reshape) # melt
library(ggplot2)

PFT_control <- data_both_PFT[which(data_both_PFT[,"Invasion"]<5),]
POT_control <- data_both_POT[which(data_both_POT[,"Invasion"]<5),]

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/Control plots comparison")

cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
functions=c("logDB","logP","logN","logC","tea_S","tea_k")
labels<- c(logDB  = "Dry biomass", logP = "soil P", logN = "soil N",logC = "soil C", tea_S = "TBI S", tea_k =  "TBI k")

meltData_PFT_functions<-melt(PFT_control[,c("Species",functions)],id.vars="Species",na.rm=T)
b <- ggplot(meltData_PFT_functions, aes(Species, value, color=Species)) +
  geom_boxplot() + 
  geom_jitter()+
  facet_wrap(~variable, scale="free", labeller=labeller(variable = labels))+
  scale_x_discrete(labels=c("I. glandulifera","S. gigantea"))+
  scale_color_manual(values=cols)+
  theme_bw()+
  theme(legend.position="none",axis.title = element_blank())
print(b)
ggsave("Boxplots_control_plots.jpeg", width = 20,height = 15,units=("cm"),dpi=600)

# Unpaired two samples t-test (parametric) or Welch t-test
# http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r
# Assumptions
library(car) # qpplot
with(PFT_control, shapiro.test(logDB[Species == "Impatiens"]))
with(PFT_control, shapiro.test(logDB[Species == "Solidago"]))
with(PFT_control, shapiro.test(tea_S[Species == "Impatiens"]))
with(PFT_control, shapiro.test(tea_S[Species == "Solidago"]))
with(PFT_control, shapiro.test(tea_k[Species == "Impatiens"])) # p = 0.008
with(PFT_control, shapiro.test(tea_k[Species == "Solidago"]))
with(PFT_control, shapiro.test(logP[Species == "Impatiens"]))
with(PFT_control, shapiro.test(logP[Species == "Solidago"])) # p = 0.03
with(PFT_control, shapiro.test(logC[Species == "Impatiens"]))
with(PFT_control, shapiro.test(logC[Species == "Solidago"]))
with(PFT_control, shapiro.test(logN[Species == "Impatiens"])) # p = 0.04
with(PFT_control, shapiro.test(logN[Species == "Solidago"]))
var.test(logDB ~ Species, data = PFT_control)
var.test(tea_S ~ Species, data = PFT_control)
var.test(tea_k ~ Species, data = PFT_control)
var.test(logP ~ Species, data = PFT_control)
var.test(logC ~ Species, data = PFT_control)
var.test(logN ~ Species, data = PFT_control)

t.test(logDB ~ Species, data = PFT_control, alternative = "two.sided", var.equal = TRUE)
t.test(logP ~ Species, data = PFT_control, alternative = "two.sided", var.equal = TRUE)
t.test(logN ~ Species, data = PFT_control, alternative = "two.sided", var.equal = TRUE)
t.test(logC ~ Species, data = PFT_control, alternative = "two.sided", var.equal = TRUE)
t.test(tea_S ~ Species, data = PFT_control, alternative = "two.sided", var.equal = TRUE)
t.test(tea_k ~ Species, data = PFT_control, alternative = "two.sided", var.equal = TRUE)

#--------------------------------------------------------------------------------
# D) VISUALISE VARIATION IN THE DATA: CORRELATIONS, PCA
#----------------------------

# ---- 1) CORRELATIONS (not reported in paper)
library(gclus) # dmat.color
panel.cor.custom <- function(x, y, digits = 2, cex.cor,  ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  # txt <- paste("r= ", txt, sep = "")
  if(missing(cex.cor))  cex <- 0.5/strwidth(txt)
  txt <- text(0.5, 0.5, txt)
  # txt <- text(0.5, 0.5, txt,cex = cex*(abs(r)+0.4))
  # txt <- text(0.5, 0.5, txt,cex = cex)
  
  #p-value calculation --> kan je ook weglaten
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  # txt2 <- paste("p = ", txt2, sep = "")
  txt2 <- paste("(", ")", sep = txt2)
  if(p<0.01) txt2 <- paste("(", "0.01)", sep = "< ")
  text(0.5, 0.2, txt2)
  # txt <- text(0.5, 0.2, txt2,cex = cex*(abs(r)+0.2))
}
corr_plot <- function(dataset, invader, name, labels) {
  
  colnames(dataset) <- labels
  
  # Examine correlation among the functional traits and functions: Scatter plot matrix
  dataset.r <- abs(cor(dataset)) # get correlations
  dataset.col <- dmat.color(dataset.r) # get colors
  # dataset.o <- order.single(dataset.r) # Optional: reorder variables so those with highest correlation are closest to the diagonal
  png(paste(invader,"_correlations.png",sep=name),width = 20,height = 20,units=("cm"),res=600)
  # x11(width=720,height=720)
  # cpairs(dataset, dataset.o, panel.colors=dataset.col, gap=.5,
  #        main="Pearson correlation" , upper.panel=panel.cor.custom)
  cpairs(dataset, order=NULL, panel.colors=dataset.col, gap=.5,
          upper.panel=panel.cor.custom) #main="Pearson correlation"
  dev.off()
}

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/Correlation_traits_function")
var_select <- c("SLA", "LDMC", "LNC_mass", "LCaMgC_mass_log", "Height","logDB","tea_S","tea_k","logP","logC","logN")
labels <- c("SLA", "LDMC", "Leaf N", "Leaf Ca+Mg", "Height","Dry biomass"," TBI S","TBI k","soil P","soil C","soil N")
corr_plot(data_both_PFT[,var_select], "both", name = "PFT_select", labels)
corr_plot(data_both_POT[,var_select], "both", name = "POT_select", labels)

var_select <- c("SLA", "LDMC", "LNC_mass", "LCaMgC_mass_log", "Height")
labels <- c("SLA", "LDMC", "Leaf N", "Leaf Ca+Mg (log)", "Height")
corr_plot(data_both_PFT[,var_select], "both", name = "PFT_select_CWM", labels)

# ---- 2) PCA
library(gginnards) # move_layers # https://cran.r-project.org/web/packages/gginnards/vignettes/user-guide-2.html
myggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                        obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                        ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                        alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                        varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                        col_arrows="gray20",var.labels=FALSE,
                        ...) {
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  library(ggrepel) #geom_text_repel
  
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable names
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    # g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
    #                                        xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
    #                                                                                               "picas")), color = muted("red"))
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = col_arrows)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.labels) {
    # g <- g + geom_text(data = df.v, aes(label = varname, 
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust), 
    #                    color = "darkred", size = varname.size)
    # g <- g + geom_text(data = df.v, aes(label = varname,
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust),
    #                    color = col_arrows, size = varname.size)
    g <- g + geom_text_repel(data = df.v, aes(label = varname,
                                              x = xvar, y = yvar, angle = angle),
                             color = col_arrows, size = varname.size)
    
  }
  return(g)
}

# ---- 2a) PCA on plot level
pca_traits <- function(data, turn, invader, traits, traits_labels, name) {
  inv_grade <- data[,"Invasion"]
  species <- data[,"Species"]

  data_pca <- data[,traits]
  colnames(data_pca) <- traits_labels
  pca<-prcomp(data_pca,scale=TRUE,center=TRUE)
  
  pca$x[,1] <- pca$x[,1] * turn[1]
  pca$x[,2] <- pca$x[,2] * turn[2]
  pca$x[,3] <- pca$x[,3] * turn[3]
  # for plotting purposes, turn axis if wanted
  pca$rotation[,1] <- pca$rotation[,1] * turn[1]
  pca$rotation[,2] <- pca$rotation[,2] * turn[2]
  pca$rotation[,3] <- pca$rotation[,3] * turn[3]
  
  pca_Imp <- cbind(pca$x[which(species == "Impatiens"),c(1,2)],inv_grade[which(species == "Impatiens")])
  colnames(pca_Imp)[3] <- "Impatiens"
  pca_Imp <- as.data.frame(pca_Imp)
  pca_Sol <- cbind(pca$x[which(species == "Solidago"),c(1,2)],inv_grade[which(species == "Solidago")])
  colnames(pca_Sol)[3] <- "Solidago"
  pca_Sol <- as.data.frame(pca_Sol)
  
  # summary(pca)
  # x11(width=1500,height=720)
  par(mfrow=c(1,2))
  screeplot(pca) # idem as plot(pca) # It seems reasonable to withhold 3 PC's (also when looking at proportion of variance they explain)
  plot(pca,type='l')
  abline(h=1, col="blue")
  summary(pca)
  # print(pca)
  
  x11(width=1500,height=720)
  g <- myggbiplot(pca,pc.biplot = TRUE, obs.scale = 1, var.scale = 1,var.labels=T,#choices = 2:3,
                  ellipse = FALSE, circle = FALSE, var.axes=T) #, groups=groups labels=label, #,groups=label2) #grouping only works PFT_intra dataset!
  g <- g + 
    geom_point(data=pca_Sol, aes(x=PC1, y=PC2, fill=Solidago), color=rep("white",nrow(pca_Sol)), size = 3.5, shape=21) +
    scale_fill_gradient(low = 'wheat1', high = 'goldenrod3', name = "S. gigantea % cover")+
    geom_point(data=pca_Imp, aes(x=PC1, y=PC2, colour=Impatiens),size = 3, shape = 16) +
    scale_color_gradient(low = 'mistyrose', high = 'deeppink3', name = "I. glandulifera % cover")+
    scale_x_continuous(breaks=seq(floor(min(pca$x[,1]))+floor(min(pca$x[,1]))%%2,ceiling(max(pca$x[,1]))+ceiling(max(pca$x[,1]))%%2,2))+
    scale_y_continuous(breaks=seq(floor(min(pca$x[,2]))+floor(min(pca$x[,2]))%%2,ceiling(max(pca$x[,2]))+ceiling(max(pca$x[,2]))%%2,2))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(legend.title = element_text(face = "italic"))
  g_reorder <- move_layers(g,"GeomPoint",position="bottom")
  print(g_reorder)
  ggsave(paste("Invasion_CWM_pca",".jpeg",sep=name),width = 15,height = 10,units=("cm"),dpi=600)
 
  return(pca)
} 

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/CWM_trait_pca")
traits_select <- c("SLA", "LDMC", "LNC_mass", "LCaMgC_mass_log", "Height")
traits_labels <- c("SLA", "LDMC", "Leaf N", "Leaf Ca+Mg (log)", "Height")
pca_PFT_select <- pca_traits(data_both_PFT_mass,turn=c(-1,-1,1), "both", traits_select, traits_labels, name="_PFT_mass")
pca_POT_select <- pca_traits(data_both_POT_mass, turn=c(1,-1,1), "both", traits_select, traits_labels, name="_POT_mass")

# Procrustes analysis between the two PCA's
library(vegan)
nb_comp = 2
proc <- procrustes(pca_PFT_select$x[,1:nb_comp], pca_POT_select$x[,1:nb_comp],symmetric=F)
prot <- protest(pca_PFT_select$x[,1:nb_comp], pca_POT_select$x[,1:nb_comp])

# Note:
#' Both 'prcomp()' and 'princomp()' give the square roots of the eigenvalues ("standard deviations").
#' Also, both functions do not give the loadings (= matrix A) but the eigenvectors (= matrix U):
#' sd >< eigenvalues --> square them
#' Eigenvectors (= matrix U): pca$rotation
#' Component matrix (= matrix A): pca$rotation %*% diag(pca$sdev)
#' Or: tcrossprod(pca$rotation, diag(pca$sdev))
pca_PFT_select$loadings <- pca_PFT_select$rotation %*% diag(pca_PFT_select$sdev)
pca_POT_select$loadings <- pca_POT_select$rotation %*% diag(pca_POT_select$sdev)


# ---- 2a) PCA on species level (not reported in paper)
pca_traits_species <- function(data, turn, traits_select, traits_labels, vis, name){
  data_pca <- data[,traits_select]
  colnames(data_pca) <- traits_labels
  
  pca<-prcomp(data_pca,scale=TRUE,center=TRUE)
  
  pca$x[,1] <- pca$x[,1] * turn[1]
  pca$x[,2] <- pca$x[,2] * turn[2]
  pca$x[,3] <- pca$x[,3] * turn[3]
  # for plotting purposes, turn axis if wanted
  pca$rotation[,1] <- pca$rotation[,1] * turn[1]
  pca$rotation[,2] <- pca$rotation[,2] * turn[2]
  pca$rotation[,3] <- pca$rotation[,3] * turn[3]
  
  palette_color<-c("#73E0B0","#73E0B0","#73E0B0","#E06762","#E06762","#E06762","#CCBAD8","#CCBAD8","#CCBAD8",
                   "#E0D64D","#E0D64D","#E0D64D","#804BDD","#804BDD","#804BDD","#D3DFCD","#D3DFCD","#D3DFCD",
                   "#8AE659","#8AE659","#8AE659","#7AC8D7","#7AC8D7","#7AC8D7","#DC8DC8","#DC8DC8","#DC8DC8",
                   "#C3DB8D","#C3DB8D","#C3DB8D","#7A88D4","#7A88D4","#7A88D4","#DA50CB","#DA50CB","#DA50CB",
                   "#D4A67E","#D4A67E","#D4A67E")
  palette_shape <- rep(c(15,16,17,15,16,17,15,16,17),5)
  # x11(width=1500,height=720)
  win.graph()
  g <- myggbiplot(pca, pc.biplot=T, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = F, circle = F, var.axes=T,alpha=0, # specify groups = vis and remove alpha = 0 when no shape differences
                  var.labels=T) + 
    geom_point(aes(colour=vis, shape=vis), size = 3) +
    scale_color_manual(name='Species',values=palette_color) +
    scale_shape_manual(name='Species',values=palette_shape)+
    theme(legend.direction = 'vertical', 
          legend.position = 'right') +
    theme(legend.position="none")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  g_reorder <- move_layers(g,idx=4L,position="bottom")
  print(g_reorder)
  
  ggsave(paste("Species_level_pca",".jpeg",sep=name),width = 20,height = 15,units=("cm"),dpi=600)
  
  return(pca)
  
}

### Conventionally measured traits 
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data")
PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
PFT$LNC_mass <- 10 * PFT$LNC # % dry mass to mg/g
PFT$LCaC_mass <- 10 * PFT$Ca # % dry mass to mg/g
PFT$LMgC_mass <- 10 * PFT$Mg # % dry mass to mg/g
PFT$LCaMgC_mass <- PFT$LCaC_mass + PFT$LMgC_mass
PFT$LCaC_mass_log <- log10(PFT$LCaC_mass)
PFT$LMgC_mass_log <- log10(PFT$LMgC_mass)
PFT$LCaMgC_mass_log <- PFT$LCaC_mass_log + PFT$LMgC_mass_log
PFT <- PFT[,c("Species.code","Site",c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height"))]

traits_select <- c("SLA", "LDMC", "LNC_mass", "LCaMgC_mass_log", "Height")
traits_labels <- c("SLA", "LDMC", "Leaf N", "Leaf Ca+Mg (log)", "Height")
vis <- PFT[,"Species.code"]
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/CWM_trait_pca")
pca_PFT_species <- pca_traits_species(PFT, turn=c(-1,-1,1), traits_select, traits_labels, vis, name="_PFTraits")


### Optically measured traits 
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/c) PLSR best models overview/")
PFT_pred_plsr <- read.csv("PLSR_best_predictions_cal.csv",sep=";",header=TRUE, row.names=1)
PFT_pred_plsr$LCaMgC_mass_log <- PFT_pred_plsr$LCaC_mass_log + PFT_pred_plsr$LMgC_mass_log
PFT_pred_plsr_acc <- read.csv("PLSR_best_predictions_accuracy.csv",sep=";",header=TRUE, row.names=1)
PFT_pred_plsr$Species.code <- substr(rownames(PFT_pred_plsr), 1, 6)
vis <- PFT_pred_plsr[,"Species.code"]

traits_select <- c("SLA", "LDMC", "LNC_mass", "LCaMgC_mass_log", "Height")
traits_labels <- c("Optically measured SLA", "Optically measured LDMC", "Optically measured Leaf N",
                   "Optically measured Lraf Ca+Mg (log)", "Optically measured Height")
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/CWM_trait_pca")
pca_POT_species <- pca_traits_species(PFT_pred_plsr, turn=c(-1,-1,1),traits_select,traits_labels,vis,name="_POTraits")

# Procrustes analysis between the two PCA's
nb_comp = 2
proc_species <- procrustes(pca_PFT_species$x[,1:nb_comp], pca_POT_species$x[,1:nb_comp],symmetric=F)
prot_species <- protest(pca_PFT_species$x[,1:nb_comp], pca_POT_species$x[,1:nb_comp])

#--------------------------------------------------------------------------------
# E) Mixed models: functions
#----------------------------

# library(nlme) # --> lme function
library(lme4) # --> lmer function
# library(lmerTest) # lmer function from lme4 package + p-values
# detach("package:lmerTest", unload=TRUE)
library(MuMIn) # r.squaredGLMM, dredge

# We will set REML=F in order to compare fixed effects with likelihood ratio test
# http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf
# https://www.researchgate.net/post/REML_FALSE_versus_REML_TRUE_lme4_package_in_R-any_thoughts

# Assumptions of residuals
#1. Linearity --> residual plot or scatter plot: no nonlinear or curvy pattern (model with fixed effects only or including random effects) 
#2. Absence of colinearity --> regression between residuals and predictor variables
#3. Homoscedasticity --> residual plot: blob-like; boxplot
#4. Normality of residuals --> histogram or QQplot of residuals
#5. Absence of influential data points (outliers) --> residual plot, boxplot, dot plot of residuals
#6. Normality of random effects
library(car) # qqPlot and Anova
LMM_assumptions_fig <- function (model, response, rand){
  # x11(width=3000,height=1700)
  jpeg(paste(paste(paste("Invasion_",response,sep="impacts_"),rand,sep="_Rand"),".jpeg",sep="_Assumptions"),width = 20,height = 15,units=("cm"),res=600)
  
  par(mfrow=c(2,3))
  
  # Histogram of residuals
  hist(residuals(model),col="darkgray",freq=TRUE,breaks=20,density=70,main="Histogram of residuals",xlab="Raw residuals",ylab="Absolute frequencies")
  # hist(residuals(model),col="darkgray",freq=FALSE,breaks=20,density=70,main="Histogram of residuals",xlab="Raw residuals",ylab="Relative density")
  # curve(dnorm(x,mean=mean(residuals(model)),sd=sd(residuals(model))),col="darkblue",lwd=2,add=TRUE,yaxt="n")
  
  # Boxplot of residuals
  boxplot(residuals(model),main="Boxplot of residuals",ylab="Residuals")
  
  # QQPlot
  # http://www.sthda.com/english/wiki/normality-test-in-r
  # print(ggqqplot(residuals(model),main="Normal Q-Q Plot for residuals",xlab="Theoretical quantiles",ylab="Sample quantiles"))
  shapiro.test(residuals(model))
  # qqnorm(residuals(model))
  # qqline(residuals(model))
  qqPlot(residuals(model))
  text(shapiro.test(residuals(model))$p)
  legend("topleft", paste("Shapiro p =", round(shapiro.test(residuals(model))$p,2)), bty="n") 
  
  # Residual plot
  plot(fitted(model),residuals(model),xlab="Predicted",ylab="Raw Residuals")
  abline(h=0,lty=2,col = "gray60") #adds horizontal line for y=0
  # lines(smooth.spline(fitted(model),residuals(model))) #Fits a cubic smoothing spline to the supplied data
  
  # Checking Colinearity 
  plot(model@frame$Invasion,residuals(model),xlab="Invasion %",ylab="Raw Residuals")
  # plot(attributes(model)$frame[,1],residuals(model),xlab=names(attributes(model)$frame)[1],ylab="Raw Residuals")
  # plot(model$data[,"Invasion"],residuals(model),xlab="Invasion %",ylab="Raw Residuals")
  
  pvalues <- Anova(model, test.statistic = c("F"))
  
  plot.new()
  legend("center",bty="n",legend=c(paste("Function =",names(attributes(model)$frame)[1]),
                                   paste("R²m =",round(r.squaredGLMM(model)[1],2)),
                                   paste("R²c =",round(r.squaredGLMM(model)[2],2)),
                                   # paste("AIC =",round(AIC(model),2)), 
                                   # paste("AICc =",round(AICc(model),2)),
                                   paste("sd Site =", round(as.data.frame(VarCorr(model))[1,"sdcor"],2)),
                                   paste("sd Residual =", round(as.data.frame(VarCorr(model))[2,"sdcor"],2)),
                                   paste("F (p) Study species =",paste(round(pvalues["Species",1],2),paste(round(pvalues["Species",4],2), ")",sep=""),sep=" (")),
                                   paste("F (p) Invasion =",paste(round(pvalues["Invasion",1],2),paste(round(pvalues["Invasion",4],2), ")",sep=""),sep=" (")),
                                   paste("F (p) Species:Invasion =",paste(round(pvalues["Species:Invasion",1],2),paste(round(pvalues["Species:Invasion",4],2), ")",sep=""),sep=" ("))),
                          cex=0.9)

  dev.off()
}

library(ggplot2)
LMM_fig <- function(modeldata, inv, response, rand, mm, xtext, ytext, letter){
  if (length(predict(mm)) != dim(modeldata)[1]){
    modeldata<- modeldata[complete.cases(modeldata),]
  }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Predictions", values=c("random"=0.2, "overall"=2)) +
    scale_color_manual(values=cols)+
    xlab(xtext)+
    ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(legend.position="none")
  print(p)
  ggsave(paste(paste(paste(inv,response,sep="_impacts_"),rand,sep="_Rand"),".jpeg",sep=""),width = 15,height = 15,units=("cm"),dpi=600)
  return(p)
}
LMM_fig_noimpact <- function(modeldata, inv, response, rand, mm, xtext, ytext, letter){
  if (length(predict(mm)) != dim(modeldata)[1]){
    modeldata<- modeldata[complete.cases(modeldata),]
  }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    # geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    # geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Predictions", values=c("random"=0.2, "overall"=2)) +
    scale_color_manual(values=cols)+
    xlab(xtext)+
    ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(legend.position="none")
  print(p)
  ggsave(paste(paste(paste(inv,response,sep="_doesnotimpact_"),rand,sep="_Rand"),".jpeg",sep=""),width = 15,height = 15,units=("cm"),dpi=600)
  return(p)
}
LMM_fig_minmax <- function(modeldata, inv, response, rand, mm, ymin, ymax, xtext, ytext, letter){
  # if (length(predict(mm)) != dim(modeldata)[1]){
  #   modeldata<- modeldata[complete.cases(modeldata),]
  # }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Predictions", values=c("random"=0.2, "overall"=2)) +
    scale_color_manual(values=cols)+
    # scale_y_continuous(limits=c(ymin, ymax), expand = c(0,0))+
    coord_cartesian(ylim=c(ymin, ymax))+
    xlab(xtext)+
    ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(legend.position="none")
  print(p)
  ggsave(paste(paste(paste(inv,response,sep="_impacts_"),rand,sep="_Rand"),".jpeg",sep=""),width = 15,height = 15,units=("cm"),dpi=600)
  # return(p)  
}
LMM_fig_minmax_noY <- function(modeldata, inv, response, rand, mm, ymin, ymax, xtext, ytext, letter){
  # if (length(predict(mm)) != dim(modeldata)[1]){
  #   modeldata<- modeldata[complete.cases(modeldata),]
  # }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Predictions", values=c("random"=0.2, "overall"=2)) +
    scale_color_manual(values=cols)+
    # scale_y_continuous(limits=c(ymin, ymax), expand = c(0,0))+
    coord_cartesian(ylim=c(ymin, ymax))+
    xlab(xtext)+
    # ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())+
    theme(legend.position="none")
  print(p)
  ggsave(paste(paste(paste(inv,response,sep="_impacts_"),rand,sep="_Rand"),"_noY.jpeg",sep=""),width = 13,height = 15,units=("cm"),dpi=600)
  # return(p)  
}
LMM_fig_minmax_noimpact <- function(modeldata, inv, response, rand, mm, ymin, ymax, xtext, ytext, letter){
  # if (length(predict(mm)) != dim(modeldata)[1]){
  #   modeldata<- modeldata[complete.cases(modeldata),]
  # }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    # geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    # geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Predictions", values=c("random"=0.2, "overall"=2)) +
    scale_color_manual(values=cols)+
    # scale_y_continuous(limits=c(ymin, ymax), expand = c(0,0))+
    coord_cartesian(ylim=c(ymin, ymax))+
    xlab(xtext)+
    ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(legend.position="none")
  print(p)
  ggsave(paste(paste(paste(inv,response,sep="_doesnotimpact_"),rand,sep="_Rand"),".jpeg",sep=""),width = 15,height = 15,units=("cm"),dpi=600)
  # return(p)  
}
LMM_fig_minmax_noimpact_noY <- function(modeldata, inv, response, rand, mm, ymin, ymax, xtext, ytext, letter){
  # if (length(predict(mm)) != dim(modeldata)[1]){
  #   modeldata<- modeldata[complete.cases(modeldata),]
  # }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    # geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    # geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Predictions", values=c("random"=0.2, "overall"=2)) +
    scale_color_manual(values=cols)+
    # scale_y_continuous(limits=c(ymin, ymax), expand = c(0,0))+
    coord_cartesian(ylim=c(ymin, ymax))+
    xlab(xtext)+
    # ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())+
    theme(legend.position="none")
  print(p)
  ggsave(paste(paste(paste(inv,response,sep="_doesnotimpact_"),rand,sep="_Rand"),"_noY.jpeg",sep=""),width = 13,height = 15,units=("cm"),dpi=600)
  # return(p)  
}


Boxplot_fig <- function(modeldata, inv, response, rand, mm, ytext, letter){
  if (length(predict(mm)) != dim(modeldata)[1]){
    modeldata<- modeldata[complete.cases(modeldata),]
  }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  q <- ggplot(modeldata, aes_string(x="Species", y=response, color="Species")) +
    geom_boxplot(size=2) +
    geom_jitter()+
    scale_color_manual(values=cols)+
    scale_x_discrete(labels=c("I. glandulifera","S. gigantea"))+
    ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(axis.title.x = element_text(colour = "white")) +
    theme(axis.text.x = element_text(face = "italic")) +
    # theme(axis.title.x = element_text(colour = "white"),axis.text.x = element_text(colour = "white"),
    #       axis.ticks.x=element_blank())+
    theme(legend.position="none")
  print(q)
  ggsave(paste(paste(paste(inv,response,sep="_impacts_"),rand,sep="_Boxplot_Rand"),".jpeg",sep=""),width = 15,height = 15,units=("cm"),dpi=600)
  return(q)  
}
Boxplot_fig_minmax <- function(modeldata, inv, response, rand, mm, ymin, ymax, ytext, letter){
  # if (length(predict(mm)) != dim(modeldata)[1]){
  #   modeldata<- modeldata[complete.cases(modeldata),]
  # }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  q <- ggplot(modeldata, aes_string(x="Species", y=response, color="Species")) +
    geom_boxplot(size=2) +
    geom_jitter()+
    scale_color_manual(values=cols)+
    coord_cartesian(ylim=c(ymin, ymax))+
    scale_x_discrete(labels=c("I. glandulifera","S. gigantea"))+
    ylab(ytext)+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(axis.title.x = element_text(colour = "white")) +
    theme(axis.text.x = element_text(face = "italic")) +
    # theme(axis.title.x = element_text(colour = "white"),axis.text.x = element_text(colour = "white"),
    #       axis.ticks.x=element_blank())+
    theme(legend.position="none")
  print(q)
  ggsave(paste(paste(paste(inv,response,sep="_impacts_"),rand,sep="_Boxplot_Rand"),".jpeg",sep=""),width = 15,height = 15,units=("cm"),dpi=600)
  return(q)  
}
Boxplot_IAS_fig_minmax <- function(data, response, ymin, ymax, letter){
  cols <- c("IMPGLA" = "deeppink3", "SOLGIG" = "gold3")
  
  q <- ggplot(data, aes_string(x="Species.code", y=response)) +
    geom_boxplot(aes(color=Species.code),size=2) +
    # geom_jitter(aes(color=Species.code),size=3,width=0.1,pch=21,fill="black")+
    geom_point(position=position_jitterdodge(jitter.width=3, dodge.width = 0),
    pch=21,size=3, aes(fill=Species.code), color="white", show.legend = F) +
    # geom_point(aes(colour=id), size=12) + 
    # geom_point(shape = 1,size = 12,colour = "black")+
    scale_color_manual(values=cols)+
    scale_fill_manual(values=cols)+
    coord_cartesian(ylim=c(ymin, ymax))+
    scale_x_discrete(labels=c("I. glandulifera","S. gigantea"))+
    annotate(geom="text", label=letter, x=-Inf, y=Inf, hjust=-1, vjust=2, cex=8)+
    theme_bw(base_size=20)+
    theme(axis.title.x = element_text(colour = "white")) +
    theme(axis.text.x = element_text(face = "italic")) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())+
    # theme(axis.title.x = element_text(colour = "white"),axis.text.x = element_text(colour = "white"),
    #       axis.ticks.x=element_blank())+
    theme(legend.position="none")
  print(q)
  ggsave(paste("Boxplots_IAS_",".jpeg",sep=response),width = 10,height = 15,units=("cm"),dpi=600)
  return(q)  
}

library(cowplot) # get_legend
LMM_fig_legend_hor <- function (modeldata, inv, response, rand, mm){
  if (length(predict(mm)) != dim(modeldata)[1]){
    modeldata<- modeldata[complete.cases(modeldata),]
  }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Regression", values=c("random"=0.5, "overall"=2)) +
    scale_color_manual(name = "Invader identity", values=cols, labels=c("I. glandulifera","S. gigantea"))+
    xlab("Invader cover")+
    theme_bw(base_size=20)+
    theme(legend.direction = "vertical", 
          legend.position = "bottom",
          legend.box = "horizontal")+
    theme(legend.text = element_text(face = "italic"))
  print(p)
  legend<-get_legend(p)
  ggdraw(plot_grid(legend,ncol=1))
  # print(p)
  ggsave("Legend.jpeg",width = 18,height = 5,units=("cm"),dpi=600)
}
LMM_fig_legend_ver <- function (modeldata, inv, response, rand, mm){
  if (length(predict(mm)) != dim(modeldata)[1]){
    modeldata<- modeldata[complete.cases(modeldata),]
  }
  
  cols <- c("Impatiens" = "deeppink3", "Solidago" = "gold3")
  
  p <- ggplot(modeldata, aes_string(x=inv, y=response, group=rand, color="Species")) +
    geom_point(size=3) +
    geom_line(aes(y=predict(mm), group=modeldata[,rand], size="random")) +
    geom_line(aes(y=predict(mm, re.form=NA), size="overall")) + # color = "black"
    scale_size_manual(name="Regression", values=c("random"=0.5, "overall"=2)) +
    scale_color_manual(name = "Invader identity", values=cols, labels=c("I. glandulifera","S. gigantea"))+
    xlab("Invader % coverage")+
    theme_bw(base_size=20)+
    theme(legend.text = element_text(face = "italic"))
  print(p)
  legend<-get_legend(p)
  ggdraw(plot_grid(legend,ncol=1))
  # print(p)
  ggsave("Legend.jpeg",width = 5,height = 10,units=("cm"),dpi=600)
}

#--------------------------------------------------------------------------------
# E-1) Mixed models: ecosystem functioning
#----------------------------
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Function ifo invasion and species")
modeldata <- data_both_PFT
modeldata$Invasion_dec <- modeldata$Invasion/100

#  ------- biomass
mm_lmer_logDB_Site <- lmer(logDB ~ Species + Invasion_dec + Species*Invasion_dec +
                             (1|Site), data=modeldata, REML=T)
summary(mm_lmer_logDB_Site)
r.squaredGLMM(mm_lmer_logDB_Site)
Anova(mm_lmer_logDB_Site, test.statistic = c("F"))
LMM_assumptions_fig(mm_lmer_logDB_Site,"logDB","Site")
F.a <- LMM_fig(modeldata, "Invasion", "logDB", "Site", mm_lmer_logDB_Site, 
               xtext = "Invader cover (%)", ytext = expression("Dry biomass (g "~m^{-2}~", log scale)"), "(a)")
Fbox.a <- Boxplot_fig(modeldata, "Invasion", "logDB", "Site", mm_lmer_logDB_Site, 
                      expression("Dry biomass (g "~m^{-2}~", log scale)"), "(a)")
#  ------- P Olsen
mm_lmer_logP_Site <- lmer(logP ~ Species + Invasion_dec + Species*Invasion_dec +
                              (1|Site), data=modeldata, REML=T)
summary(mm_lmer_logP_Site)
r.squaredGLMM(mm_lmer_logP_Site)
Anova(mm_lmer_logP_Site, test.statistic = c("F"))
LMM_assumptions_fig(mm_lmer_logP_Site,"logP","Site")
F.b <- LMM_fig(modeldata, "Invasion", "logP", "Site", mm_lmer_logP_Site,
               xtext = "Invader cover (%)", expression("Soil P (mg "~kg^{-1}~", log scale)"), "(b)")
Fbox.b <- Boxplot_fig(modeldata, "Invasion", "logP", "Site", mm_lmer_logP_Site,
                      expression("Soil P (mg "~kg^{-1}~", log scale)"), "(b)")

#  ------- soil N and C
mm_lmer_logN_Site <- lmer(logN ~ Species + Invasion_dec + Species*Invasion_dec +
                            (1|Site), data=modeldata, REML=T)
summary(mm_lmer_logN_Site)
r.squaredGLMM(mm_lmer_logN_Site)
Anova(mm_lmer_logN_Site, test.statistic = c("F"))
LMM_assumptions_fig(mm_lmer_logN_Site,"logN","Site")
F.c <- LMM_fig_noimpact(modeldata, "Invasion", "logN", "Site", mm_lmer_logN_Site,
                        xtext = "Invader cover (%)", expression("Soil N (mg"~g^{-1}~", log scale)"),"(c)")
Fbox.c <- Boxplot_fig(modeldata, "Invasion", "logN", "Site", mm_lmer_logN_Site,
                      expression("Soil N (mg"~g^{-1}~", log scale)"),"(c)")

mm_lmer_logC_Site <- lmer(logC ~ Species + Invasion_dec + Species*Invasion_dec +
                          (1|Site), data=modeldata, REML=T)
summary(mm_lmer_logC_Site)
r.squaredGLMM(mm_lmer_logC_Site)
Anova(mm_lmer_logC_Site, test.statistic = c("F"))
LMM_assumptions_fig(mm_lmer_logC_Site,"logC","Site")
F.d <- LMM_fig_noimpact(modeldata, "Invasion", "logC", "Site", mm_lmer_logC_Site,
                        xtext = "Invader cover (%)", expression("Soil C (mg"~g^{-1}~", log scale)"),"(d)")
Fbox.d <- Boxplot_fig(modeldata, "Invasion", "logC", "Site", mm_lmer_logC_Site,
                      expression("Soil C (mg"~g^{-1}~", log scale)"),"(d)")

#  ------ decomposition 
mm_lmer_tea_S_Site <- lmer(tea_S ~ Species + Invasion_dec + Species*Invasion_dec +
                             (1|Site), data=modeldata, REML=T) # na.omit not necessay after applying mice
summary(mm_lmer_tea_S_Site)
r.squaredGLMM(mm_lmer_tea_S_Site)
Anova(mm_lmer_tea_S_Site, test.statistic = c("F"))
LMM_assumptions_fig(mm_lmer_tea_S_Site, "tea_S", "Site")
F.e <- LMM_fig(modeldata, "Invasion", "tea_S", "Site", mm_lmer_tea_S_Site, 
               xtext = "Invader cover (%)", ytext = expression(paste("Litter stabilisation factor ", italic("S"), "")), "(e)")
Fbox.e <-Boxplot_fig(modeldata, "Invasion", "tea_S", "Site", mm_lmer_tea_S_Site, 
                     expression(paste("Litter stabilisation factor ", italic("S"), "")), "(b)")

mm_lmer_tea_k_Site <- lmer(tea_k ~ Species + Invasion_dec + Species*Invasion_dec +
                             (1|Site), data=modeldata, REML=T)
summary(mm_lmer_tea_k_Site)
r.squaredGLMM(mm_lmer_tea_k_Site)
Anova(mm_lmer_tea_k_Site, test.statistic = c("F"))
LMM_assumptions_fig(mm_lmer_tea_k_Site, "tea_k", "Site")
F.f <- LMM_fig_noimpact(modeldata, "Invasion", "tea_k", "Site", mm_lmer_tea_k_Site, 
                        xtext = "Invader cover (%)", expression(paste("Litter decomposition rate ", italic("k"), "")), "(f)")
Fbox.f <- Boxplot_fig(modeldata, "Invasion", "tea_k", "Site", mm_lmer_tea_k_Site, 
                      expression(paste("Litter decomposition rate ", italic("k"), "")), "(f)")


#  ------- legend
F.legend <- LMM_fig_legend_hor(modeldata, "Invasion", "logN", "Site", mm_lmer_logN_Site)

# To extract information, consult following functions for each model
# summary(mm_lmer_logP_Site)
# r.squaredGLMM(mm_lmer_logP_Site)
# Anova(mm_lmer_logP_Site, test.statistic = c("F"))


#--------------------------------------------------------------------------------
# E-2) Mixed models: CWM of all species
#----------------------------
# We want the plots of the conventionally and optically measured traits to have the same Y-scale
traits <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")
letters <- c("(a)","(b)","(c)","(d)","(e)")
# letters <- c("(a)","(d)","(g)","(j)","(m)")
ymins <- rep(0,length(traits))
ymaxs <- rep(0,length(traits))
# for (i in 1:length(traits)){
#   ymins[i] <-min(data_both_PFT_mass[,traits[i]],data_both_POT_mass[,traits[i]])
#   ymaxs[i] <-max(data_both_PFT_mass[,traits[i]],data_both_POT_mass[,traits[i]])
# }
for (i in 1:length(traits)){
  ymins[i] <-min(data_both_PFT[,traits[i]],data_both_POT[,traits[i]],
                 data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
  ymaxs[i] <-max(data_both_PFT[,traits[i]],data_both_POT[,traits[i]],
                 data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
}

# Conventionally measured traits
modeldata <- data_both_PFT
modeldata$Invasion_dec <- modeldata$Invasion/100
ytexts <- c(expression("SLA (mm "~mg^{-1}~")"),
            expression("LDMC (mg "~g^{-1}~")"),
            expression("Leaf N (mg "~g^{-1}~")"), 
            expression("Leaf Ca+Mg (mg "~g^{-1}~", log scale)"),
            "vegetative height (m)")
xtext <- "Invader cover (%)"
  
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Trait ifo invasion and species")
for (i in 1:length(traits)){
  lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion_dec + Species*Invasion_dec+ 
                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
  LMM_assumptions_fig(lmer_trait,traits[i],"Site")
  LMM_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  LMM_fig_minmax_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  Boxplot_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], ytexts[i], letters[i])

  # Print results
  print(summary(lmer_trait))
  print(r.squaredGLMM(lmer_trait))
  print(Anova(lmer_trait, test.statistic = c("F")))
}
CWM.legend <- LMM_fig_legend_ver(modeldata, "Invasion", "logN", "Site", mm_lmer_logN_Site)

# Conventioanlly measured Leaf Ca+Mg was not significantly impacted by invader cover
i=4
lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion_dec + Species*Invasion_dec+ 
                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
LMM_assumptions_fig(lmer_trait,traits[i],"Site")
LMM_fig_minmax_noimpact(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
LMM_fig_minmax_noimpact_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])

# Optically measured traits
modeldata <- data_both_POT
modeldata$Invasion_dec <- modeldata$Invasion/100
# ytexts <- c(expression("Optical SLA (mm "~mg^{-1}~")"),
#             expression("Optical  LDMC (mg "~g^{-1}~")"),
#             expression("Optical LNC (mg "~g^{-1}~")"), 
#             expression("Optical LCa+MgC (mg "~g^{-1}~", log scale)"),
#             "Optical vegetative height (m)")
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Trait ifo invasion and species - optical")
for (i in 1:length(traits)){
  lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion_dec + Species*Invasion_dec+ 
                                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
  LMM_assumptions_fig(lmer_trait,traits[i],"Site")
  LMM_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  LMM_fig_minmax_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  Boxplot_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], ytexts[i], letters[i])

  # Print results
  print(summary(lmer_trait))
  print(r.squaredGLMM(lmer_trait))
  print(Anova(lmer_trait, test.statistic = c("F")))
}

# summary(lmer_trait)
# r.squaredGLMM(lmer_trait)
# Anova(lmer_trait, test.statistic = c("F"))

#--------------------------------------------------------------------------------
# E-3) Mixed models: CWM of native species only
#----------------------------
# We want the plots of the conventionally and optically measured traits to have the same Y-scale
traits <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")
# letters <- c("(a)","(b)","(c)","(d)","(e)")
letters <- c("(f)","(g)","(h)","(i)","(j)")
ymins <- rep(0,length(traits))
ymaxs <- rep(0,length(traits))
# for (i in 1:length(traits)){
#   ymins[i] <-min(data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
#   ymaxs[i] <-max(data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
# }
for (i in 1:length(traits)){
  ymins[i] <-min(data_both_PFT[,traits[i]],data_both_POT[,traits[i]],
                 data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
  ymaxs[i] <-max(data_both_PFT[,traits[i]],data_both_POT[,traits[i]],
                 data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
}

# Conventionally measured traits 
modeldata <- data_both_PFT_NAT
modeldata$Invasion_dec <- modeldata$Invasion/100
ytexts <- c(expression("native SLA (mm "~mg^{-1}~")"),
            expression("native LDMC (mg "~g^{-1}~")"),
            expression("native Leaf N (mg "~g^{-1}~")"), 
            expression("native Leaf Ca+Mg (mg "~g^{-1}~", log scale)"),
            "native vegetative height (m)")
xtext <- "Invader cover (%)"

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Trait_NAT ifo invasion and species")
for (i in 1:length(traits)){
  lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion_dec + Species*Invasion_dec+ 
                                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
  LMM_assumptions_fig(lmer_trait,traits[i],"Site")
  LMM_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  LMM_fig_minmax_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  Boxplot_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], ytexts[i], letters[i])
  
  # Print results
  print(summary(lmer_trait))
  print(r.squaredGLMM(lmer_trait))
  print(Anova(lmer_trait, test.statistic = c("F")))
}

# Conventionally measured LDMC was not significantly impacted by invader cover
i=2
lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion + Species*Invasion+
                                    (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
LMM_assumptions_fig(lmer_trait,traits[i],"Site")
LMM_fig_minmax_noimpact(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
LMM_fig_minmax_noimpact_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])

# Optically measured traits
modeldata <- data_both_POT_NAT
modeldata$Invasion_dec <- modeldata$Invasion/100

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Trait_NAT ifo invasion and species - optical")
for (i in 1:length(traits)){
  lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion + Species*Invasion+ 
                                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
  LMM_assumptions_fig(lmer_trait,traits[i],"Site")
  LMM_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  LMM_fig_minmax_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  Boxplot_fig_minmax(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], ytexts[i], letters[i])
  
  # Print results
  print(summary(lmer_trait))
  print(r.squaredGLMM(lmer_trait))
  print(Anova(lmer_trait, test.statistic = c("F")))
}

# Optically measured LDMC and Leaf N were not significantly impacted by invader cover
i=2
lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion + Species*Invasion+
                                    (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
LMM_assumptions_fig(lmer_trait,traits[i],"Site")
LMM_fig_minmax_noimpact(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
LMM_fig_minmax_noimpact_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])

i=3
lmer_trait <- lmer(as.formula(paste(traits[i], "Species + Invasion + Species*Invasion+
                                    (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
LMM_assumptions_fig(lmer_trait,traits[i],"Site")
LMM_fig_minmax_noimpact(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
LMM_fig_minmax_noimpact_noY(modeldata, "Invasion", traits[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])

# summary(lmer_trait)
# r.squaredGLMM(lmer_trait)
# Anova(lmer_trait, test.statistic = c("F"))

#--------------------------------------------------------------------------------
# E-4) Mixed models: FD
#----------------------------
FD <- c("FDis_1","FDis_2","RaoQ_1","RaoQ_2")
letters <- c("(a)","(b)","(a)","(b)")

ymins <- rep(0,length(FD))
ymaxs <- rep(0,length(FD))
for (i in 1:length(FD)){
  ymins[i] <-min(data_both_PFT[,FD[i]],data_both_POT[,FD[i]])
  ymaxs[i] <-max(data_both_PFT[,FD[i]],data_both_POT[,FD[i]])
}
ymaxs[1] <- max(ymaxs[1],ymaxs[2]) # set FDis1 and FDis2 on the same scale
ymaxs[2] <- max(ymaxs[1],ymaxs[2]) 
ymaxs[3] <- max(ymaxs[3],ymaxs[4]) # set RaoQ1 and RaoQ2 on the same scale
ymaxs[4] <- max(ymaxs[3],ymaxs[4])
xtext <- "Invader cover (%)"
  
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Diversity ifo invasion and species")
modeldata <- data_both_PFT
modeldata$Invasion_dec <- modeldata$Invasion/100
ytexts <- c(expression("FDis"[1]),
            expression("FDis"[2]),
            expression("Rao's Q"[1]),
            expression("Rao's Q"[2]))
for (i in 1:length(FD)){
  lmer_trait <- lmer(as.formula(paste(FD[i], "Species + Invasion_dec + Species*Invasion_dec+ 
                                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
  LMM_assumptions_fig(lmer_trait,FD[i],"Site")
  LMM_fig_minmax(modeldata, "Invasion", FD[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  
  # Print results
  print(summary(lmer_trait))
  print(r.squaredGLMM(lmer_trait))
  print(Anova(lmer_trait, test.statistic = c("F")))
}
CWM.legend <- LMM_fig_legend_ver(modeldata, "Invasion", "logN", "Site", mm_lmer_logN_Site)


setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/LMM_Diversity ifo invasion and species - optical")
modeldata <- data_both_POT
modeldata$Invasion_dec <- modeldata$Invasion/100
ytexts <- c(expression("FDis"[1]~" of optically measured traits"),
            expression("FDis"[2]~" of optically measured traits"),
            expression("Rao's Q"[1]~" of optically measured traits"),
            expression("Rao's Q"[2])~" of optically measured traits")
for (i in 1:length(FD)){
  lmer_trait <- lmer(as.formula(paste(FD[i], "Species + Invasion_dec + Species*Invasion_dec+ 
                                      (1|Site)", sep=" ~ ")), data=modeldata, REML=T)
  LMM_assumptions_fig(lmer_trait,FD[i],"Site")
  LMM_fig_minmax(modeldata, "Invasion", FD[i], "Site", lmer_trait, ymins[i], ymaxs[i], xtext, ytexts[i], letters[i])
  
  # Print results
  print(summary(lmer_trait))
  print(r.squaredGLMM(lmer_trait))
  print(Anova(lmer_trait, test.statistic = c("F")))
}
#--------------------------------------------------------------------------------
# E-5) Extra boxplot figures for interpretation of E-3
#----------------------------
# For each trait we will also fit trait values of the two IAS

# We want the plots of the conventionally and optically measured traits to have the same Y-scale
traits <- c("SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")
letters <- c("(k)","(l)","(m)","(n)","(o)")
ymins <- rep(0,length(traits))
ymaxs <- rep(0,length(traits))
for (i in 1:length(traits)){
  ymins[i] <-min(data_both_PFT[,traits[i]],data_both_POT[,traits[i]],
                 data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
  ymaxs[i] <-max(data_both_PFT[,traits[i]],data_both_POT[,traits[i]],
                 data_both_PFT_NAT[,traits[i]],data_both_POT_NAT[,traits[i]])
}

### Conventionally measured traits 
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/Raw data")
PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
PFT$LNC_mass <- 10 * PFT$LNC # % dry mass to mg/g
PFT$LCaC_mass <- 10 * PFT$Ca # % dry mass to mg/g
PFT$LMgC_mass <- 10 * PFT$Mg # % dry mass to mg/g
PFT$LCaMgC_mass <- PFT$LCaC_mass + PFT$LMgC_mass
PFT$LCaC_mass_log <- log10(PFT$LCaC_mass)
PFT$LMgC_mass_log <- log10(PFT$LMgC_mass)
PFT$LCaMgC_mass_log <- PFT$LCaC_mass_log + PFT$LMgC_mass_log
PFT <- PFT[,c("Species.code","SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")]
PFT_IAS <- PFT[which(PFT$Species.code =="IMPGLA" | PFT$Species.code == "SOLGIG"),]

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/Boxplots_Traits_IAS")
for (i in 1:length(traits)){
  Boxplot_IAS_fig_minmax(PFT_IAS, response = traits[i], ymins[i], ymaxs[i], letters[i])
}


### Optically measured traits 
setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Data and analyses/PLSR results/c) PLSR best models overview/")
POT_pred_plsr <- read.csv("PLSR_best_predictions_cal.csv",sep=";",header=TRUE, row.names=1)
POT_pred_plsr$LCaMgC_mass_log <- POT_pred_plsr$LCaC_mass_log + POT_pred_plsr$LMgC_mass_log
POT_pred_plsr$Species.code <- substr(rownames(POT_pred_plsr), 1, 6)
POT <- POT_pred_plsr[,c("Species.code","SLA","LDMC","LNC_mass","LCaMgC_mass_log","Height")]
POT_IAS <- POT[which(POT$Species.code =="IMPGLA" | POT$Species.code == "SOLGIG"),]
POT_IAS$Species.code <- as.factor(POT_IAS$Species.code)

setwd("C:/Users/u0091812/Box Sync/05. Invasion impact/Figures/Boxplots_Traits_IAS - optical")
for (i in 1:length(traits)){
  Boxplot_IAS_fig_minmax(POT_IAS, response = traits[i], ymins[i], ymaxs[i], letters[i])
}

