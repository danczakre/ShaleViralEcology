### This script will perform spearman mantel correlations 

library(vegan) # Needed for permuted Procrustes
library(ecodist)
library(reshape2) # Needed to melt data for plotting
library(ggplot2) # Plotting
library(gridExtra) # Helps in the plotting process

########################
### Loading in data ###
########################

# Set directory
setwd("//path/to/ecological_modeling_results")

# S3 data
S3.bNTI = read.csv("S3_Community_bNTI_TotalCounts.csv",
                   row.names = 1)
S3.RCBC = read.csv("S3_Community_RCBC_TotalCounts.csv",
                   row.names = 1)

# Viral data
Vir.bNTI = read.csv("Viral_bNTI_NormCounts_NormUPGMA.csv",
                    row.names = 1)
Vir.U3RC = read.csv("U3_RCBC_NormCounts_NormUPGMA.csv",
                     row.names = 1)
Vir.U2RC = read.csv("U2_RCBC_NormCounts_NormUPGMA.csv",
                     row.names = 1)

# Factors
S3.factors = read.csv("S3_Community_Factors_Sheet.csv")
Vir.factors = read.csv("Viral_Model_Factors.csv")

############################
### Cleaning up the data ###
############################

# Removing unnecessaary comparisons
S3.bNTI = S3.bNTI[-grep("HT|LW", row.names(S3.bNTI)), -grep("HT|LW", row.names(S3.bNTI))]
S3.RCBC = S3.RCBC[-grep("HT|LW", row.names(S3.RCBC)), -grep("HT|LW", row.names(S3.RCBC))]
S3.factors = S3.factors[-grep("HT|LW", S3.factors$Sample_Name), ]

Vir.bNTI = Vir.bNTI[-grep("HT|LW", row.names(Vir.bNTI)), -grep("HT|LW", row.names(Vir.bNTI))]
if(length(grep("HT|LW", row.names(Vir.U2RC))) > 1){
  Vir.U2RC = Vir.U2RC[-grep("HT|LW", row.names(Vir.U2RC)), -grep("HT|LW", row.names(Vir.U2RC))]
}
Vir.factors = Vir.factors[-grep("HT|LW", Vir.factors$Sample_Name), ]

####################################
### Setting up some initial data ###
####################################

Wells = c("Utica-3", "Utica-2")
cols = data.frame(light = c("#abc2d6", "#b8e0c8"),
                  dark = c("#2171B5" , "#006D2C"))


##############################
### Comparing bNTI results ###
##############################

S3.bNTI[upper.tri(S3.bNTI)] = t(S3.bNTI)[upper.tri(S3.bNTI)] # Need to mirror both bNTI datasets
Vir.bNTI[upper.tri(Vir.bNTI)] = t(Vir.bNTI)[upper.tri(Vir.bNTI)]

# Creating dummy datasets
mant = NULL
prot = NULL
names = NULL

for(i in 1:length(Wells)){
  temp.SbNTI = S3.bNTI[grep(Wells[i], S3.factors$Well), grep(Wells[i], S3.factors$Well)]
  temp.VbNTI = Vir.bNTI[grep(Wells[i], Vir.factors$Alt_Sample), grep(Wells[i], Vir.factors$Alt_Sample)]
  
  if(Wells[i] == "Utica-3"){
    w = c(grep("FT|RT", row.names(temp.VbNTI), invert = T), grep("FT|RT", row.names(temp.VbNTI)))
    temp.VbNTI = temp.VbNTI[w, w]
  }
  
  if(any(apply(temp.VbNTI, 1, function(x) length(which(is.na(x) %in% TRUE))) == length(temp.VbNTI[1,]))){
    w = which(apply(temp.VbNTI, 1, function(x) length(which(is.na(x) %in% TRUE))) == length(temp.VbNTI[1,]))
    
    temp.SbNTI = temp.SbNTI[-w,-w]
    temp.VbNTI = temp.VbNTI[-w,-w]
  }
  
  x = vegan::mantel(xdis = as.dist(temp.SbNTI), ydis = as.dist(temp.VbNTI), permutations = 9999, method = "spearman", na.rm = T) # Running the Mantel test...
  y = protest(X = as.dist(temp.SbNTI), Y = as.dist(temp.VbNTI), permutations = 9999) # Running a permuted Procrustes analyses
  
  mant = rbind(mant, c(x$statistic, x$signif)) # Storing mantel results
  prot = rbind(prot, c(y$ss, y$scale, y$signif)) # Storing the procrustes results
  names = c(names, paste(Wells[i], " bNTI", sep = ""))
  
  temp.SbNTI = melt(temp.SbNTI, na.rm = T)
  temp.VbNTI = melt(temp.VbNTI, na.rm = T)
  
  temp.bNTI = data.frame(SbNTI = temp.SbNTI$value, VbNTI = temp.VbNTI$value)

  # Generating a plot to demonstrate the relationship (using Mantel statistics, however)
  adj = data.frame(xpos = Inf, ypos = -Inf, text = paste("Mantel r: ", round(x$statistic, digits = 3),
                                                         " p-value:", round(x$signif, digits = 5), sep = ""), hjustvar = 1, vjustvar = 0) # Creating mapping information for labels
  print(
    ggplot(data = temp.bNTI, aes(x = SbNTI, y = VbNTI))+
      geom_point(color = cols$dark[i])+
      geom_smooth(colour = "red", linetype = "longdash", size = 0.8, se = F)+
      geom_text(data = adj, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = text))+
      theme_bw()
  )
  
  # rm("x", "y", "temp.SbNTI", "temp.VbNTI")
}



##############################
### Comparing RCBC results ###
##############################

diag(S3.RCBC) = NA
diag(Vir.U3RC) = NA
diag(Vir.U2RC) = NA

for(i in 1:length(Wells)){
  temp.SRCBC = S3.RCBC[grep(Wells[i], S3.factors$Well), grep(Wells[i], S3.factors$Well)]
  
  temp.SbNTI = S3.bNTI[grep(Wells[i], S3.factors$Well), grep(Wells[i], S3.factors$Well)]
  temp.VbNTI = Vir.bNTI[grep(Wells[i], Vir.factors$Alt_Sample), grep(Wells[i], Vir.factors$Alt_Sample)]
  
  
  if(Wells[i] == "Utica-3"){
    temp.VRCBC = Vir.U3RC
    
    w = c(grep("FT|RT", row.names(temp.VRCBC), invert = T), grep("FT|RT", row.names(temp.VRCBC)))
    temp.VRCBC = temp.VRCBC[w, w]
    
    w = c(grep("FT|RT", row.names(temp.VbNTI), invert = T), grep("FT|RT", row.names(temp.VbNTI)))
    temp.VbNTI = temp.VbNTI[w, w]
  }
  
  if(Wells[i] == "Utica-2"){
    temp.VRCBC = Vir.U2RC
  }
  
  if(any(apply(temp.VRCBC, 1, function(x) length(which(is.na(x) %in% TRUE))) == length(temp.VRCBC[1,]))){
    w = which(apply(temp.VRCBC, 1, function(x) length(which(is.na(x) %in% TRUE))) == length(temp.VRCBC[1,]))
    
    temp.SRCBC = temp.SRCBC[-w,-w]
    temp.SbNTI = temp.SbNTI[-w,-w]
    
    temp.VRCBC = temp.VRCBC[-w,-w]
    temp.VbNTI = temp.VbNTI[-w,-w]
    
  }
  
  x = vegan::mantel(xdis = as.dist(temp.SRCBC), ydis = as.dist(temp.VRCBC), permutations = 9999, method = "spearman", na.rm = T) # Running the Mantel test...
  y = protest(X = as.dist(temp.SRCBC), Y = as.dist(temp.VRCBC), permutations = 9999) # Running a permuted Procrustes analyses
  
  mant = rbind(mant, c(x$statistic, x$signif)) # Storing mantel results
  prot = rbind(prot, c(y$ss, y$scale, y$signif)) # Storing the procrustes results
  names = c(names, paste(Wells[i], " RCBC", sep = ""))
  
  # Looking at only comparisons with insignificant bNTI values
  bak.S = temp.SRCBC
  bak.V = temp.VRCBC
  
  temp.VRCBC[temp.VbNTI > 2] = NA
  temp.SRCBC[temp.SbNTI > 2] = NA
  
  temp.VRCBC[is.na(temp.SRCBC)] = NA
  temp.SRCBC[is.na(temp.VRCBC)] = NA
  
  z = vegan::mantel(xdis = as.dist(temp.SRCBC), ydis = as.dist(temp.VRCBC), permutations = 9999, method = "spearman", na.rm = T) # Running the Mantel test...
  mant = rbind(mant, c(z$statistic, z$signif)) # Storing mantel results
  names = c(names, paste(Wells[i], " important RCBC", sep = ""))
  
  temp.SRCBC = bak.S
  temp.VRCBC = bak.V
  
  # Setting up plots
  temp.SRCBC = melt(temp.SRCBC, na.rm = T)
  temp.VRCBC = melt(temp.VRCBC, na.rm = T)
  
  temp.SbNTI = melt(temp.SbNTI, na.rm = T)
  temp.VbNTI = melt(temp.VbNTI, na.rm = T)
  
  temp.RCBC = data.frame(SRCBC = temp.SRCBC$value, VRCBC = temp.VRCBC$value, 
                         SbNTI = temp.SbNTI$value, VbNTI = temp.VbNTI$value)
  
  
  ##############
  ### Visualizing the complete dataset
  
  int.slo = coef(lm(temp.VRCBC$value~temp.SRCBC$value)) # Storing intercept/slope
  adj = data.frame(xpos = Inf, ypos = -Inf, text = paste("Mantel r: ", round(x$statistic, digits = 3),
                                                         " p-value:", round(x$signif, digits = 5), sep = ""), hjustvar = 1, vjustvar = 0)
  
  temp.RCBC$sig = as.character(cols$light[i])
  temp.RCBC$sig[which(abs(temp.RCBC$SbNTI) < 2 & abs(temp.RCBC$VbNTI) < 2)] = as.character(cols$dark[i])
  
  # Top histogram
  hist_top = ggplot()+geom_point(data = temp.RCBC, aes(x = SRCBC, y = SbNTI), color = temp.RCBC$sig)+
    geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
    theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank())
  
  # I'm not sure, guessing it deals with graph spacing
  empty = ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())
  
  # Actual scatter plot
  scatter = ggplot()+geom_point(data = temp.RCBC, aes(x = SRCBC, y = VRCBC), color = temp.RCBC$sig)+
    geom_abline(intercept = int.slo[1], slope = int.slo[2], colour = "blue", size = 1)+
    geom_smooth(data = temp.RCBC, aes(x = SRCBC, y = VRCBC), colour = "red", linetype = "longdash", size = 0.8, se = F)+
    geom_text(data = adj, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = text))+
    theme_bw()
  
  # Right histogram
  hist_right = ggplot()+geom_point(data = temp.RCBC, aes(x = VRCBC, y = VbNTI), color = temp.RCBC$sig)+
    geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
    theme_bw()+
    theme(legend.position = "none", axis.title.y = element_blank())+
    coord_flip()
  
  # Plotting the graphs
  grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  
  
  #########
  ### Plotting only significant values
  
  temp.RCBC = temp.RCBC[which(abs(temp.RCBC$SbNTI) < 2 & abs(temp.RCBC$VbNTI) < 2),]
  
  # Top histogram
  hist_top = ggplot()+geom_point(data = temp.RCBC, aes(x = SRCBC, y = SbNTI))+
    geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_blank())
  
  # I'm not sure, guessing it deals with graph spacing
  empty = ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(), 
          panel.background=element_blank(), 
          axis.text.x=element_blank(), axis.text.y=element_blank(),           
          axis.title.x=element_blank(), axis.title.y=element_blank())
  
  # Actual scatter plot
  scatter = ggplot()+geom_point(data = temp.RCBC, aes(x = SRCBC, y = VRCBC), color = cols$dark[i])+
    geom_abline(intercept = int.slo[1], slope = int.slo[2], colour = "red")+
    theme_bw()
  
  # Right histogram
  hist_right = ggplot()+geom_point(data = temp.RCBC, aes(x = VRCBC, y = VbNTI))+
    geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
    theme_bw()+
    theme(legend.position = "none", axis.text.y = element_blank())+
    coord_flip()
  
  # Plotting the graphs
  grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  
  rm("x", "y")
}

colnames(mant) = c("Mantel r", "Significance")
row.names(mant) = names
