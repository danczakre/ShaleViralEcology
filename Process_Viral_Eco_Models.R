### Script to plot boxplots and heatmaps for viral bNTI and RCBC results
# RED 2019

bNTI.fix = "U2" # Processes only one well at a time; acceptable flags are "U2" or "U3"

library(ggplot2) # For plotting
library(reshape2) # Formatting data for ggplot

#######################
#### Load in files ####
#######################
setwd("~/Documents/Shale Ecological Modeling/Documents/Supplemental/Data with only U2-U3/")
bNTI = read.csv("Viral_bNTI_NormCounts_NormUPGMA.csv", row.names = 1)
factors = read.csv("Viral_Model_Factors.csv")


############################
#### Some Data Cleaning ####
############################
if(bNTI.fix == "U2"){
  
  # Keeping only data from Utica-2
  w = grep("U2", row.names(bNTI))
  bNTI = bNTI[w, w]
  RCBC = read.csv("U2_RCBC_NormCounts_NormUPGMA.csv",row.names = 1)
  
  w = grep("U2", factors$Sample_Name)
  factors = factors[w,]
  
  RCBC = RCBC[-grep("HT|LW", row.names(RCBC)), -grep("HT|LW", row.names(RCBC))]
  bNTI = bNTI[-grep("HT|LW", row.names(bNTI)), -grep("HT|LW", row.names(bNTI))]
  factors = factors[-grep("HT|LW", factors$Sample_Name),]
}

if(bNTI.fix == "U3"){
  # Keeping only data from Utica-3
  w = grep("U3_", row.names(bNTI))
  bNTI = bNTI[w, w]
  RCBC = read.csv("U3_RCBC_NormCounts_NormUPGMA.csv", row.names = 1)
  
  w = grep("U3_", factors$Sample_Name)
  factors = factors[w,]
  
  if(identical(x = row.names(RCBC), y = row.names(bNTI)) == FALSE){
    w = which(!(row.names(bNTI) %in% row.names(RCBC)) %in% TRUE)
    RCBC = rbind(RCBC[1:(w-1),], NA, RCBC[w:length(RCBC[,1]),])
    RCBC = cbind(RCBC[,1:(w-1)], NA, RCBC[,w:length(RCBC[1,])])
    
    colnames(RCBC)[w] = colnames(bNTI)[w]
    row.names(RCBC)[w] = row.names(bNTI)[w]
    
    if(identical(x = row.names(RCBC), y = row.names(bNTI)) == FALSE){
      stop("Something went wrong when renaming the dimensions")
    }
  }
  
}

rm("w")

# Need to do some processing on the bNTI data
bNTI[upper.tri(bNTI)] = t(bNTI)[upper.tri(bNTI)]


############################
#### Plot RCBC boxplots ####
############################
RCBC[abs(bNTI) > 2] = NA # Replacing RCBC values where the bNTI is significant
diag(RCBC) = NA # Self comparisons as zero could inflate/deflate differences

heat = melt(as.matrix(RCBC))
ggplot(data = heat, aes(x = Var1, y = Var2, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#79539B", mid = "#F9F9F9", high = "#306338")+
  xlab(label = "")+
  ylab(label = "")+
  theme_bw()+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

RCBC$Day = factor(factors$Day)
RCBC.melt = melt(data = RCBC, id.vars = c("Day"))

ggplot(data = RCBC.melt, aes(x = Day, y = value))+
  geom_boxplot()+
  geom_hline(yintercept = c(0.95,-0.95), color = "red", linetype = 2)+
  theme_bw()


############################
#### Plot bNTI boxplots ####
############################
diag(bNTI) = NA # Self comparisons as zero could inflate/deflate differences

heat = melt(as.matrix(bNTI))
ggplot(data = heat, aes(x = Var1, y = Var2, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "#7797BE", mid = "#F9F9F9", high = "#DD4E39")+
  xlab(label = "")+
  ylab(label = "")+
  theme_bw()+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())
  

bNTI$Day = factor(factors$Day)
bNTI.melt = melt(data = bNTI, id.vars = c("Day"))

ggplot(data = bNTI.melt, aes(x = Day, y = value))+
  geom_boxplot()+
  geom_hline(yintercept = c(2,-2), color = "red", linetype = 2)+
  theme_bw()
