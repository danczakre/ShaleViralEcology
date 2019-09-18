### Plotting organisms identified by rpS3 across time

# Cleaning up the work environment
if(length(dev.list()) > 1){
  dev.off()
} # Clears plots - if loop needed if there are no plots present
rm(list=ls()) # Clears the environment

library(ggplot2); library(ggthemes) # Plotting
library(reshape2) # Melting data for ggplot2

taxa.level = "Family" # Specifies which level of taxonomy will be plotted
Wells = c("Utica-3", "Utica-2") # Specifies which wells to plot (should only be U3 and U2)

#######################
### Loading in data ###
#######################

# Setting work directory and loading in data
setwd("//path/to/directory/")

reads = read.csv("S3_Community_RelAbund.csv", row.names = 1)
taxa = read.csv("S3_Taxonomy.csv")
factors = read.csv("S3_Community_Factors_Sheet.csv")

### Cleaning up data ###
reads = reads[,-grep("HT|LW", factors$Sample_Name)] # Removing the HT and LW samples as they are no longer needed
factors = factors[-grep("HT|LW", factors$Sample_Name),]

##################################
### Assigning taxonomy to data ###
##################################

reads.taxa = reads # Creating a new object for modification purposes

reads.taxa$taxa = NA
level.loc = grep(taxa.level, colnames(taxa))

for(f in 1:length(reads.taxa$taxa)){
  reads.taxa$taxa[f] = as.character(taxa[grep(row.names(reads.taxa)[f], taxa$Sequence), level.loc])
}

reads.taxa = aggregate(. ~ taxa, data = reads.taxa, FUN = "sum")
row.names(reads.taxa) = reads.taxa[,1]
reads.taxa = reads.taxa[,-1]

#############################
### Generating some plots ###
#############################

# Adding some information to my data frame to make melting more efficient
reads.taxa = as.data.frame(t(reads.taxa)) # Transposing for simplicity
temp = reads.taxa[,which(apply(reads.taxa, 2, max) < 0.05)] # Storing community members that never go above a certain threshold
reads.taxa = reads.taxa[,-which(apply(reads.taxa, 2, max) < 0.05)] # Removing said members
reads.taxa$Other = rowSums(temp) # Adding the abundances of those organisms back into the data frame

rm("temp") # Removing temporary storage

# Adding some information to the data frame to make melting more informative
reads.taxa$Day = factors$Day
reads.taxa$Well = factors$Well

reads.taxa = melt(data = reads.taxa, id.vars = c("Day", "Well")) # Melting the data for freindly play with ggplot2
reads.taxa$value = reads.taxa$value*100

# Plotting wells through time
for(f in 1:length(Wells)){
  
  temp = reads.taxa[grep(Wells[f], reads.taxa$Well),]
  
  # Need to turn days into factors to allow bar chart plotting
  if(Wells[f] == "Utica-3"){
    temp$Day = factor(as.character(temp$Day), levels = c("-25", "-24", "38", "46", "54", 
                                                         "60", "96", "124", "159"))}
  
  if(Wells[f] == "Utica-2"){
    temp$Day = factor(as.character(temp$Day), levels = c("-5", "0", "86", "87", "88", 
                                                         "92", "94", "96", "98", "105", 
                                                         "112", "126", "140", "154", 
                                                         "175", "204", "302"))}
  
  # Actually plotting...
  print(
    ggplot(data = temp, aes(x = Day, y = value, fill = variable))+
      geom_bar(stat = "identity")+
      ylab("Relative Abundance (%)")+
      ggtitle(Wells[f])+
      theme_bw()+
      theme(text = element_text(size = 14),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 1),
            axis.text.y = element_text(colour = "black"),
            panel.border = element_rect(size = 1, colour = "black"),
            panel.grid = element_blank())
  )
}

hala = reads.taxa[grep("Halanaero", reads.taxa$variable),]
hala$Day = factor(as.character(hala$Day), levels = c("-25", "-24", "-5", "0", "38", "46", "54", 
                                                     "60", "86", "87", "88", 
                                                     "92", "94", "96", "98", "105", 
                                                     "112", "124", "126", "140", "154", "159",
                                                     "175", "204", "302"))
hala$Well = gsub("Wagner", "Utica-2", hala$Well)
hala$Well = gsub("Shriver", "Utica-3", hala$Well)

print(
  ggplot(data = hala, aes(x = Day, y = value))+
    geom_bar(stat = "identity", aes(fill = Well))+
    ylab("Proportion Halanaerobiaceae (%)")+
    scale_fill_manual(values = c("darkgreen", "dodgerblue3"))+
    facet_grid(.~Well, scales = "free_x")+
    theme_bw()+
    theme(text = element_text(size = 14),
          axis.text.x = element_text(colour = "black", angle = 37.5, hjust = 1, vjust = 1),
          axis.text.y = element_text(colour = "black"),
          panel.border = element_rect(size = 1, colour = "black"),
          panel.grid = element_blank())
)

rm("f", "temp", "level.loc")