# Modified bNTI script for the Shale S3 (and probably other) datasets

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# RED 2019; danczak.6@osu.edu

library(vegan)
library(picante)

###################################
#### Data Loading and cleaning ####
###################################

setwd("//path/to/microbial_community_data") # Working directory for 16S comm on my local computer
data = read.csv("S3_Community_TotalCounts.csv", row.names = 1) # Importing the organismal data
tree = read.tree("rps3_trimmed.RAxML.tre") # Importing the tree


####################################
#### Beginning the bNTI Process ####
####################################

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree, data)

# Calculating bMNTD for my samples
bMNTD = as.matrix(comdistnt(t(phylo$data), cophenetic(phylo$phy), abundance.weighted = T))

# Calculating the bMNTD for 999 random distributions
bMNTD.rand = array(c(-999), dim = c(ncol(phylo$data), ncol(phylo$data), 999)) # Creating 999 'dummy' matrices to put random results into

for(i in 1:999){
  bMNTD.rand[,,i] = as.matrix(comdistnt(t(phylo$data), taxaShuffle(cophenetic(phylo$phy)), abundance.weighted = T, exclude.conspecifics = F))
  print(c(date(),i)) # Measuring how far the loop has gotten
} # Performing the calculations on using the OTU table but with randomized taxonomic affiliations

# Calculating the bNTI between these randomized communities and the empirical dataset
bNTI = matrix(c(NA), nrow = ncol(phylo$data), ncol = ncol(phylo$data))

for(i in 1:(ncol(phylo$data)-1)){ 
  for(j in (i+1):ncol(phylo$data)){
    m = bMNTD.rand[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

rownames(bNTI) = colnames(phylo$data)
colnames(bNTI) = colnames(phylo$data)

write.csv(bNTI, "S3_Community_bNTI_TotalCounts.csv", quote = F)