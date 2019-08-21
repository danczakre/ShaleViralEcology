### Script to perform viral bNTI
# RED 2019

local = T # Flag if it was run on server or not
norm = T # Flag whether norm counts were used
raw = F # Flag if raw counts were used

library(vegan)
library(picante)

###################################
#### Data Loading and cleaning ####
###################################
if(local == T){
  setwd("//path/to/directory/") # Working directory for 16S comm on my local computer
  data = read.csv("Viral_Normalized_Counts.csv", row.names = 1) # Importing the organismal data
  tree = read.tree("vContact_UPGMA_Normalized.tre") # Importing the tree
} else {
  setwd("//path/to/directory/") # Working directory for 16S comm on my local computer
  data = read.csv("Viral_Normalized_Counts.csv", row.names = 1) # Importing the organismal data
  tree = read.tree("vContact_UPGMA_Normalized.tre") # Importing the tree
}


###################################
#### Running the bNTI pipeline ####
###################################

# If the data comes in the form of raw counts, I need to rarefy it to match the RCBC analysis
if(raw == T){
  data = t(rrarefy(t(data), 50000))
}

# Given that the normalized data is pseudo-relative abundance,it needs to be multiplied by a
# common factor to allow RCBC to work (current formulation doesn't like non-integers)
# I'm working on a different way to accept rel. abund. data.
if(norm == T){
  data[data>0] = data[data>0]*10000
}

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree, data)

# Calculating bMNTD for my samples
bMNTD = as.matrix(comdistnt(t(phylo$data), cophenetic(phylo$phy), abundance.weighted = T))

# Calculating the bMNTD for 999 random distributions
bMNTD.rand = array(c(-999), dim = c(ncol(phylo$data), ncol(phylo$data), 999)) # Creating 999 'dummy' matrices to put random results into

for(i in 1:999){
  bMNTD.rand[,,i] = as.matrix(comdistnt(t(phylo$data), taxaShuffle(cophenetic(phylo$phy)), abundance.weighted = T))
  
  if(length(which(bMNTD.rand[,,i] == 0)) > 1){
    w = which(bMNTD.rand[,,i] == 0) # Finding zeroes in random dataset
    x = sample(x = seq(0.0000001, 0.0000002,by = 0.0000001), size = length(w), replace = T)
    bMNTD.rand[,,i][w] = ((bMNTD.rand[,,i][w])+x)
    diag(bMNTD.rand[,,i]) = 0
    
    rm("w", "x")
  } # Used to add more noise to the data
  
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

bNTI[upper.tri(bNTI)] = t(bNTI)[upper.tri(bNTI)]
boxplot(bNTI[24:40, 24:40])
abline(h = c(-2,2))

bNTI.melt = reshape2::melt(bNTI[24:40, 24:40])

ggplot2::ggplot(data = bNTI.melt, ggplot2::aes(x = Var1, y = Var2, fill = value))+
  ggplot2::geom_tile()+
  ggplot2::scale_fill_gradient2(low = "darkblue", mid = "gray80", high = "darkred")+
  ggplot2::theme_bw()

write.csv(bNTI, "Viral_bNTI_NormCounts_NormTree.csv", quote = F)
