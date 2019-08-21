# Script to perform a Raup-Crick(Bray-Curtis) analysis on viral data
# This code prevents cross-well calculations becuase they take a prohibitively long time if utilized with the raw counts configuration

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# RED 2019

local = T # Whether on the server or not
U2 = T # Flag to remove all data aside from U_2
U3 = F # Flag to remove all data aside from U_3
FINAL = T # Relative abudances do not play well with Raup-Crick. This flag will help alleviate these problems.
match.tree = T # Ensures only viruses present on the tree are analyzed

library(vegan)
library(picante)
library(doMC)

###################################
#### Data Loading and cleaning ####
###################################
if(local == T){
  setwd("~/Documents/Shale Ecological Modeling/Documents/Supplemental/Data with only U2-U3/") # Working directory for 16S comm on my local computer
  data = read.csv("Viral_Counts_Normalized.csv", row.names = 1) # Importing the organismal data
  tree = read.tree("vContact_UPGMA_Normalized.tre") # Importing the tree
  
  registerDoMC(cores = 1) # Setting the number of processors
  rc.reps = 9999
} else {
  setwd("//path/to/directory/") # Working directory for 16S comm on my local computer
  data = read.csv("Viral_Counts_Normalized.csv", row.names = 1) # Importing the organismal data
  tree = read.tree("vContact_UPGMA_Normalized.tre") # Importing the tree
  
  registerDoMC(cores = 20) # Setting the number of processors
  rc.reps = 9999
}


##############################
#### Cleaning up the data ####
##############################

if(U2 == T){
  data = data[,grep("U2_", colnames(data))]
  data = data[-which(rowSums(data) == 0),]
}

if(U3 == T){
  data = data[,grep("U3_", colnames(data))]
  data = data[-which(rowSums(data) == 0),]
}


# ######################## #
#### Matching bNTI Tree ####
# ######################## #

# Added January 28th, 2019
if(match.tree == T){
  # Matching the provided tree
  data = match.phylo.data(phy = tree, data = data)
  data = data$data
  
  # Removing absent viruses
  if(length(which(rowSums(data) == 0)) > 0){
    data = data[-which(rowSums(data) == 0),]
  }
} # Ensures only viruses on the tree are used in RCBC analyses (like-for-like community comparisons)

############################
#### Rarefying the data ####
############################

if(!(FINAL == T)){
  data = t(rrarefy(t(data), 50000))
}

############################
#### Running Raup-Crick ####
############################

# Setting up the data to mesh with the Stegen et al. code
spXsite = t(data)

if(FINAL == T){
  spXsite[spXsite>0] = spXsite[spXsite>0]*1000
} # Adding one to my RPKM value to generate an "FINAL+1" stat. This is necessary to work around Raup-Crick

# Count number of sites and total species richness across all plots (gamma)
n_sites = nrow(spXsite)
gamma = ncol(spXsite)

# Build a site by site matrix for the results, with the names of the sites in the row and col names:
results = matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))

# Make the spXsite matrix into a new, pres/abs. matrix:
spXsite.inc = ceiling(spXsite/max(spXsite))

# Create an occurrence vector- used to give more weight to widely distributed species in the null model
occur = apply(spXsite.inc, MARGIN=2, FUN=sum)

# Create an abundance vector- used to give more weight to abundant species in the second step of the null model
abundance = apply(spXsite, MARGIN=2, FUN=sum)

# Loops through every pairwise comparison, generating null results

for(null.one in 1:(nrow(spXsite)-1)){
  for(null.two in (null.one+1):nrow(spXsite)){
    
    null_bray_curtis<-NULL
    null_bray_curtis = foreach(i=1:rc.reps, .packages = c("vegan","picante")) %dopar% {
      
      # Generates two empty communities of size gamma
      com1<-rep(0,gamma)
      com2<-rep(0,gamma)
      
      # Add observed number of species to com1, weighting by species occurrence frequencies
      com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
      com1.samp.sp = sample(which(com1>0), (sum(spXsite[null.one,])-sum(com1)), replace=TRUE, prob=abundance[which(com1>0)]);
      com1.samp.sp = cbind(com1.samp.sp,1);
      com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2], com1.samp.sp[,1], FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
      com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts));
      com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts;
      #sum(com1) - sum(spXsite[null.one,]); # This should be zero if everything worked properly
      rm('com1.samp.sp','com1.sp.counts');			
      
      # Again for com2
      com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
      com2.samp.sp = sample(which(com2>0), (sum(spXsite[null.two,])-sum(com2)), replace=TRUE, prob=abundance[which(com2>0)]);
      com2.samp.sp = cbind(com2.samp.sp,1);
      com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2], com2.samp.sp[,1], FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
      com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts));
      com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts;
      # sum(com2) - sum(spXsite[null.two,]); # This should be zero if everything worked properly
      rm('com2.samp.sp','com2.sp.counts');
      
      null.spXsite = rbind(com1,com2); # Null.spXsite
      
      # Calculates the null Bray-Curtis
      null_bray_curtis[i] = vegdist(null.spXsite, method='bray');
      
    }; # End of the null loop
    
    # Unlisting the parallel list
    null_bray_curtis = unlist(null_bray_curtis)
    
    # Calculates the observed Bray-Curtis
    obs.bray = vegdist(spXsite[c(null.one,null.two),], method='bray');
    
    # How many null observations is the observed value tied with?
    num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
    
    # How many null values are smaller than the observed *dissimilarity*?
    num_less_than_in_null = sum(null_bray_curtis<obs.bray);
    
    rc = ((num_less_than_in_null + (num_exact_matching_in_null)/2)/rc.reps) # This variation of rc splits ties
    
    rc = (rc-.5)*2 # Adjusts the range of the  Raup-Crick caclulation to -1 to 1
    
    results[null.two,null.one] = round(rc,digits=2); # Stores rc into the results matrix
    
    print(c(null.one,null.two,date())); # Keeps track of position
    
  }; # End of inner loop
  
}; # End of outer loop

rc.results = as.dist(results) # Converts results into a distance matrix

if(U_2 == T){
  write.csv(as.matrix(rc.results), "U2_RCBC_NormCounts_NormUPGMA.csv", quote = F)
}

if(U_3 == T){
  write.csv(as.matrix(rc.results), "U3_RCBC_NormCounts_NormUPGMA.csv", quote = F)
}

rm('spXsite.inc')
