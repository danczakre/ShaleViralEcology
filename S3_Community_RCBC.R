# Modified RCBC script for the S3 microbial community

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# RED 2019


local = F # Whether on the server or not
rpkm = F # So, relative abudances do not play well with Raup-Crick. This loop will help alleviate these problems.

library(vegan)
library(picante)
library(doMC)

rc.reps = 9999

###################################
#### Data Loading and cleaning ####
###################################
if(local == T){
  setwd("//path/to/community_data/") # Working directory for 16S comm on my local computer
  data = read.csv("S3_Community_TotalCounts.csv", row.names = 1) # Importing the organismal data
  tree = read.tree("S3_Community_trimmed.RAxML.tre") # Importing the tree
  registerDoMC(cores = 1) # Setting the number of processors
} else {
  setwd("//path/to/community_data/") # Working directory for 16S comm on my local computer
  data = read.csv("S3_Community_TotalCounts.csv", row.names = 1) # Importing the organismal data
  tree = read.tree("S3_Community_trimmed.RAxML.tre") # Importing the tree
  registerDoMC(cores = 20) # Setting the number of processors
}

############################
#### Running Raup-Crick ####
############################

# Setting up the data to mesh with the Stegen et al. code
spXsite = t(data)

if(rpkm == T){
  spXsite[spXsite>0] = spXsite[spXsite>0]+1
} # Adding one to my RPKM value to generate an "RPKM+1" stat. This is necessary to work around Raup-Crick

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

write.csv(as.matrix(rc.results), "S3_RCBC_TotalCounts.csv", quote = F)

rm('spXsite.inc')