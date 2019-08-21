### Converting vCONtact data into a relational tree
# RED 2019

library(reshape2) # Needed to melting/recasting
library(ape) # Necessary for some tree-based functions
library(phangorn) # Necessary for some tree-based functions

First.Time = F # Flag to create a list of unique samples if this is the first time running.
               # This was necesssary because underscores were dropped at some point and 
               # they are necessary to match the tree to viral data
perform.norm = T # This flag will perform a range normalization on the vContact similarities

setwd("//path/to/vContact_data")

# Load in vContact file
data = read.csv(file = "vContact_Data.csv", header = 1)

# Remove reference viral genomes
data = data[-grep("NC_", data$genome_1),]
data = data[-grep("NC_", data$genome_2),]

# Renaming viral genomes so that they match the counts dataset
# Need to create a list of unique identifiers (this only needs to be run once)
uniq = unique(data$genome_1, data$genome_2)
uniq = gsub(pattern = "aff.*$", replacement = "", x = uniq)
uniq = unique(uniq)
uniq = uniq[order(uniq)]

if(First.Time == T){
  write.csv(uniq, "vContact_to_Counts.csv", quote = F, row.names = F) # This creates a file that makes renaming easier
  # This file must be editted to match the names of the other data files
  stop("Given that this is a first time run, the script is stopping")
} else {
  names = read.csv("vContact_to_Counts.csv", header = 1)
  
  # Renaming viruses to match the abundance table
  for(i in 1:length(names[,1])){
    data$genome_1 = gsub(pattern = names$vContact_Names[i], replacement = names$Counts_Names[i], data$genome_1)
    data$genome_2 = gsub(pattern = names$vContact_Names[i], replacement = names$Counts_Names[i], data$genome_2)
  }
  
  # Fixing some naming problems not addressed above
  data$genome_1 = gsub(pattern = "trim", replacement = "_trim", data$genome_1)
  data$genome_2 = gsub(pattern = "trim", replacement = "_trim", data$genome_2)
  
  data$genome_1 = gsub(pattern = "scaffoldd", replacement = "scaffold", data$genome_1)
  data$genome_2 = gsub(pattern = "scaffoldd", replacement = "scaffold", data$genome_2)
  
}

# Casting data into similarity matrix
data = dcast(data, genome_1~genome_2)

# Setting row names
row.names(data) = data[,1]
data = data[,-1]

# Changing NA values to 0; 0 implicates that viruses are not related
data[is.na(data)] = 0

# Range normalization of dissimilarity
if(perform.norm == T){
  norm = (data-min(data))/(max(data)-min(data)) # Normalizing the similarity to range from 0 to 1
  dis.norm = as.dist(1-norm) # Converting the similarity measurement to dissimilarity
  
  # Creating dendrogram
  av = upgma(dis.norm)
  av = as.phylo(av)
  
  write.tree(av, "vContact_UPGMA_Normalized.tre")
} else {
  # Creating dendrogram
  av = upgma(data)
  av = as.phylo(av)
  
  write.tree(av, "vContact_UPGMA.tre")
}