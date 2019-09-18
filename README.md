# Studying the ecological assembly of shale microbial/viral communities
The set of scripts used to process data with in the paper entitled "Ecological assembly processes are coordinated between bacterial and viral communities in fractured shale ecosystems"

Each script contains a distinct analysis:
- <b>Correlating_S3-Viral_Models.R</b> will allow you to relate the null models used throughout this study according to a Spearman-based Mantel test. It will also generate figures.
- <b>Process_Viral_Eco_Models.R</b> generates viral bNTI and RCBC boxplots
- <b>S3_Community_bNTI.R</b> is the bNTI implementation for the microbial community
- <b>S3_Community_RCBC.R</b> is the RCBC implementation for the microbial community
- <b>S3_Phylogeny_Plot.R</b> generates a stacked-bar chart for each of the fractured shale wells studied, as well as a plot illustrating the proportion represented by <i>Halanaerobium</i> through time
- <b>vContact_Processing.R</b> converts vConTACT viral similarity data into a UPGMA dendrogram
- <b>Viral_Community_bNTI.R</b> is the bNTI implementation for the viral community
- <b>Viral_Community_RCBC.R</b> is the RCBC implementation for the viral community

All necessary files should be located in the "Date Files" folder. 
If you have any questions, please contact me at danczak.6@gmail.com
