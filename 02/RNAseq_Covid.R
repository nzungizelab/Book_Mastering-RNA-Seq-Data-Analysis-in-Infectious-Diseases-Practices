#Create a new R script, open Rstudio then click File -> New File -> R script and save it as RnaseqCovid.R

# 1. Set your R working directory to the kallisto_output by locating your project folder into your local computer
#click Session -> Set Working Directory -> Choose Directory... then locate your working directory in your local computer. click on the folder immediately new line will appear in console (copy the line started by "setwd(your directory path)"
setwd("D:/RNAseq_Covid/")

# 2. Installation of sleuth package
# if it isn't installed
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
install.packages("devtools") 
library("devtools")
devtools::install_github('pachterlab/sleuth', ref = 'devel')

# 3. Loading required R libraries

#Load Sleuth and other packages 
#suppress warning messages when loading a library
suppressMessages({library("sleuth")})

suppressMessages({library("dplyr")})
suppressMessages({library("rlang")})
suppressMessages({library("tidyverse")})
suppressMessages({library("ggplot2")})
suppressMessages({library("pheatmap")})

#4. File and data management
outpath <- "D:/RNAseq_Covid/"  #location to save the outputs

#Locate your kallisto output files and if ne slash "/" doesn't work use two slashes "//" 
#Get the sample id
sample_id <- dir(file.path("D:/RNAseq_Covid/kallisto_output"))

#Check the sample id names
sample_id

#Add the specific file path to each files
kal_dirs <- file.path("kallisto_output",sample_id) 

kal_dirs <- sapply(sample_id, function(id) file.path("kallisto_output", id))

#Save kal_dirs into data frame
#The column sample should be in the same order as the associated path and also remember that the column path is required. 
kal_dirs1 <- data.frame(sample = names(kal_dirs), path = kal_dirs)

#Check if the path for each file is added
kal_dirs1

# 5. Affiliate the path to the experiment design (expdesign.txt)
#Note: #open notepad and type your experiment design , condition, treated or untreated (control) and replicates then save as .txt
       #remember to check if "s2c" has five obs. and 2 variables if not open excel sheet and add your experiment design and save as Text (Tab delimited) 
       #if Rstudio doesn't recognize the expdesign.txt at the end create an empty line by pressed enter at the last row of inside expdesign.txt 
       
s2c <- read.table(file.path("D:/RNAseq_Covid/expdesign.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")

#Check if the experiment design table is added correctly 
s2c 

#Attach file paths directories to our data to the design experiment table
# Add file paths to s2c table
# Note: this step will adding each of the files' filepath to the table as a separate column called "path" and turns expdesign.txt table into an object
# Note: make sure that the number of file paths equal to the number of observation in the expdesign.txt table
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# check the table if kallisto output is correctly matched with the samples in experiment design table.
print(s2c)


# 6. Importing gene names into transcripts analysis

# Install cowplot and BioMart
#cowplot package support the sleuth ability to build visuals
#install.packages("cowplot") # if it isn't installed
library(cowplot)

#Install BioMart
#Biomart is an online database which allows the extraction of useful datasets for quantification and data analysis.

#Download and install biomaRt 
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")  # if it isn't installed
biocLite("biomaRt")

#Load biomaRt library
library(biomaRt)

# 7. Choose mart to use
#Check the list mart available and choose which mart we are interested in our project
listMarts()


#We want to use mart "ENSEMBL_MART_ENSEMBL in our project.
ensembl = useMart("ENSEMBL_MART_ENSEMBL") #this line tell BioMart to connect to this specific mart

#In each specific mart under BioMart there are different species by each dataset. 
#The list of datasets with descriptions and version.
listDatasets(ensembl)


#Load the dataset of choice and collect gene names from homo sapiens dataset
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")


# 8. Extract the targeted transcripts

#Extract data from the hsapiens_gene_ensembl dataset via BioMart and get a list of transcripts.
#We have the access to the genomic data for hsapiens_gene_ensembl provided by BioMart.
#add genes names into the table
t2g <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", 
                 "transcript_version", 
                 "ensembl_gene_id",
                 "external_gene_name",
                 "description",
                 "transcript_biotype"), mart = mart)

#Renaming table for sleuth analysis
#Combine two columns (ensembl_transcript_id and transcript_version) both separate by "." into one column called "target_id
t2g <- dplyr::mutate(t2g, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))

#t2g <- dplyr::mutate(t2g, target_id = paste(ensembl_transcript_id)

#Rename the table into three column such as target_id, ens_gene and Gene_Name.
t2g <- dplyr::select(t2g, target_id, ens_gene = ensembl_gene_id, Gene_Name = external_gene_name)

#Check the table "t2g" contains 'ens_gene'(Ensembl gene names) and target_id (associated transcripts from Ensembl).
head(t2g)



 
# 9. Creation of Sleuth object for transcript-level analysis

#Object will help us to get all the information and data necessary to perform sleuth analysis to our kallisto output data.
#This process might take a few minutes
# Create the Sleuth object (a data structure that contains our results). 

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

#so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE, target_mapping = t2g)
#reading in kallisto results
#dropping unused factor levels
#... 
#normalizing est_counts
#74208 targets passed the filter
#normalizing tpm
#merging in metadata
#summarizing bootstraps



#so <- sleuth_prep(s2c,target_mapping = t2g, aggregation_column = 'ens_gene',gene_mode = T, extra_bootstrap_summary = TRUE)



#Note: if the warning message appear it: "Warning message:  In check_num_cores(num_cores):"



#reading in kallisto results
#dropping unused factor levels
#...
#normalizing est_counts
#74208 targets passed the filter
#normalizing tpm
#merging in metadata
#summarizing bootstraps


#so <- sleuth_prep(s2c, target_mapping = t2g,aggregation_column = 'ens_gene')






# 10. Visualizations of the expression data using PCA 

#PCA can help us to check the accuracy of our data
png(paste0(outpath,"PCA_Visualization.png"))
plot_pca(so, color_by = 'condition') + 
  theme_bw()
dev.off()

png(paste0(outpath,"PCA Visualization11.png"))
plot_pca(so, color_by = "condition", text_labels = TRUE, units = "est_counts") +
  theme_bw()
dev.off()

png(paste0(outpath,"PCA Visualization1.png", width=7, height=7, units = "in", res = 300))
plot_pc_variance(so, use_filtered = TRUE, units = "est_counts") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

png(paste0(outpath,"PCA Visualization5.png"))
plot_loadings(so) + theme_bw()
dev.off()

#Plot sample heatmap using the Jensen-Shannon divergence
png(paste0(outpath,"PCA Visualization7.png"))
plot_sample_heatmap(so)
dev.off()

png(paste0(outpath,"PCA Visualization17.png"))
plot_sample_heatmap(so, use_filtered = FALSE)
dev.off()



#adjust the size of point
#png(paste0(outpath,"PCA_data.png"))
#PCA <- plot_pca(so, color_by = 'condition', point_size = 2) +
#  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10)) + 
#  theme(axis.text.y = element_text(face = "bold", color = "black", size = 10)) + 
#  theme(axis.title.x = element_text(face = "bold", color = "black", size = 10)) + 
#  theme(axis.title.y = element_text(face = "bold", color = "black", size = 16))
#dev.off()


# Plot the density of a some grouping
png(paste0(outpath,"Plot density.png"))
plot_group_density(so, use_filtered = TRUE, units = "est_counts",trans = "log", grouping = "condition", offset = 1) +
  theme_bw()
dev.off()

#If we use argument "use_filtered = TRUE", can help us to look the distributions of the filtered genes used for differential expression analysis.

png(paste0(outpath,"Plot density1.png"))
plot_group_density(so, use_filtered = FALSE, units = "est_counts",trans = "log", grouping = "condition") + 
  theme_bw()
dev.off()

#The count distributions on the y-axis represent the proportion of genes associated with the number of counts on the x-axis


  

# 11. Creation of data models

#Create two model (full model and condition model)
#Model will used each other to calculate the gene expression levels in certain tests
#With the model of our data we can start the differential expression analysis of individual genes


#Create a full model that contains all covariates 
#After full model ,then create the additional model where both models will be compared to identify differential expression 

#fit the full model
#so <- sleuth_fit(so, ~ condition) 
so <- sleuth_fit(so, ~condition, 'full')

#Alternative command line
so <- sleuth_fit(so, ~condition, fit_name = "full")

#Note:#You cannot perform differential expression analysis under sleuth package without replicates.
      #sleuth_fit is treating both samples (treated and untreated) as replicates

#fit the reduced model
so <- sleuth_fit(so, ~1, 'reduced')

#Alternative command line
so <-  sleuth_fit(so, ~1, fit_name = "reduced")

#fitting measurement error models
#shrinkage estimation
#computing variance of betas

#Finally perform the test by comparing both models
so <- sleuth_lrt(so, 'reduced', 'full')

#Checking the models created by which have been fit use model function
models(so)


#[  full  ]
#formula:  ~condition 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#	(Intercept)
# 	conditionpatient
#[  reduced  ]
#formula:  ~1 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#  (Intercept)

#Alternative command line
so <- sleuth_lrt(so, null_model = "reduced", alt_model = "full")


#check the test
tests(so)

#~likelihood ratio tests:
#  reduced:full

#~wald tests:
#  no tests found.


# extract the results of the likelihood
lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(lrt_results[,"qval"] < 0.05)

#Output
#FALSE  TRUE 
#73740    10 
#In total of 73740  transcripts only 10 have significantly expression pattern.


#lrt_results %>% head(n = 20) %>% dplyr::select(target_id, qval, ensembl_gene_id, ext_gene)




#Generate the results table for analysis 

#summarize the sleuth results and view 10 most significant DE transcripts with q-value <= 0.05.
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 10)

dim(sleuth_table)

#[1] 73750    12 

dim(sleuth_significant)
#[1] 10 12

# 10 transcripts out of 73750 was DE

#Save the DTE into .csv table
write.csv(sleuth_significant,"sleuth_DE_Covid19.csv",row.names=FALSE,quote=FALSE)





#12. Plot the normalized bootstraps across all samples

#12.1. Plot variation in transcripts per est_counts


# Using the list of DE transcript in the table created you can plot variation among samples.
#We selected among top 10 genes DE
png(paste0(outpath,"Plot.png"))
plot_bootstrap(so, "NM_001193370.2", units = "est_counts", color_by = "condition") +
  theme_bw()
dev.off() 


#Alternative option to save it in single pdf file
pdf(paste0(outpath,"Plot.pdf"))
plot_bootstrap(so, "NM_001193370.2", units = "est_counts", color_by = "condition") + theme_bw()
dev.off()


png(paste0(outpath,"Plot02.png"))
plot_bootstrap(so, "NR_146118.1", units = "est_counts", color_by = "condition") 
dev.off()




#12.2. Plot variation in transcripts per million (TPM)

png(paste0(outpath,"TPM_plot.png"))
plot_bootstrap(so, "NM_001193370.2", units = "tpm", color_by = "condition") +
  theme_bw()
dev.off()

png(paste0(outpath,"TPM_plot1.png"))
plot_bootstrap(so, "NR_146118.1", units = "tpm", color_by = "condition") +
  theme_bw()
dev.off()

#clustered heatmap of these top genes
png(paste0(outpath,"heatmap_top101.png"))
plot_transcript_heatmap(so, head(lrt_results, n = 10)$target_id, 'est_counts')
dev.off()

png(paste0(outpath,"heatmap_top102.png"))
plot_transcript_heatmap(so, head(lrt_results, n = 20)$target_id, 'est_counts')
dev.off()






#Note : for your command line sleuth_prep() if you didn't specified boostraps in tpm units , you cant generated plot with TPM



# 3. perform differential analysis (testing)

#Explore expression variation among samples using shiny.
#Sleuth live function gives you a web application for interactive visualization of differential expression analysis.
library(Rcpp)
sleuth_live(so)

sessionInfo()
