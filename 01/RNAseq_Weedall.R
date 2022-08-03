##  Author : Lambert Nzungize
##  Email  : nzulapa@outlook.com
## 
###################################################################


# Set working directory 
setwd("D:/RNAseq_Weedall/7.DE")


# 1. Reading in the Data

# load the required library
library(edgeR)
library(limma)
library(pheatmap)

# 2. Import raw_data

#sep argument says what character is used toseparate items
#The way to tell R R that you mean tab character is to use "\t" (\ mean treat what follow it as special)
#data_raw1 <- read.table("counts.matrix", sep="\t", header = TRUE)

#Note: when using "read.table" you get "data.frame" where the first column is a factor with geneids are levels and the other columns contain integer
#It shows the error with "DGEList" like <The count matrix is a data.frame instead of a matrix and the first column is of class character>
#Solution: use "read.lim" because DGEList is expecting to get matrix rather than "data.frame"

data_raw2 <-read.delim("counts.matrix", row.names = 1, sep="\t", header = TRUE) # counts can be read into a data.frame
expdesign <- read.table("expdesign", header=T)
outpath <- "D:/RNAseq_Weedall/7.DE/"  #location to save the outputs

# Import .gtf file via rtracklayer::
AfunF3_gtf <- rtracklayer::import("Anopheles-funestus-FUMOZ_BASEFEATURES_AfunF3.gtf")
AfunF3_gtf_df = as.data.frame(AfunF3_gtf)   

# Save table in .csv format or in .txt format
write.table (AfunF3_gtf_df, "AfunF3_gtf_df.csv", sep="\t", row.names=F)
write.table (AfunF3_gtf_df, "AfunF3_gtf_df.txt", sep="\t", row.names=F)

ls()
dim(data_raw2)  # Check dimensions
head(data_raw2)
head(expdesign)
head(AfunF3_gtf_df)
summary(data_raw2)


# 3. Quality control 

#check the library size for each sample
colSums(data_raw2[,1:ncol(data_raw2)]) #remember to specify the column "data[,1]"

#check the total detected genes
apply(data_raw2[,2:ncol(data_raw2), drop=F], 2, function(c)sum(c!=0)) #starting from col 3
apply(data_raw2[,1:ncol(data_raw2), drop=F], 2, function(c)sum(c!=0)) 
apply(data_raw2[,1, drop=F], 2, function(c)sum(c!=0)) # only select column 1

# 4.  Creating a DGEList object

# edgeR stores data in a simple list-based data object called a DGEList
# Create a DGEList object to hold our read counts
# generate DGE object
# Merging count reads with experiment design table
#DGEList object has "group" which is the basic experimental group for each sample.
y <- DGEList(data_raw2, samples=expdesign, group=expdesign$condition) #Converting counts to DGEList object

#Alternative option:
groups <- factor(expdesign$condition)
y1 <- DGEList(counts=data_raw2, samples=expdesign, group=groups)

#Explore the dataset
y
y$samples
head(y$counts) #Many rows!

# 4.Filtering

#There are approximately 14176 genes in this dataset, then filtering low-expression genes can improve DEG detection sensitivity.
# Many of them will not be expressed, or will not be represented by enough reads to contribute to the analysis. 
# we retain only those genes that are represented at least 1cpm reads in at least two samples (cpm=counts per million).
#filtering is based CPM values to avoid genes that are expressed in larger libraries over those expressed in smaller libraries.
countsPerMillion <- cpm(y) #Calculate the Counts Per Million measure
summary(countsPerMillion) ##'summary' is a useful function for exploring numeric data
head(countsPerMillion)

countCheck <- countsPerMillion > 1 #Identify genes with at least 1 cpm in at least 2 samples
head(countCheck)
tail(countCheck)
table(countCheck) #Checking before filtering step

keep <- which(rowSums(countCheck) >= 2)
y_filtered <- y[keep,] #keep the more highly expressed genes
dim(y_filtered) #Checking filtered data set
summary(cpm(y_filtered)) #compare this to the original summary

#Average-log-CPM for each gene
AveLogCPM <- aveLogCPM(y)
png(paste0(outpath,"average-log-CPM.png"))
hist(AveLogCPM)  #histogram can help to choose a cutoff value heuristically
dev.off()

# 5. Normalization

#Trimmed mean of M values (TMM) normalization estimates sequencing depth after excluding genes for which the ratio of counts between a pair of experiments is too extreme or for which the average expression is too extreme. 
# edgeR implements the Trimmed Mean of M-values (TMM) method.
#uses TMM method to eliminate RNA composition effect (Robinson and Oshlack 2010)
#Estimate normalization factors using
y_filtered_normalized <- calcNormFactors(y_filtered, method="TMM") #Perform TMM normalisation
y_filtered_normalized$samples


# 6. Data exploration

#Plotting correlation of samples with normalized counts
png(paste0(outpath,"normalized_sample_correlation.png"))
pheatmap(cor(log2(cpm(y_filtered_normalized)+1)), main="Correlation")
dev.off()

#Plotting library sizes 
#Number of reads assigned to all genes from the annotation .gtf file 
png(paste0(outpath,"library_sizes.png"))
barplot(y_filtered_normalized$samples$lib.size,las=2,names=colnames(y_filtered_normalized),
        cex.names=0.6,horiz=T, col=groups, main="Library Sizes")
dev.off()

#We can examine inter-sample relationships by producing a plot based on mutlidimensional scaling.
png(paste0(outpath,"sample_interelationships.png"))
plotMDS(y_filtered_normalized) 
dev.off()

#  MDS plot showing distances between expression profiles
png(paste0(outpath,"sample_interelationships1.png"))
pch <- c(9,3,6,8,13)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(y_filtered_normalized, col=colors[groups]) 
legend("topleft", legend=levels(groups), pch=pch, col=colors, ncol=2)
dev.off()
#The figure shows differences between groups are much larger than those within groups( indicate that there is statistically significant differences between the groups) 

# MD plot visualizes individual samples
# Each point represents a gene, and the red line indicates a log-ratio of zero.
#The majority of points cluster around the red
# It is good practice to make MD plots for all the samples as a quality check.
# MD plot shows many genes that are both highly expressed and highly up-regulated.
# If a number of points are in the upper right of MD plot, meaning that genes are both highly expressed and highly up-regulated
# In our case CMR-1 and GHA-1 are highly expressed only
head(y_filtered_normalized)

# sample CMR-1
png(paste0(outpath,"MDplot_compares_sample1_to_ref_library.png"))
plotMD(y_filtered_normalized, column=1, main = "DEG within sample CMR-1")
abline(h=0, col="red", lty=2, lwd=2)
dev.off()
# The Figure indicates Mean-difference plot of log2-expression in CMR-1 versus the average log2-expression across all other samples

# sample GHA-1
png(paste0(outpath,"MDplot_compares_sample2_to_ref_library.png"))
plotMD(y_filtered_normalized, column=2, main = "DEG within sample GHA-1")
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

# Sample FANG-4
png(paste0(outpath,"MDplot_compares_sample3_to_ref_library.png"))
plotMD(y_filtered_normalized, column=3, main = "DEG within sample FANG-4")
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

#View an MDS plot
png(paste0(outpath,"MDS_plot.png"))
plotMDS(y_filtered_normalized, labels=groups, main="Differential expressed genes") # error: Only 2 columns of data: need at least 3
dev.off()

# 7. Create the design matrix ( Setting up the Model)

# Design matrix records which conditions/treatment were applied to each sample.
# Specify a design matrix versus an experimental design for this study 
# Design matrix assign each group to the sampleamples$group)
#designs that belong to it

design <- model.matrix(~y_filtered_normalized$samples$group)
design

# 8. Estimate the dispersion (variance)

#Estimate dispersion Using classic mode.
y_et3 <- estimateDisp(y_filtered_normalized) #combining both steps

y_et3$common.dispersion #The common dispersion is found to be 0.1499

#Plot the Estimate the dispersion (BCV: Biological coefficient of variation)
png(paste0(outpath,"est_disp.png"))
plotBCV(y_et3, main = "RNAseq experiment") #coefficient of biological variation is around 0.4
dev.off()

#The common dispersion is found to be 0.1499 and the coefficient of biological variation is around 0.4

#Estimating the overall dispersion (using CommonDisp) 
#we need to estimate the dispersion parameter for our negative binomial model
#As there are only a few samples, it is difficult to estimate the dispersion accurately for each gene

y_et <- estimateCommonDisp(y_filtered_normalized, verbose=TRUE) 
y_et1 <- estimateTagwiseDisp(y_et)
y_et$tagwise.dispersion

#Plot the Estimate Common dispersion 
png(paste0(outpath,"est_Comm_disp.png")) #No dispersion value found
plotBCV(y_et, main = "RNAseq experiment")
dev.off()

# Estimating gene-wise dispersion estimates using GLM (Generalized Linear Model) fitting
# Help us to get an overall level of biological variablility
y_et2 <- estimateGLMCommonDisp(y_et, design, verbose=T)
#output: Disp = 0.14502 , BCV = 0.3808 

y_et2 <- estimateGLMTagwiseDisp(y_et, design) 
y_et2$tagwise.dispersion

#Plot the Estimate Tagwise dispersion 
png(paste0(outpath,"est_tagwise_disp.png"))
plotBCV(y_et2, main = "RNAseq experiment")
dev.off()

# 9. Differential expression analysis (between libraries)

# Fit the linear mode using glmFit()
fit <- glmFit(y_et2, design) 
fit
names(fit)
head(coef(fit))

# Two test in edgeR (LRT and QLF)

# Conduct likelihood ratio test (lrt: use it when there is no replicate)
# Perform a test and extract the top hits
# we test for significant DE in each gene using lrt
lrt <- glmLRT(fit, coef = 2) 
edgeR_result_2 <- topTags(lrt) # top 10 DEG
edgeR_result_2  # Print top genes
topTags(lrt, n=100000) # top 100000 differential expression genes

summary(edgeR_result_2) 
write.csv(edgeR_result_2$table, file="edgeR_result_2") # extract table

# Extract table from lrt object and save it
results <- lrt$table
results$ENSEMBL <- rownames(results)
View(results)
save(results, file='edgeR_results_1') # save in compressed format
write.csv(results$table, file="edgeR_results_1") # extract table
write.table(results,file="edgeR_results_1.tsv",sep="\t",row.names=FALSE)

# Check the plot of p-values
png(paste0(outpath,"results_histogram.png"))
hist(results$PValue, breaks=50, col="blue", main = "Plot of p-values", )
dev.off()

#Quasi-likehood F-test(QLF; use it when there are replicates)
#we test for significant DE in each gene using QL F-test
fit2 <-glmQLFit(y_et2, design)
qlf <- glmQLFTest(fit2, coef=2)
edgeR_result_3 <- topTags(qlf) # top 10 DEG
edgeR_result_3  # Print top genes
topTags(qlf, n=100000)  # top 100k DEG
summary(edgeR_result_3)

# Alternative for Quasi-likehood F-test
# top genes individual by counts-per-million  
edgeR_result_4 <- rownames(topTags(qlf)) # top 10 in rows
edgeR_result_4 
cpm(y_et2)[edgeR_result_4,] # top 10 in table with cpm

#Extract table from qlf object and save it
results <- lrt$table
results$ENSEMBL <- rownames(results)
View(results)
save(results, file='results_DE_edgeR-1') #save in compressed format
write.csv(results$table, file="results_DE_edgeR") #extract table

#Note: sometimes lrt reveal top GE level more interested than QLF because p-values are somehow lower than in QLF

#plot the log-fold changes of all the genes, and the highlight those that are differentially expressed.
de2 <- decideTestsDGE(lrt, p=0.001)
de2 <- rownames(lrt)[as.logical(de2)]
de2

#Visualize the results using plotSmear
png(paste0(outpath,"log-fold_changes_genes.png"))
plotSmear(lrt, de.tags=de2, main = "RNAseq experiment")
abline(h=c(-1, 1), col=2)
dev.off()

#Adjust p-value
results$Padj <- p.adjust(results$PValue, method = "BH") #Benjamini-Hochberg method  (BH)
results <-results[order(results$Padj),]
View(results)

#Volcano plot
deg <- which(results$Padj<=0.05)
png(paste0(outpath,"volcano_plot.png"))
plot(results$logFC, -log10(results$Padj), pch=21, col="black", bg="black", xlab="log2(FoldChange)", ylab="-log10(adjusted p-value)", main = "Volcano plot")
points(results$logFC[deg], -log10(results$Padj[deg]), pch=21, col="black", bg="red", )
dev.off()

#additional option :
?decideTests

# Print number of up/down significant genes at FDR = 0.05  significance level (lrt method)
# Total number of genes significantly up-regulated (2 genes) or down-regulated (3 genes) at 5% FDR
de_lrt <- decideTestsDGE(lrt, adjust.method="BH", p=.05) # from lrt table 
summary(de_lrt)
summary(decideTests(lrt)) # produce the same results

#Plot log-fold change against log-counts per million, with DE genes highlighted
png(paste0(outpath,"genes_significantly_FDR_5_percentage_lrt.png"))
plotMD(lrt, main = "RNAseq experiment")
abline(h=c(-1, 1), col="blue")
dev.off()

# Print number of up/down significant genes at FDR = 0.05  significance level (lrt method)
summary(de_qlf <- decideTestsDGE(qlf, adjust.method="BH", p=.05)) # from qlf table 
summary(decideTests(qlf)) # produce the same results

#Plot log-fold change against log-counts per million, with DE genes highlighted
png(paste0(outpath,"genes_significantly_FDR_5_percentage_qlf.png"))
plotMD(qlf, main = "RNAseq experiment")
abline(h=c(-1, 1), col="blue")
dev.off()

#Get output with BH-adjusted FDR values - all genes, any p-value, unsorted
out_lrt <- topTags(lrt, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out_lrt
out_qlf <- topTags(qlf, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table
out_qlf

#Reference:
#Robinson, MD.et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data, Bioinformatics, 26 (1) 139-140.

#Alternative option
#Differential expression between Treated and control 
#Hypothesis testing
de <- exactTest(y_et)
de

#To display the most significant tags
topTags(de, n=10)

#Selecting deferentially expressed genes with FDR < 0.05
de_05 <- decideTestsDGE(de,p.value=0.05)

#explore the results
de_05 

summary(de_05)

#Generating a dataframe containing DE genes with FDR < 0.05
isDE <- as.logical(de_05)
de_05name <-rownames(y_et)[isDE]
de_05.table <- de[de_05name, ]

#Exporting data
#The blue lines indicate 2 fold-changes
write.csv(de_05.table$table, file="de_05")
write.csv(de$table, file="de")

#Alternative option
#Differential expression between Treated and control 

# Testing for DE between groups
de1 <- exactTest(y_et1, pair=c("Treated","Control")) 

#Extracting DEG (with topTags(de1) option will show ten DEGs)
res <- topTags(de1, n=length(de1))

#making rownames of the results to be first column of the results with column name 'Geneid' (easy for saving)
res$table[,"Geneid"] <- rownames(res$table) 

#Verify the columns of results
colnames(res$table)

#Storing selective columns to a variable 'resF'
resF <- res$table[,c("Geneid","logFC","logCPM","PValue","FDR")] 
resF

# 10. Pathway analysis

# Two option: Gene ontology analysis and KEGG pathway analysis

# Create a matrix of gene log2 fold changes
gene_matrix <- results$logFC

## Add the ENSEMBLID's as names for each logFC entry
names(gene_matrix) <- results$ENSEMBL

# View the format of the gene matrix
head(gene_matrix)

# Output: 
# AFUN020895 AFUN017558 AFUN017694 AFUN002422 AFUN015094 AFUN016451 
# -4.808196   7.613668  -5.806523  -2.571662   3.080024   5.362686 

# 10.1. Enrich genes using the KEGG database
library(KEGG.db)
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'anophelesf',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)


#This analysis was conducted on:
sessionInfo()

R version 4.0.5 (2021-03-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
  [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] KEGG.db_3.2.4        AnnotationDbi_1.52.0 IRanges_2.24.1       S4Vectors_0.28.1    
[5] Biobase_2.50.0       BiocGenerics_0.36.1  pheatmap_1.0.12      edgeR_3.32.1        
[9] limma_3.46.0        

loaded via a namespace (and not attached):
  [1] Rcpp_1.0.6                  locfit_1.5-9.4              lattice_0.20-41            
[4] Rsamtools_2.6.0             Biostrings_2.58.0           R6_2.5.0                   
[7] GenomeInfoDb_1.26.7         RSQLite_2.2.7               httr_1.4.2                 
[10] zlibbioc_1.36.0             rlang_0.4.11                rstudioapi_0.13            
[13] annotate_1.68.0             blob_1.2.1                  Matrix_1.3-2               
[16] splines_4.0.5               BiocParallel_1.24.1         RCurl_1.98-1.3             
[19] bit_4.0.4                   munsell_0.5.0               tinytex_0.31               
[22] DelayedArray_0.16.3         compiler_4.0.5              rtracklayer_1.50.0         
[25] xfun_0.22                   pkgconfig_2.0.3             SummarizedExperiment_1.20.0
[28] GenomeInfoDbData_1.2.4      matrixStats_0.58.0          XML_3.99-0.6               
[31] crayon_1.4.1                withr_2.4.2                 GenomicAlignments_1.26.0   
[34] bitops_1.0-7                grid_4.0.5                  xtable_1.8-4               
[37] gtable_0.3.0                lifecycle_1.0.0             DBI_1.1.1                  
[40] scales_1.1.1                cli_2.5.0                   cachem_1.0.4               
[43] XVector_0.30.0              genefilter_1.72.1           vctrs_0.3.8                
[46] RColorBrewer_1.1-2          tools_4.0.5                 bit64_4.0.5                
[49] MatrixGenerics_1.2.1        fastmap_1.1.0               survival_3.2-10            
[52] colorspace_2.0-0            BiocManager_1.30.12         GenomicRanges_1.42.0       
[55] sessioninfo_1.1.1           memoise_2.0.0              
> 



