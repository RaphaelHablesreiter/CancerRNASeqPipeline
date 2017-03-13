#!/usr/bin/env Rscript

#*****************************************************************************************
# Libraries
#*****************************************************************************************

#Loading necassary libraries
library(cummeRbund)
library(plyr)

#*****************************************************************************************
# Initialisation
#*****************************************************************************************

#cmdl_input for custom use of this R-script
cmdl_input <- array(1:12)
cmdl_input[6] = "Normal"
cmdl_input[7] = "Adenocarcinoma"
cmdl_input[8] = "/Volumes/HablesExtern/PRJNA163279PE/tmp/cuffdiff_output"
cmdl_input[9] = "/Volumes/HablesExtern/PRJNA163279PE/tmp"
cmdl_input[10] = "/Volumes/HablesExtern/PRJNA163279PE/Graphics"
cmdl_input[11] = "/Volumes/HablesExtern/PRJNA163279PE/Results"


#Setting input parameters
#cmdl_input[6] is the name of sequence 1
#cmdl_input[7] is the name of sequence 2
#cmdl_input[8] is the folder of the cuffdiff output
#cmdl_input[9] is the folder of the wd(=temporary folder)
#cmdl_input[10] is the folder of the graphics
#cmdl_input[11] is the folder of the results

#cmdl_input = commandArgs(comtrailingOnly = FALSE)

#Global parameters
barplot_height = 6
prefix = paste(cmdl_input[6],"Vs",cmdl_input[7],sep = "")
setwd(cmdl_input[9])

#Open Cuffdiff output
cuffdata <- readCufflinks(cmdl_input[8])
number_rep = max(replicates(cuffdata)$replicate) + 1

#Generating a DispersionPlot
dispersionPlot <- dispersionPlot(genes(cuffdata)) + xlab("Count") + ylab("Dispersion") + labs(colour = "Condition")
plotname = paste(cmdl_input[10], "/", prefix, "_DispersionPlot.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(dispersionPlot)
dev.off()


#*****************************************************************************************
# DensityPlots
#*****************************************************************************************

#Genes densityplot from both conditions
densityPlotGenes <- csDensity(genes(cuffdata)) + xlab("Log10(FPKM)") + ylab("Density") + ggtitle("Genes")
densityPlotGenes$labels$fill = "Condition"
densityPlotGenes$labels$colour = "Condition"
plotname = paste(cmdl_input[10], "/", prefix, "_DensityPlotGenes.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(densityPlotGenes)
dev.off()

#Genes densityplot from all replicates
densityPlotGenesReplicates <- csDensity(genes(cuffdata), replicates = TRUE) + xlab("Log10(FPKM)") + ylab("Density") + ggtitle("Genes")
densityPlotGenesReplicates$labels$fill = "Condition"
densityPlotGenesReplicates$labels$colour = "Condition"
plotname = paste(cmdl_input[10], "/", prefix, "_DensityPlotGenesReplicates.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(densityPlotGenesReplicates)
dev.off()

#Isoform densityplot from both conditions
densityPlotIsoforms <- csDensity(isoforms(cuffdata)) + xlab("Log10(FPKM)") + ylab("Density") + ggtitle("Isoforms")
densityPlotIsoforms$labels$fill = "Condition"
densityPlotIsoforms$labels$colour = "Condition"
plotname = paste(cmdl_input[10], "/", prefix, "_DensityPlotIsoforms.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(densityPlotIsoforms)
dev.off()

#Isoform densityplot from all replicates
densityPlotIsoformsReplicates <- csDensity(isoforms(cuffdata), replicates = TRUE) + xlab("Log10(FPKM)") + ylab("Density") + ggtitle("Isoforms")
densityPlotIsoformsReplicates$labels$fill = "Condition"
densityPlotIsoformsReplicates$labels$colour = "Condition"
plotname = paste(cmdl_input[10], "/", prefix, "_DensityPlotIsoformsReplicates.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(densityPlotIsoformsReplicates)
dev.off()

#*****************************************************************************************
# BoxPlots
#*****************************************************************************************

#Genes boxplot from both conditions
boxPlotGenes <- csBoxplot(genes(cuffdata)) + xlab("Condition") + ylab("Log10(FPKM)") + ggtitle("Genes") + scale_fill_hue("Condition",l=50,h.start=200)
plotname = paste(cmdl_input[10], "/", prefix, "_BoxPlotGenes.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(boxPlotGenes)
dev.off()

#Genes boxplot from all replicates
boxPlotGenesReplicates <- csBoxplot(genes(cuffdata), replicates = TRUE) + xlab("Condition") + ylab("Log10(FPKM)") + ggtitle("Genes") + scale_fill_hue("Condition",l=50,h.start=200)
plotname = paste(cmdl_input[10], "/", prefix, "_BoxPlotGenesReplicates.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(boxPlotGenesReplicates)
dev.off()

#Isoform densityplot from both conditions
boxPlotIsoforms <- csBoxplot(genes(cuffdata)) + xlab("Condition") + xlab("Condition") + ylab("Log10(FPKM)") + ggtitle("Genes") + scale_fill_hue("Condition",l=50,h.start=200)
plotname = paste(cmdl_input[10], "/", prefix, "_BoxPlotIsoforms.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(boxPlotIsoforms)
dev.off()

#Isoform boxplot from all replicates
boxPlotIsoformsReplicates <- csBoxplot(genes(cuffdata), replicates = TRUE) + xlab("Condition") + ylab("Log10(FPKM)") + ggtitle("Genes") + scale_fill_hue("Condition",l=50,h.start=200)
plotname = paste(cmdl_input[10], "/", prefix, "_BoxPlotIsoformsReplicates.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(boxPlotIsoformsReplicates)
dev.off()

#*****************************************************************************************
# MA-Plot (Genes and Isoforms)
#*****************************************************************************************


maPlotGenes <- MAplot(genes(cuffdata), cmdl_input[6], cmdl_input[7], smooth = TRUE)
plotname = paste(cmdl_input[10], "/", prefix, "_MAPlotGenes.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(maPlotGenes)
dev.off()

maPlotIsoforms <- MAplot(isoforms(cuffdata), cmdl_input[6], cmdl_input[7], smooth = TRUE)
plotname = paste(cmdl_input[10], "/", prefix, "_MAPlotIsoforms.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(maPlotIsoforms)
dev.off()

#*****************************************************************************************
# Volcano-Plot (Genes and Isoforms)
#*****************************************************************************************


volcanoPlotGenes <- csVolcano(genes(cuffdata), cmdl_input[6], cmdl_input[7], showSignificant = TRUE, alpha = .05)
volcanoPlotGenes$labels$colour = "Significant"
plotname = paste(cmdl_input[10], "/", prefix, "_VolcanoPlotGenes.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
#png(file = plotname, width = 1000, height = 1000, res = 300)
plot(volcanoPlotGenes)
dev.off()

volcanoPlotIsoforms <- csVolcano(isoforms(cuffdata), cmdl_input[6], cmdl_input[7], showSignificant = TRUE, alpha = .05)
volcanoPlotIsoforms$labels$colour = "Significant"
plotname = paste(cmdl_input[10], "/", prefix, "_VolcanoPlotIsoforms.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
#png(file = plotname, width = 1000, height = 1000, res = 300)
plot(volcanoPlotIsoforms)
dev.off()

#*****************************************************************************************
# Scatter-Plot (Genes and Isoforms)
#*****************************************************************************************


scatterPlotGenes <- csScatter(genes(cuffdata), cmdl_input[6], cmdl_input[7], smooth = TRUE) + ggtitle("Genes")
plotname = paste(cmdl_input[10], "/", prefix, "_ScatterPlotGenes.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
#png(file = plotname, width = 1000, height = 1000, res = 300)
plot(scatterPlotGenes)
dev.off()

scatterPlotIsoforms <- csScatter(isoforms(cuffdata), cmdl_input[6], cmdl_input[7], smooth = TRUE) + ggtitle("Isoforms ")
plotname = paste(cmdl_input[10], "/", prefix, "_ScatterPlotIsoforms.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
#png(file = plotname, width = 1000, height = 1000, res = 300)
plot(scatterPlotIsoforms)
dev.off()

# #Generating genes and isoforms ScattermatrixPlots
# scatterMatrixGenesReplicates <- csScatterMatrix(genes(cuffdata), replicates = TRUE)
# scatterMatrixGenesReplicates$labels$x = "log10(FPKM)"
# scatterMatrixGenesReplicates$labels$y = "log10(FPKM)"
# scatterMatrixIsoformsReplicates <- csScatterMatrix(isoforms(cuffdata), replicates = TRUE)
# scatterMatrixIsoformsReplicates$labels$x = "log10(FPKM)"
# scatterMatrixIsoformsReplicates$labels$y = "log10(FPKM)"
# 
# plotname = paste(cmdl_input[10], "/", prefix, "_ScatterMatrixGenesReplicates.pdf", sep = "")
# png(file = plotname,  width = 10*number_rep, height = 10*number_rep, res = 300)
# plot(scatterMatrixGenesReplicates)
# dev.off()
# plotname = paste(cmdl_input[10], "/", prefix, "_ScatterMatrixIsoformsReplicates.pdf", sep = "")
# png(file = plotname, width = 10*number_rep, height = 10*number_rep, res = 300)
# plot(scatterMatrixIsoformsReplicates)
# dev.off()


#*****************************************************************************************
# Significantly DE genes and isoforms
#*****************************************************************************************

#Get significantly DE genes
diffGeneIDs <- getSig(cuffdata, alpha = 0.05, level = 'genes')
diffGenes <- getGenes(cuffdata,diffGeneIDs)
GeneNames <- featureNames(diffGenes)

row.names(GeneNames) = GeneNames$tracking_id

diffGenesNames <- as.matrix(GeneNames)
diffGenesNames <- diffGenesNames[,-1]

diffGenesData <- diffData(diffGenes)
row.names(diffGenesData) = diffGenesData$gene_id
diffGenesData <- diffGenesData[,-1]

diffGenesOutput <- merge(diffGenesNames,diffGenesData,by="row.names")
names(diffGenesOutput)[1] <- "gene_id"
names(diffGenesOutput)[2] <- "gene_name"

diffGenesOutput <- diffGenesOutput[!is.na(diffGenesOutput$test_stat),]
diffGenesOutput <- diffGenesOutput[!is.na(diffGenesOutput$gene_name),]

#Introduction of an edited log2-FC to avoid infinite values 
diffGenesOutput$log2_fold_change_finite <- 0
diffGenesOutput$log2_fold_change_finite = log2((diffGenesOutput$value_2 + .1)/(diffGenesOutput$value_1 + .1))

names(diffGenesOutput)[6] <- paste("FPKM_", cmdl_input[6])
names(diffGenesOutput)[7] <- paste("FPKM_", cmdl_input[7])
names(diffGenesOutput)[11] <- paste("adjusted_p_value")

#Get FPKM values for all replicates in both condtions
diffGenesFPKMreplicates <- repFpkmMatrix(genes(cuffdata))
diffGenesFPKMreplicatesOutput <- merge(diffGenesNames, diffGenesFPKMreplicates, by="row.names")
names(diffGenesFPKMreplicatesOutput)[1] <- "gene_id"
names(diffGenesFPKMreplicatesOutput)[2] <- "gene_name"

#Adding gene-locus to the DE genes tables
geneInfo <- annotation(diffGenes)
geneInfo <- subset(geneInfo, select = c(gene_id, locus))
names(geneInfo)[2] <- "gene_locus"

diffGenesOutputWithLocus <- merge(diffGenesOutput, geneInfo, by="gene_id")
diffGenesOutputWithLocus <- diffGenesOutputWithLocus[, c(1,2,14,3,4,6:8,13,10,11,9,5,12)]
diffGenesFPKMreplicatesOutputWithLocus <- merge(diffGenesFPKMreplicatesOutput, geneInfo, by="gene_id")
diffGenesFPKMreplicatesOutputWithLocus <- diffGenesFPKMreplicatesOutputWithLocus[, c(1:2,ncol(diffGenesFPKMreplicatesOutputWithLocus),3:(ncol(diffGenesFPKMreplicatesOutputWithLocus)-1))]


#Get significantly DE isoforms
diffIsoformsIDs <- getSig(cuffdata, alpha = 0.05, level = 'isoforms')
diffIsoforms <- getFeatures(cuffdata, diffIsoformsIDs, level = 'isoforms')
IsoformNames <- featureNames(diffIsoforms)

row.names(IsoformNames) = IsoformNames$tracking_id
diffIsoformsNames <- as.matrix(IsoformNames)
diffIsoformsNames <- diffIsoformsNames[,-1]

diffIsoformsData <- diffData(diffIsoforms)
row.names(diffIsoformsData) = diffIsoformsData$isoform_id
diffIsoformsData <- diffIsoformsData[,-1]

diffIsoformsOutput <- merge(diffIsoformsNames,diffIsoformsData,by="row.names")
names(diffIsoformsOutput)[1] <- "isoform_id"
names(diffIsoformsOutput)[2] <- "isoform_name"

diffIsoformsOutput <- diffIsoformsOutput[!is.na(diffIsoformsOutput$test_stat),]

#Introduction of an edited log2-FC to avoid infinite values 
diffIsoformsOutput$log2_fold_change_finite <- 0
diffIsoformsOutput$log2_fold_change_finite = log2((diffIsoformsOutput$value_2 + .1)/(diffIsoformsOutput$value_1 + .1))

names(diffIsoformsOutput)[6] <- paste("FPKM_", cmdl_input[6])
names(diffIsoformsOutput)[7] <- paste("FPKM_", cmdl_input[7])
names(diffIsoformsOutput)[11] <- paste("adjusted_p_value")

#Adding isoform-locus to the DE isoforms tables
isoformInfo <- annotation(diffIsoforms)
isoformInfo <- subset(isoformInfo, select = c(isoform_id, locus))
names(isoformInfo)[2] <- "isoform_locus"
diffIsoformsOutputWithLocus <- merge(diffIsoformsOutput, isoformInfo, by="isoform_id")
diffIsoformsOutputWithLocus <- diffIsoformsOutputWithLocus[, c(1,2,14,3,4,6:8,13,10,11,9,5,12)]

#Export tables as *.csv-files
csvname = paste(cmdl_input[11], "/", prefix, "_genes.csv", sep = "")
write.csv(diffGenesOutputWithLocus, file = csvname, row.names = FALSE, na = "NA")
csvname = paste(cmdl_input[11], "/", prefix, "_isoforms.csv", sep = "")
write.csv(diffIsoformsOutputWithLocus, file = csvname, row.names = FALSE, na = "NA")
csvname = paste(cmdl_input[11], "/", prefix, "_genesFPKMreplicates.csv", sep = "")
write.csv(diffGenesFPKMreplicatesOutputWithLocus, file = csvname, row.names = FALSE, na = "NA")

#*****************************************************************************************
# Information for Barplots and Heatmaps of DE genes and isoforms
#*****************************************************************************************

#Create table diffGenesOutput with criteria of an adjusted(p-value)<0.05 and log2-FC > 1
diffGenesOutput <- diffGenesOutput[!is.na(diffGenesOutput$gene_name),]
diffGenesOutput_std = diffGenesOutput
diffGenesOutput = diffGenesOutput[diffGenesOutput$adjusted_p_value < 0.05 & abs(diffGenesOutput$log2_fold_change) > 1,]

diffGenesOutputWithLocus <- merge(diffGenesOutput, geneInfo, by="gene_id")
diffGenesOutputWithLocus <- diffGenesOutputWithLocus[, c(1,2,14,3,4,6:8,13,10,11,9,5,12)]

#Export table of genes that meet the criteria of an adjusted(p-value)<0.05 and log2-FC > 1
csvname = paste(cmdl_input[11], "/", prefix, "_genes_adjPvalueLT05_Log2FcGT1.csv", sep = "")
write.csv(diffGenesOutputWithLocus, file = csvname, row.names = FALSE, na = "NA")

#Separate genes (STD) into two groups (log2-FC < 0 and log2-FC > 0)
diffGenesOutputSortedAsc_std = arrange(diffGenesOutput_std,log2_fold_change_finite)
diffGenesOutputSortedAsc_std = diffGenesOutputSortedAsc_std[diffGenesOutputSortedAsc_std[8] < 0,]
diffGenesOutputSortedDesc_std = arrange(diffGenesOutput_std,desc(log2_fold_change_finite))
diffGenesOutputSortedDesc_std = diffGenesOutputSortedDesc_std[diffGenesOutputSortedDesc_std[8] > 0,]

#Get TOP50 genes (STD)
if(nrow(diffGenesOutputSortedAsc_std) < 50 ){
  diffGenesOutputTop50UpCondA_std <- diffGenesOutputSortedAsc_std[1:nrow(diffGenesOutputSortedAsc_std),]
}else{
  diffGenesOutputTop50UpCondA_std <- diffGenesOutputSortedAsc_std[1:50,] 
}  
if(nrow(diffGenesOutputSortedDesc_std) < 50 ){
  diffGenesOutputTop50UpCondB_std <- diffGenesOutputSortedDesc_std[1:nrow(diffGenesOutputSortedDesc_std),]
}else{
  diffGenesOutputTop50UpCondB_std <- diffGenesOutputSortedDesc_std[1:50,] 
}  

diffGenesOutputTop50UpCondA_std = arrange(diffGenesOutputTop50UpCondA_std,desc(diffGenesOutputTop50UpCondA_std[6]))
diffGenesOutputTop50UpCondB_std = arrange(diffGenesOutputTop50UpCondB_std,desc(diffGenesOutputTop50UpCondB_std[7]))
diffGenesOutput_std = arrange(diffGenesOutput_std,diffGenesOutput_std[2])

#Separate isoforms (STD) into two groups (log2-FC < 0 and log2-FC > 0)
diffIsoformsOutput <- diffIsoformsOutput[!is.na(diffIsoformsOutput$isoform_name),]

diffIsoformsOutputSortedAsc_std = arrange(diffIsoformsOutput,log2_fold_change_finite)
diffIsoformsOutputSortedAsc_std = diffIsoformsOutputSortedAsc_std[diffIsoformsOutputSortedAsc_std[8] < 0,]
diffIsoformsOutputSortedDesc_std = arrange(diffIsoformsOutput,desc(log2_fold_change_finite))
diffIsoformsOutputSortedDesc_std = diffIsoformsOutputSortedDesc_std[diffIsoformsOutputSortedDesc_std[8] > 0,]

#Get TOP50 isoforms (STD)
if(nrow(diffIsoformsOutputSortedAsc_std) < 50 ){
  diffIsoformsOutputTop50UpCondA_std <- diffIsoformsOutputSortedAsc_std[1:nrow(diffIsoformsOutputSortedAsc_std),]
}else{
  diffIsoformsOutputTop50UpCondA_std <- diffIsoformsOutputSortedAsc_std[1:50,] 
}  
if(nrow(diffIsoformsOutputSortedDesc_std) < 50 ){
  diffIsoformsOutputTop50UpCondB_std <- diffIsoformsOutputSortedDesc_std[1:nrow(diffIsoformsOutputSortedDesc),]
}else{
  diffIsoformsOutputTop50UpCondB_std <- diffIsoformsOutputSortedDesc_std[1:50,] 
}  

diffIsoformsOutputTop50UpCondA_std = arrange(diffIsoformsOutputTop50UpCondA_std,desc(diffIsoformsOutputTop50UpCondA_std[6]))
diffIsoformsOutputTop50UpCondB_std = arrange(diffIsoformsOutputTop50UpCondB_std,desc(diffIsoformsOutputTop50UpCondB_std[7]))
diffIsoformsOutput = arrange(diffIsoformsOutput,diffIsoformsOutput[2])
diffIsoformsOutput_std = diffIsoformsOutput

#Get information for barplot
sigGenesPlot_std = getGenes(cuffdata, diffGenesOutput_std$gene_id)
sigGenesCondAUpPlot_std = getGenes(cuffdata, diffGenesOutputTop50UpCondA_std$gene_id)
sigGenesCondBUpPlot_std = getGenes(cuffdata, diffGenesOutputTop50UpCondB_std$gene_id)

sigIsoformsPlot_std = getFeatures(cuffdata,diffIsoformsOutput$isoform_id,level='isoforms')
sigIsoformsCondAUpPlot_std = getFeatures(cuffdata,diffIsoformsOutputTop50UpCondA_std$isoform_id,level='isoforms')
sigIsoformsCondBUpPlot_std = getFeatures(cuffdata,diffIsoformsOutputTop50UpCondB_std$isoform_id,level='isoforms')


#Separate genes (adjusted p-value < 0.05 and abs(log2FC) > 1) 
#    into two groups (log2-FC < 0 and log2-FC > 0)
diffGenesOutputSortedAsc = arrange(diffGenesOutput,log2_fold_change_finite)
diffGenesOutputSortedAsc = diffGenesOutputSortedAsc[diffGenesOutputSortedAsc[8] < 0,]
diffGenesOutputSortedDesc = arrange(diffGenesOutput,desc(log2_fold_change_finite))
diffGenesOutputSortedDesc = diffGenesOutputSortedDesc[diffGenesOutputSortedDesc[8] > 0,]

#Get TOP50 genes (adjusted p-value < 0.05 and abs(log2FC) > 1)
if(nrow(diffGenesOutputSortedAsc) < 50 ){
  diffGenesOutputTop50UpCondA <- diffGenesOutputSortedAsc[1:nrow(diffGenesOutputSortedAsc),]
}else{
  diffGenesOutputTop50UpCondA <- diffGenesOutputSortedAsc[1:50,] 
}  
if(nrow(diffGenesOutputSortedDesc) < 50 ){
  diffGenesOutputTop50UpCondB <- diffGenesOutputSortedDesc[1:nrow(diffGenesOutputSortedDesc),]
}else{
  diffGenesOutputTop50UpCondB <- diffGenesOutputSortedDesc[1:50,] 
}  

diffGenesOutputTop50UpCondA = arrange(diffGenesOutputTop50UpCondA,desc(diffGenesOutputTop50UpCondA[6]))
diffGenesOutputTop50UpCondB = arrange(diffGenesOutputTop50UpCondB,desc(diffGenesOutputTop50UpCondB[7]))
diffGenesOutput = arrange(diffGenesOutput,diffGenesOutput[2])

sigGenesPlot = getGenes(cuffdata, diffGenesOutput$gene_id)
sigGenesCondAUpPlot = getGenes(cuffdata, diffGenesOutputTop50UpCondA$gene_id)
sigGenesCondBUpPlot = getGenes(cuffdata, diffGenesOutputTop50UpCondB$gene_id)

#Get isoforms meeting the criteria of adjusted p-value < 0.05 and abs(log2FC) > 1
diffIsoformsOutput <- diffIsoformsOutput[!is.na(diffIsoformsOutput$isoform_name),]
diffIsoformsOutput = diffIsoformsOutput[diffIsoformsOutput$adjusted_p_value < 0.05 & abs(diffIsoformsOutput$log2_fold_change) > 1,]

diffIsoformsOutputWithLocus <- merge(diffIsoformsOutput, isoformInfo, by="isoform_id")
diffIsoformsOutputWithLocus <- diffIsoformsOutputWithLocus[, c(1,2,14,3,4,6:8,13,10,11,9,5,12)]

#Export as *.csv-file
csvname = paste(cmdl_input[11], "/", prefix, "_isoforms_adjPvalueLT05_Log2FcGT1.csv", sep = "")
write.csv(diffIsoformsOutputWithLocus, file = csvname, row.names = FALSE, na = "NA")

#Separate isoforms (adjusted p-value < 0.05 and abs(log2FC) > 1) 
#    into two groups (log2-FC < 0 and log2-FC > 0)
diffIsoformsOutputSortedAsc = arrange(diffIsoformsOutput,log2_fold_change_finite)
diffIsoformsOutputSortedAsc = diffIsoformsOutputSortedAsc[diffIsoformsOutputSortedAsc[8] < 0,]
diffIsoformsOutputSortedDesc = arrange(diffIsoformsOutput,desc(log2_fold_change_finite))
diffIsoformsOutputSortedDesc = diffIsoformsOutputSortedDesc[diffIsoformsOutputSortedDesc[8] > 0,]

#Get TOP50 isoforms (adjusted p-value < 0.05 and abs(log2FC) > 1)
if(nrow(diffIsoformsOutputSortedAsc) < 50 ){
  diffIsoformsOutputTop50UpCondA <- diffIsoformsOutputSortedAsc[1:nrow(diffIsoformsOutputSortedAsc),]
}else{
  diffIsoformsOutputTop50UpCondA <- diffIsoformsOutputSortedAsc[1:50,] 
}  
if(nrow(diffIsoformsOutputSortedDesc) < 50 ){
  diffIsoformsOutputTop50UpCondB <- diffIsoformsOutputSortedDesc[1:nrow(diffIsoformsOutputSortedDesc),]
}else{
  diffIsoformsOutputTop50UpCondB <- diffIsoformsOutputSortedDesc[1:50,] 
}  

diffIsoformsOutputTop50UpCondA = arrange(diffIsoformsOutputTop50UpCondA,desc(diffIsoformsOutputTop50UpCondA[6]))
diffIsoformsOutputTop50UpCondB = arrange(diffIsoformsOutputTop50UpCondB,desc(diffIsoformsOutputTop50UpCondB[7]))
diffIsoformsOutput = arrange(diffIsoformsOutput,diffIsoformsOutput[2])

#Get information for barplot
sigIsoformsPlot = getFeatures(cuffdata,diffIsoformsOutput$isoform_id,level='isoforms')
sigIsoformsCondAUpPlot = getFeatures(cuffdata,diffIsoformsOutputTop50UpCondA$isoform_id,level='isoforms')
sigIsoformsCondBUpPlot = getFeatures(cuffdata,diffIsoformsOutputTop50UpCondB$isoform_id,level='isoforms')

#*****************************************************************************************
# Barplots and Heatmaps of DE genes and isoforms
#*****************************************************************************************

#Heatmaps and BarPlots for significant Genes and significant Isoforms STD
genesHeatmap <- csHeatmap(sigGenesPlot_std, clustering = 'none')
genesHeatmapCondAUp <- csHeatmap(sigGenesCondAUpPlot_std, clustering = 'both')
genesHeatmapCondBUp <- csHeatmap(sigGenesCondBUpPlot_std, clustering = 'both')

genesBarplot <- expressionBarplot(sigGenesPlot_std, pseudocount = 1, logMode = TRUE)
genesBarplot$labels$fill = "Sample Name"
genesBarplot$labels$x = "Gene Name"
genesBarplot$labels$y = "log(FPKM + 1)"
genesBarplot <- genesBarplot + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffGenesOutput_std$gene_id, labels = sub("*\\,.*", "", diffGenesOutput_std$gene_name))

genesBarplotCondAUp <- expressionBarplot(sigGenesCondAUpPlot_std, pseudocount = 1, logMode = TRUE)
genesBarplotCondAUp$labels$fill = "Sample Name"
genesBarplotCondAUp$labels$x = "Gene Name"
genesBarplotCondAUp$labels$y = "log(FPKM + 1)"
genesBarplotCondAUp <- genesBarplotCondAUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffGenesOutputTop50UpCondA_std$gene_id, labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondA_std$gene_name))

genesBarplotLog2FCCondAUp <- ggplot(data = diffGenesOutputTop50UpCondA_std, aes(x = reorder(gene_id,log2_fold_change_finite), y = -log2_fold_change_finite)) + geom_bar(stat="identity", fill = "steelblue3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondA_std$gene_name))

genesBarplotCondBUp <- expressionBarplot(sigGenesCondBUpPlot_std, pseudocount = 1, logMode = TRUE)
genesBarplotCondBUp$labels$fill = "Sample Name"
genesBarplotCondBUp$labels$x = "Gene Name"
genesBarplotCondBUp$labels$y = "log(FPKM + 1)"
genesBarplotCondBUp <- genesBarplotCondBUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffGenesOutputTop50UpCondB_std$gene_id, labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondB_std$gene_name))

genesBarplotLog2FCCondBUp <- ggplot(data = diffGenesOutputTop50UpCondB_std, aes(x = reorder(gene_id,desc(log2_fold_change_finite)), y = log2_fold_change_finite)) + geom_bar(stat="identity", fill = "chocolate3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondB_std$gene_name))


isoformsBarplot <- expressionBarplot(sigIsoformsPlot_std, pseudocount = 1, logMode = TRUE)
isoformsBarplot$labels$fill = "Sample Name"
isoformsBarplot$labels$x = "Gene Name"
isoformsBarplot$labels$y = "log(FPKM + 1)"
isoformsBarplot <- isoformsBarplot + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffIsoformsOutput_std$isoform_id, labels = sub("*\\,.*", "", diffIsoformsOutput_std$isoform_name))

isoformsBarplotCondAUp <- expressionBarplot(sigIsoformsCondAUpPlot_std, pseudocount = 1, logMode = TRUE)
isoformsBarplotCondAUp$labels$fill = "Sample Name"
isoformsBarplotCondAUp$labels$x = "Gene Name"
isoformsBarplotCondAUp$labels$y = "log(FPKM + 1)"
isoformsBarplotCondAUp <- isoformsBarplotCondAUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffIsoformsOutputTop50UpCondA_std$isoform_id, labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondA_std$isoform_name))

isoformsBarplotLog2FCCondAUp <- ggplot(data = diffIsoformsOutputTop50UpCondA_std, aes(x = reorder(isoform_id,log2_fold_change_finite), y = -log2_fold_change_finite)) + geom_bar(stat="identity", fill = "steelblue3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondA_std$isoform_name))

isoformsBarplotCondBUp <- expressionBarplot(sigIsoformsCondBUpPlot_std, pseudocount = 1, logMode = TRUE)
isoformsBarplotCondBUp$labels$fill = "Sample Name"
isoformsBarplotCondBUp$labels$x = "Gene Name"
isoformsBarplotCondBUp$labels$y = "log(FPKM + 1)"
isoformsBarplotCondBUp <- isoformsBarplotCondBUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffIsoformsOutputTop50UpCondB_std$isoform_id, labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondB_std$isoform_name))

isoformsBarplotLog2FCCondBUp <- ggplot(data = diffIsoformsOutputTop50UpCondB_std, aes(x = reorder(isoform_id,desc(log2_fold_change_finite)), y = log2_fold_change_finite)) + geom_bar(stat="identity", fill = "chocolate3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondB_std$isoform_name))


lengthgenes = nrow(diffGenesOutput_std)
lengthisoforms = nrow(diffIsoformsOutput)

# plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesHeatmap.pdf", sep = "")
# pdf(file = plotname, width = 10, height = .3*lengthgenes)
# plot(genesHeatmap)
# dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesHeatmap_TOP50_", cmdl_input[6], "_UP.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(genesHeatmapCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesHeatmap_TOP50_", cmdl_input[7], "_UP.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(genesHeatmapCondBUp)
dev.off()

plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot.pdf", sep = "")
pdf(file = plotname, width = .25*lengthgenes , height = 15)
plot(genesBarplot)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_", cmdl_input[6], "_UP.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_", cmdl_input[7], "_UP.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotCondBUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_Log2FC_", cmdl_input[6], "_UP.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotLog2FCCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_Log2FC_", cmdl_input[7], "_UP.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotLog2FCCondBUp)
dev.off()


plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot.pdf", sep = "")
pdf(file = plotname, width = .25*lengthgenes , height = 10)
plot(isoformsBarplot)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_",cmdl_input[6],"_UP.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_",cmdl_input[7],"_UP.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotCondBUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_Log2FC_",cmdl_input[6],"_UP.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotLog2FCCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_Log2FC_",cmdl_input[7],"_UP.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotLog2FCCondBUp)
dev.off()

#Heatmaps and BarPlots for significant Genes and significant Isoforms
###adjusted p-value < 0.05 and abs(log2FC) > 1
genesHeatmap <- csHeatmap(sigGenesPlot, clustering = 'none')
genesHeatmapCondAUp <- csHeatmap(sigGenesCondAUpPlot, clustering = 'both')
genesHeatmapCondBUp <- csHeatmap(sigGenesCondBUpPlot, clustering = 'both')

genesBarplot <- expressionBarplot(sigGenesPlot, pseudocount = 1, logMode = TRUE)
genesBarplot$labels$fill = "Sample Name"
genesBarplot$labels$x = "Gene Name"
genesBarplot$labels$y = "log(FPKM + 1)"
genesBarplot <- genesBarplot + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffGenesOutput$gene_id, labels = sub("*\\,.*", "", diffGenesOutput$gene_name))

genesBarplotCondAUp <- expressionBarplot(sigGenesCondAUpPlot)
genesBarplotCondAUp$labels$fill = "Sample Name"
genesBarplotCondAUp$labels$x = "Gene Name"
genesBarplotCondAUp$labels$y = "log(FPKM + 1)"
genesBarplotCondAUp <- genesBarplotCondAUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffGenesOutputTop50UpCondA$gene_id, labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondA$gene_name))

genesBarplotLog2FCCondAUp <- ggplot(data = diffGenesOutputTop50UpCondA, aes(x = reorder(gene_id,log2_fold_change_finite), y = -log2_fold_change_finite)) + geom_bar(stat="identity", fill = "steelblue3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondA$gene_name))

genesBarplotCondBUp <- expressionBarplot(sigGenesCondBUpPlot)
genesBarplotCondBUp$labels$fill = "Sample Name"
genesBarplotCondBUp$labels$x = "Gene Name"
genesBarplotCondBUp$labels$y = "log(FPKM + 1)"
genesBarplotCondBUp <- genesBarplotCondBUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffGenesOutputTop50UpCondB$gene_id, labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondB$gene_name))

genesBarplotLog2FCCondBUp <- ggplot(data = diffGenesOutputTop50UpCondB, aes(x = reorder(gene_id,desc(log2_fold_change_finite)), y = log2_fold_change_finite)) + geom_bar(stat="identity", fill = "chocolate3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffGenesOutputTop50UpCondB$gene_name))


isoformsBarplot <- expressionBarplot(sigIsoformsPlot, pseudocount = 1, logMode = TRUE)
isoformsBarplot$labels$fill = "Sample Name"
isoformsBarplot$labels$x = "Gene Name"
isoformsBarplot$labels$y = "log(FPKM + 1)"
isoformsBarplot <- isoformsBarplot + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffIsoformsOutput$isoform_id, labels = sub("*\\,.*", "", diffIsoformsOutput$isoform_name))

isoformsBarplotCondAUp <- expressionBarplot(sigIsoformsCondAUpPlot, pseudocount = 1, logMode = TRUE)
isoformsBarplotCondAUp$labels$fill = "Sample Name"
isoformsBarplotCondAUp$labels$x = "Gene Name"
isoformsBarplotCondAUp$labels$y = "log(FPKM + 1)"
isoformsBarplotCondAUp <- isoformsBarplotCondAUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffIsoformsOutputTop50UpCondA$isoform_id, labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondA$isoform_name))

isoformsBarplotLog2FCCondAUp <- ggplot(data = diffIsoformsOutputTop50UpCondA, aes(x = reorder(isoform_id,log2_fold_change_finite), y = -log2_fold_change_finite)) + geom_bar(stat="identity", fill = "steelblue3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondA$isoform_name))

isoformsBarplotCondBUp <- expressionBarplot(sigIsoformsCondBUpPlot, pseudocount = 1, logMode = TRUE)
isoformsBarplotCondBUp$labels$fill = "Sample Name"
isoformsBarplotCondBUp$labels$x = "Gene Name"
isoformsBarplotCondBUp$labels$y = "log(FPKM + 1)"
isoformsBarplotCondBUp <- isoformsBarplotCondBUp + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_x_discrete(limits = diffIsoformsOutputTop50UpCondB$isoform_id, labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondB$isoform_name))

isoformsBarplotLog2FCCondBUp <- ggplot(data = diffIsoformsOutputTop50UpCondB, aes(x = reorder(isoform_id,desc(log2_fold_change_finite)), y = log2_fold_change_finite)) + geom_bar(stat="identity", fill = "chocolate3", width = 0.7) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Gene Name") + ylab("Log2-FC") + scale_x_discrete(labels = sub("*\\,.*", "", diffIsoformsOutputTop50UpCondB$isoform_name))


lengthgenes = nrow(diffGenesOutput)
lengthisoforms = nrow(diffIsoformsOutput)

# plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesHeatmap_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
# pdf(file = plotname, width = 10, height = .3*lengthgenes)
# plot(genesHeatmap)
# dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesHeatmap_TOP50_", cmdl_input[6], "_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 10, height = 10)
plot(genesHeatmapCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_", cmdl_input[7], "_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 10, height = barplot_height)
plot(genesHeatmapCondBUp)
dev.off()

plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = .25*lengthgenes , height = 15)
plot(genesBarplot)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_", cmdl_input[6], "_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_", cmdl_input[7], "_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotCondBUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_Log2FC_", cmdl_input[6], "_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotLog2FCCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigGenesBarplot_TOP50_Log2FC_", cmdl_input[7], "_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20, height = barplot_height)
plot(genesBarplotLog2FCCondBUp)
dev.off()


plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = .25*lengthgenes , height = 15)
plot(isoformsBarplot)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_",cmdl_input[6],"_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_",cmdl_input[7],"_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotCondBUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_Log2FC_",cmdl_input[6],"_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotLog2FCCondAUp)
dev.off()
plotname = paste(cmdl_input[10], "/", prefix, "_SigIsoformsBarplot_TOP50_Log2FC_",cmdl_input[7],"_UP_adjPvalueLT05_Log2FcGT1.pdf", sep = "")
pdf(file = plotname, width = 20 , height = barplot_height)
plot(isoformsBarplotLog2FCCondBUp)
dev.off()




#Save Cuffdiff into a file

file_runInfo = sprintf("%s/%s",cmdl_input[11],"ReplicatesInfo.txt")
file_runInfo

file_line = "################################################################################\n"
file_head = sprintf("%s%s%s%s\n","Replicates information for ", cmdl_input[6], " vs. ", cmdl_input[7])

sink(file = file_runInfo, append = FALSE, type = "output", split = TRUE)
cat(file_line)
cat(file_line)
cat("Replicates:\n")
print(replicates(cuffdata))
cat("\n\n")
cat(file_line)
sink()


