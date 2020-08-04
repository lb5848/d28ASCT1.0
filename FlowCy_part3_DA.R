rm(list = ls())

# Load packages
library(rstudioapi)
library(devtools)
library("flowCore")
library("flowWorkspace")
library(cytofCore)
library(FlowSOM)
# library(cytofkit)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(flowViz)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(Rphenograph)
library(diffcyt)
library(SummarizedExperiment)
library(stringr)
library(ggcyto)
library(SingleCellExperiment)
library(scran)
library(scater)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory
# Define workingDirectory
wdName <- "Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")

setwd(workingDirectory)

# Load workspace and SCEobject
load("workspaceSCEclusterDR.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "mean")

plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")

plotExprHeatmap(sce, features = type_markers(sce), k = "meta12", by = "cluster_id", bars = TRUE, perc = TRUE,
                scale = "last")



FDR_cutoff <- 0.05
ei <- sce@metadata$experiment_info

# DA using GLMM

(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition",
                              cols_random = c("patient_id")))
contrast2 <- createContrast(c(0,0,1,1,1,1))
da_res2 <- diffcyt(sce,
                   formula = da_formula1, contrast = contrast2,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "meta12", verbose = TRUE, subsampling = TRUE,
                   transform = FALSE, normalize = FALSE) 
da2 <- rowData(da_res2$res)

plotDiffHeatmap(sce, da2, top_n = 12, all = TRUE, fdr = FDR_cutoff)

p <- plotMedExprs(sce, facet = "antigen", k = "meta8", features = "type")
p$facet$params$ncol <- 4
p


# Save current workspace
save(list = ls(), file = "workspaceSCEDA.rds")
