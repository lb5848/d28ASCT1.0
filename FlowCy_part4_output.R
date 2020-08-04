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
load("workspaceSCEDA.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "mean")

# Save .fcs per cluster
outputDirectory <- getwd()
outputDirectory <- paste(outputDirectory, "output", sep = "/")
dir.create(outputDirectory)
setwd(outputDirectory)
sce$Clusters <- cluster_ids(sce, "meta12")
# keep_dr = TRUE not all cells have DR
flowSet <- sce2fcs(sce, split_by = "Clusters", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(flowSet, outdir = outputDirectory, filename = "bycluster")
merged <- sce2fcs(sce, split_by = NULL, keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.FCS(merged, filename = "merged.fcs")
by_condition <- sce2fcs(sce, split_by = "condition", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(by_condition, outdir = outputDirectory, filename = "bycondition")

# Save final workspace
save(list = ls(), file = "workspaceFinal.rds")