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
load("workspaceSCE.rds")
sce
CATALYST::pbMDS(sce, by = "sample_id", color_by = "condition", features = type_markers(sce), fun = "mean")
# Run FlowSOM and ConsensusClusterPlus
seed <- 123456
set.seed(seed)
sce <- cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 20, 
               verbose = TRUE, seed = seed)
delta_area(sce)
# Run dimensionality reduction
n_cells <- 3000
max_cells <- 1000
n_cells < n_events
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
start_runDR <- Sys.time()
sce <- runDR(sce, dr = "TSNE", cells = max_cells, features = "type", theta = 0.5, max_iter = 1000, 
             distMethod = "euclidean",
             PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = max_cells, features = "type")
sce <- runDR(sce, dr = "DiffusionMap", cells = max_cells, features = "type", assay = "exprs")
end_runDR <- Sys.time()
timestamp()
difftime(end_runDR, start_runDR, "secs")

save(list = ls(), file = "workspaceSCEclusterDR.rds")

# Plots
display.brewer.all(colorblindFriendly = TRUE)
delta_area(sce)
cdx <- type_markers(sce)


plotAbundances(sce, k = "meta20", by = "cluster_id", group_by = "condition")
plotExprHeatmap(sce, features = type_markers(sce), k = "meta20", by = "cluster_id",  fun = "mean", scale = "last",
                bars = TRUE, perc = TRUE)

plotAbundances(sce, k = "meta12", by = "cluster_id", group_by = "condition")
plotExprHeatmap(sce, features = type_markers(sce), k = "meta12", by = "cluster_id",  fun = "mean", scale = "last",
                bars = TRUE, perc = TRUE)

# UMAP plot color_by = "meta8", facet_by = "condition"
plot <- plotDR(sce, dr = "UMAP", color_by = "meta12", facet_by = "condition") + 
  geom_density_2d(binwidth = 0.006, colour = "black")
plot$facet$params$ncol <- 2
plot

plot <- plotDR(sce, dr = "DiffusionMap", color_by = "meta12", facet_by = "condition")
plot$facet$params$ncol <- 2
plot


markers <- c("TCF1", "CD69", "CD127")
plot1 <- plotDR(sce, dr = "DiffusionMap", color_by = markers)
plot2 <- plotDR(sce, dr = "DiffusionMap", color_by = "condition")

plotDR(sce, dr = "DiffusionMap", color_by = "condition")
markers <- c("TCF1", "CD69", "CD127", "Tbet", "Ki67", "PD1")
plotDR(sce, dr = "DiffusionMap", color_by = markers)
plotDR(sce, dr = "UMAP", color_by = markers)
plotDR(sce, dr = "TSNE", color_by = markers)
# UMAP plot color_by = "meta8", facet_by = "sample_id"
plot <- plotDR(sce, dr = "UMAP", color_by = "meta12", facet_by = "sample_id")
plot$facet$params$ncol <- 3
plot

# UMAP color_by CD27 DNAM1
plotDR(sce, dr = "UMAP", color_by = c("CD27", "DNAM1"), facet_by = "condition")

plotDR(sce, dr = "UMAP", color_by = "meta8", facet_by = "sample_id")
