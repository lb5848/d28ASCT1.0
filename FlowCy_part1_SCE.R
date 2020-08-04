rm(list = ls())

# Install required packages
if(!require('rstudioapi')) {
  install.packages('rstudioapi')
}
if(!require('devtools')){
  install.packages("devtools")
}
if(!require('flowCore')){
  BiocManager::install("flowCore")
}
if(!require('cytofCore')){
  devtools::install_github("nolanlab/cytofCore")
}
# if(!require('JinmiaoChenLab/cytofkit')){
#   #do not update hexbin
#   remotes::install_github("JinmiaoChenLab/cytofkit")
# }
# if(!require("CytoML")){
#   BiocManager::install("CytoML")
# }
if(!require('FlowSOM')){
  BiocManager::install("FlowSOM")
}
if(!require('cluster')){
  install.packages("cluster")
}
if(!require('Rtsne')){
  install.packages("Rtsne")
}
if(!require('ggplot2')){
  install.packages("ggplot2")
}
if(!require('dplyr')){
  install.packages("dplyr")
}
if(!require('ggthemes')){
  install.packages("ggthemes")
}
if(!require('RColorBrewer')){
  install.packages('RColorBrewer')
}
if(!require("uwot")){
  install.packages("uwot")
}
if(!require("CATALYST"))
  BiocManager::install("CATALYST")
if(!require("diffcyt"))
  BiocManager::install("diffcyt")
if(!require("stringr"))
  BiocManager::install("stringr")

if(!require("JinmiaoChenLab/Rphenograph")){
  remotes::install_github("JinmiaoChenLab/Rphenograph")
}
BiocManager::install("Rphenograph")
if(!require("scran"))
  BiocManager::install("scran")
if(!require("scater"))
  BiocManager::install("scater")
if(!require("ggcyto"))
  BiocManager::install("ggcyto")
if(!require("SingleCellExperiment"))
  BiocManager::install("SingleCellExperiment")
if(!require("Rphenograph"))
  BiocManager::install("Rphenograph")
if(!require("flowWorkspace"))
  BiocManager::install("flowWorkspace")
if(!require("flowVS"))
  install.packages(file.choose(), repos = NULL, type = "source")
if(!require("flowStats")){
  BiocManager::install("flowStats")
}

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
library(flowVS)
library(readxl)
library(flowStats)
library(FlowSOMworkshop)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory


# Define fcs directory
fcsName <- "fcs_files"
fcsDirectory <- paste(PrimaryDirectory, fcsName, sep = "/")
dir.create(fcsDirectory)


# Define workingDirectory
wdName <- "Working_DirectoryFCS"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

# List files
FCSfiles <- list.files(fcsDirectory, pattern = ".fcs$", full = FALSE)

# Sample datasheet
sample_mdName <- "sample_md.xls"
sample_md <- read_excel(path = paste(PrimaryDirectory, sample_mdName, sep = "/"), col_names = TRUE) 
sample_md <- data.frame(sample_md)


# Define groups and panel marker_classes
marker <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-W", "Event #","TIGIT FITC", "Live.Dead", "TCF1 PE",
            "PD1 PE-Dazzle 594", "CD8 PerCP-Cy5.5", "Ki67 PE-Cy7", "Tbet APC",
            "CD69 APC-R700", "CD3 APC-Fire 750", "CD127 BV421", "TIME")

# Prepare panel and metadata for SCE object
setwd(workingDirectory)

# Create flowSet from FCSfiles
flowSet <- read.flowSet(files = FCSfiles, path = fcsDirectory, truncate_max_range = FALSE)
colnames(flowSet)
colnames(flowSet) <- marker

# Introduce keyword $CYT = "FACS"
flowSet[[1]]@description$`$CYT`
flowSet[[1]]@description$`$CYT` <- "FACS"

# Build panel data.frame
fcs_colname <- colnames(flowSet)
marker_class <- c(rep("state", 6), "type", "state", "type", "type", "state", "type",
                  "type", "type", "state", "type", "state")

length(marker_class) == length(fcs_colname)
# antigen: marker_name w/o fluorochrome
antigen <- fcs_colname
antigen[7] <- "TIGIT"
antigen[9] <- "TCF1"
antigen[10] <- "PD1"
antigen[11] <- "CD8"
antigen[12] <- "Ki67"
antigen[13] <- "Tbet"
antigen[14] <- "CD69"
antigen[15] <- "CD3"
antigen[16] <- "CD127"
panel <- cbind(fcs_colname, antigen, marker_class)
panel <- as.data.frame(panel)
all(panel$fcs_colname %in% colnames(flowSet))
panel
# Build md data.frame (metadata)
file_name <- sample_md$file_name
tissue <- sample_md$tissue
condition <- paste(sample_md$condition, sample_md$tissue, sep = "_")
patient_id <- sample_md$patient_id
unique_identifier <- paste(condition, patient_id, sep = "_")
sample_id <- unique_identifier

md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)
md
md$condition <- factor(md$condition)
md$sample_id <- factor(md$sample_id)

write.csv(md, file = "metadata.csv", row.names = FALSE)
write.csv(panel, file = "panel.csv", row.names = FALSE)

# apply compensation
ff <- flowSet[[1]]
head(ff@exprs)
head(ff@description)
ff@description$`$SPILLOVER`

comp <- ff@description$`$SPILLOVER`
colnames(comp) <- marker[7:16]
colnames(comp)

for(i in 1:length(flowSet)){
  tmp <- flowSet[[i]]
  tmp <- compensate(tmp, comp)
  flowSet[[i]] <- tmp
}

# asinh transfomr
cofactor <- c(rep(1,6), 1000, 2512, 3980, 5012, 700, 1000, 7943, 1600, 250, 1000, 1)
names(cofactor) <- marker
cf <- cofactor[7:16]
names(cf) <- names(cofactor)[7:16]

asinhFlowSet <- transFlowVS(flowSet, channels = names(cf), cofactors = cf)
chs_of_interest <- colnames(flowSet)[7:16]
plot_aggregate(asinhFlowSet, channels = chs_of_interest, output_image = "FCSpreNormAsinh.png")
normFlowSet <- warpSet(asinhFlowSet, stains = colnames(asinhFlowSet)[c(8, 11, 15, 16)])
plot_aggregate(normFlowSet, channels = chs_of_interest, output_image = "FCSpostNormAsinh.png")

normFlowSet[[1]]@description$`$CYT`
normFlowSet[[1]]@description$`$CYT` <- "FACS"


# names(cfEstimate) <- marker
sce <- CATALYST::prepData(normFlowSet, panel = panel, md = md, transform = FALSE)
assay(sce, "exprs") <- assay(sce, "counts")
sce_tmp <- prepData(flowSet, panel = panel, md = md, transform = FALSE)
assay(sce, "counts") <- assay(sce_tmp, "counts")
# QC
p <- plotExprs(sce, features = NULL, color_by = "condition")
p$facet$params$ncol <- 4
p
n_events <- min(n_cells(sce))
n_events
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

CATALYST::pbMDS(sce, color_by = "condition", features = type_markers(sce), fun = "mean")

plotNRS(sce, features = type_markers(sce), color_by = "condition")


rm(sce_tmp)
# save workspace 
save(list = ls(), file = "workspaceSCE.rds")


