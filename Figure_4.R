# # install spata 2.0
# install.packages("devtools")
# library(BiocManager)
# 
# if (!base::requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# 
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))
# 
# install.packages("Seurat")
# 
# # install tensorflow
devtools::install_github(repo = "kueckelj/confuns", force = TRUE)
devtools::install_github(repo = "theMILOlab/SPATAData", force = TRUE)
devtools::install_github(repo = "theMILOlab/SPATA2")
BiocManager::install("glmGamPoi")
# 
# # install.packages("remotes")
# remotes::install_github("rstudio/tensorflow")
# reticulate::install_python()
# 
# library(tensorflow)
# install_tensorflow(envname = "r-tensorflow")
# 
# install.packages("keras")
# library(keras)
# install_keras()
# 
# library(tensorflow)

# # if you want to use monocle3 related wrappers 
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")

# spata v 3.1.2
install.packages('/Users/......../SPATA2-3.1.2.tar.gz', repos = NULL, type = "source")

# load library
library(SPATA2)
library(openxlsx)
library(ggplot2)
library(devtools)
#library(monocle3)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(infercnv)
library(Seurat)
#library(SeuratData)
library(patchwork)
library(dplyr)
#library(dplyr)


# install BiocManager::install('glmGamPoi')
#BiocManager::install("glmGamPoi")

# ############################ SEURAT OBJECT ####################################
# 
# # for S1 
# #seurat from 10x
# data_dir <- '/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/Spatial_data_Nov_2024/processed/AS-1360737-LR-78146_S1/outs'
# list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
# S1 <- Load10X_Spatial(data.dir = data_dir,  filename = "filtered_feature_bc_matrix.h5")
# 
# #Data preprocessing
# plot1 <- VlnPlot(S1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# plot2 <- SpatialFeaturePlot(S1, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)
# 
# # scTransform
# S1 <- SCTransform(S1, assay = "Spatial", verbose = FALSE)
# 
# # rerun normalization to store sctransform residuals for all genes
# S1 <- SCTransform(S1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# # also run standard log normalization for comparison
# S1 <- NormalizeData(S1, verbose = FALSE, assay = "Spatial")
# # plot
# SpatialFeaturePlot(S1, features = c("CCL2", "CCL4", "IL1B"))
# 
# # assay
# Assays(S1)
# DefaultAssay(S1) 
# 
# # set active assay
# DefaultAssay(S1) <- "SCT"

# 
# ##################################### START OF SPATA 2.0 V3.1.2 ################################################################
# 
# # load spata object from 10x directory
# directory_10X <- '/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/Spatial_data_Nov_2024/processed/AS-1360737-LR-78146_S1/outs'
# 
# S1 <- 
#   initiateSpataObjectVisium(
#     sample_name = "10xmel", 
#     directory_visium = directory_10X # adjust to your liking 
#   )
# 
# # rename to spata_obj
# spata_obj <- S1
# 
# # all expression matrices before denoising
# activeMatrix(object = spata_obj)
# 
# # image handeling
# 
# # load SPATA2 inbuilt example data
# spata_obj <- setDefault(spata_obj, display_image = TRUE, pt_size = 1.5)
# 
# # the image 
# getImage(spata_obj)
# 
# # image dimensions in width, height and colors
# getImageDims(spata_obj)
# ## [1] 576 600   3
# 
# # image range in terms of data coordinates
# getImageRange(spata_obj)
# 
# # by default, the active image is chosen
# #Image processing
# spata_obj <- identifyPixelContent(spata_obj)
# plotImageMask(spata_obj)
# plotPixelContent(spata_obj)
# plotImage(spata_obj)
# plotImage(spata_obj, outline = TRUE, line_size = 1)
# 
# # this is the default input for the visium platform and has already been 
# # called in initiateSpataObjectVisium(). 
# # if the results do not satisfy you, you can run it over and over again with 
# # different parameter inputs 
# spata_obj <- identifyTissueOutline(spata_obj, method = "obs", eps = "125um", minPts = 3)
# plotSurface(spata_obj, color_by = "tissue_section", pt_clrp = "tab20")
# 
# #Tissue outline parameters
# spata_obj <- identifyTissueOutline(spata_obj, eps = "125um", minPts = 3)
# plotSurface(spata_obj)
# plotSurface(spata_obj, color_by = "tissue_section")
# 
# #Spatial outliers
# # uses the results of identifyTissueOutline() to create a logical variable called sp_outlier
# spata_obj <- identifySpatialOutliers(spata_obj, method = "obs")
# 
# # plot_with_outliers <- plotSurface(object, color_by = "sp_outlier", clrp_adjust = c("TRUE" = "blue"))
# # 
# # # remove where sp_outlier == TRUE
# # object <- removeSpatialOutliers(object)
# # 
# # plot_without_outliers <- plotSurface(object, color_by = "sp_outlier")
# # 
# # # left plot
# # plot_with_outliers
# # 
# # # right plot
# # plot_without_outliers
# 
# # Data processing
# # before
# nGenes(spata_obj)
# ## [1] 33538
# 
# # removes stress genes
# object <- removeGenesStress(spata_obj)
# 
# # removes genes that were not detected in any of the observations
# object <- removeGenesZeroCounts(object)
# 
# # afterwards
# nGenes(object)
# ## [1] 21556
# 
# # spatially variable genes
# spata_obj <- runSPARKX(spata_obj, verbose = FALSE)
# # get genes with a p-value < 0.01
# sparkx_genes <- getSparkxGenes(spata_obj, threshold_pval = 0.01)
# str(sparkx_genes)
# 
# # visualize in space
# plotSurfaceComparison(spata_obj, color_by = head(sparkx_genes, 6), nrow = 2)

######################################### LOAD SEURAT OBJECT  ############################################################

# # laod spata object 
# S1 <- loadSpataObject("/Users/sajibchakraborty/Downloads/AS-1153289-LR-70634_S1.rds")
# S2 <- loadSpataObject("/Users/sajibchakraborty/Downloads/AS-1153291-LR-70634_S2.rds")
# S3 <- loadSpataObject("/Users/sajibchakraborty/Downloads/AS-1153293-LR-70634_S3.rds")
# S4 <- loadSpataObject("/Users/sajibchakraborty/Downloads/AS-1153295-LR-70634_S4.rds")
# mel10x <- loadSpataObject("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_Spatial_Trans and APOE variant/10X/melSpataObj_Ucell_ALL.rds")
# #old <- loadSpataObject("/Users/sajibchakraborty/Documents/GBM_5ALA_Astrocyte/SPATIAL_Seg_MLA_DM_MALIG/MALDI_IMC/Tumor/248_T_SPATA_CNV_Pred.RDS")
# 
# spata_object <- mel10x
# 
# plotSurfaceInteractive(object = spata_object)

# load seurat_object list from DT
seurat_list <- readRDS('/Users/s......../proximap_runs.rds')

# proximap all samples 

proximap_list <- readRDS('/Users/......./melanoma_runs.rds')

# # from DT
# dataDir <- "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/Seurat_Object"
# fileFh <- list.files(dataDir, pattern = ".rds", full.names = TRUE, recursive = TRUE)
# outDir <- file.path(dataDir, "spata_files")
# dir.create(outDir, recursive = TRUE)
# 
# 
# for (i in fileFh) {
# 
#   # ##################################################
#   # Get the sample name
#   sample_name <- gsub(".rds", "", basename(i))
#   
#   # ##################################################
#   # Load Seurat object
#   # seurat_object <- Load10X_Spatial(data.dir = i)
#   seurat_object <- readRDS(i)
#   seurat_object <- NormalizeData(seurat_object)
#   all.genes <- rownames(seurat_object)
#   seurat_object <- ScaleData(seurat_object, features = all.genes)
# 
#   # ##################################################
#   # Fix the layers of the Seurat object
#   # object <- seurat_object
#   # object[["RNA3"]] <- as(object = object[["Spatial"]], Class = "Assay")
#   # DefaultAssay(object) <- "RNA3"
#   # object[["Spatial"]] <- NULL
#   # object <- RenameAssays(object = object, RNA3 = 'Spatial')
#   
#   #break
#   
#   # ##################################################
#   # Convert the Seurat to SPATA
#   spata_object <- SPATA2::asSPATA2(
#     object = seurat_object,
#     assay_name = "Spatial",
#     sample_name = "S1",
#     image_name = "slice1", 
#     spatial_method = "Visium"
#   )
#   # break
#   # ##################################################
#   # Save the SPATA object
#   saveRDS(spata_object, file.path(outDir, paste0(sample_name, ".rds")))
#   
# }

# extract proximap data frames
S1_df <- seurat_list[["AS-1360737-LR-78146_S1"]]@misc[["PROXIMAP"]]@colocalization@proximap_out[["gg_df"]]
S1_gsVal <- seurat_list[["AS-1360737-LR-78146_S1"]]@misc[["PROXIMAP"]]@colocalization@geneset_validity_estimation[["gs_thresholds"]]
S2_df <- seurat_list[["AS-1360739-LR-78146_S2"]]@misc[["PROXIMAP"]]@colocalization@proximap_out[["gg_df"]]
S2_gsVal <- seurat_list[["AS-1360739-LR-78146_S2"]]@misc[["PROXIMAP"]]@colocalization@geneset_validity_estimation[["gs_thresholds"]]
S3_df <- seurat_list[["AS-1360741-LR-78146_S3"]]@misc[["PROXIMAP"]]@colocalization@proximap_out[["gg_df"]]
S3_gsVal <- seurat_list[["AS-1360741-LR-78146_S3"]]@misc[["PROXIMAP"]]@colocalization@geneset_validity_estimation[["gs_thresholds"]]

#save files
# Load necessary library
library(openxlsx)

# Define the save path
save_path <- "/Users/....../PROXIMAP_NEW_SAMPLES"

# List of data frames and their corresponding file names
data_frames <- list(
  S1 = list(df = S1_df, gsVal = S1_gsVal),
  S2 = list(df = S2_df, gsVal = S2_gsVal),
  S3 = list(df = S3_df, gsVal = S3_gsVal)
)

# Loop through each dataset and write to Excel
for (name in names(data_frames)) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add sheets to the workbook
  addWorksheet(wb, "gg_df")
  writeData(wb, "gg_df", data_frames[[name]]$df)
  
  addWorksheet(wb, "gs_thresholds")
  writeData(wb, "gs_thresholds", data_frames[[name]]$gsVal)
  
  # Save the workbook
  file_path <- file.path(save_path, paste0(name, "_PROXIMAP.xlsx"))
  saveWorkbook(wb, file_path, overwrite = TRUE)
  
  # Print confirmation
  cat("Saved:", file_path, "\n")
}


######################################### CONVERT SEURAT OBJECTS TO SPATA OBJECTS #######################################################################
# Assuming your seurat_list is already loaded
# Initialize an empty list to store SPATA objects
spata_list <- list()

# Loop through the Seurat list
for (name in names(seurat_list)) {
  # Extract the Seurat object
  seurat_object <- seurat_list[[name]]
  
  # Extract the sample suffix (e.g., S1, S2, etc.)
  sample_name <- sub(".*_", "", name)
  
  # Convert to SPATA object
  spata_object <- asSPATA2(
    object = seurat_object,
    sample_name = sample_name,
    platform = "VisiumSmall",
    img_scale_fct = "lowres",
    assay_name = "Spatial",
    assay_modality = "gene"
  )
  
  # Store in the list with the new name
  spata_list[[sample_name]] <- spata_object
}

# Check SPATA objects
for (name in names(spata_list)) {
  cat("SPATA Object:", name, "\n")
  show(spata_list[[name]])
}

# Define the save path
save_path <- "/Users/....../SPATA_Object"

# Loop through the SPATA objects in the spata_list
for (name in names(spata_list)) {
  # Create the full file path
  file_path <- file.path(save_path, paste0(name, "_SPATA_object.RDS"))
  
  # Save the SPATA object
  saveRDS(spata_list[[name]], file = file_path)
  
  # Print confirmation
  cat("Saved:", file_path, "\n")
}

##################################### START OF DOWNSTREAM ANALYSIS WITH SPATA 2.0 V3.1.2 ################################################################

# now load spata objects 

S1 <- loadSpataObject('/Users/......../S1.RDS')
#S2 <- loadSpataObject('/Users/......../S2.RDS')
#S3 <- loadSpataObject('/Users/......../S3.RDS')

#S4 <- loadSpataObject('/Users/......../S4.RDS')
#S5 <- loadSpataObject('/Users/......../S5.RDS')
#S6 <- loadSpataObject('/Users/......../S6.RDS')


################################################## Matrix processing alternative to denoising  ############################

# convert the name of spata object 
spata_obj <- S1

# interactive plot
plotSurfaceInteractive(object = spata_obj)

# obtain matrix names prior to normalization
getMatrixNames(spata_obj)
## [1] "counts"

plot_before <- 
  plotSurface(spata_obj, color_by = "CCL2") + labs(color = "MAG\n(Counts)")

# create log normalized matrix
spata_obj <- normalizeCounts(spata_obj, method = "CLR", overwrite = TRUE)
## Normalizing layer: counts
## 01:10:34 Active matrix in assay 'gene': 'LogNormalize'

# obtain matrix names after normalization
getMatrixNames(spata_obj)
## [1] "counts"       "LogNormalize"

# check active matrix 
activeMatrix(spata_obj)
## [1] "LogNormalize"

plot_afterwards <- 
  plotSurface(spata_obj, color_by = "GZMB") + labs(color = "MAG\n(logNorm)")

# left plot
plot_before

# right plot
plot_afterwards


############################# OBJECT CTEATIONA, DENOISING AND SAVE COMPLETE ###############################################

############################### CLUSTERING ################################################################################

# S6 
plotImage(object = spata_obj)

# # run the pipeline
spata_obj <- 
  runBayesSpaceClustering(
    object = spata_obj, 
    name = "bayes_space", # the name of the output grouping variable
    overwrite = TRUE
  )
# 
# # results are immediately stored in the objects feature data
getGroupingOptions(spata_obj)
getGroupNames(spata_obj)
# run new clustering
# run PCA based on which clustering is conducted
spata_obj <- runPCA(spata_obj, n_pcs = 20)

#alternative clustering
spata_obj <- 
  runKmeansClustering(
    object = spata_obj, 
    ks = c(7, 8), 
    methods_kmeans = "Lloyd"
  )

# right plot
plotSurface(
  object = spata_obj, 
  color_by = "Lloyd_k7",
  pt_clrp = "jco"
)


########################### END OF CLUSTERING ######################################

################################## DEA across clusters #############################################################

# check grouping option
getGroupingOptions(object = spata_obj)

#Running the analysis
spata_obj <- runDEA(object = spata_obj, across = "bayes_space", method_de = "wilcox")


#check results 
getDeaResultsDf(spata_obj, across = "bayes_space")

# Extracting results

# extract the complete data.frame

#spata_obj
spata_obj_DEA <- getDeaResultsDf(
  object = spata_obj, 
  across = "bayes_space", 
  method_de = "wilcox",
  n_highest_lfc = 500, # top 500 genes
  max_adj_pval = 0.01
)


######################################  VISALIZATION OF DEA ##################################################################

# Heatmap
# Check the expression matrix and define unique breaks
expr_mtr <- getMatrix(spata_obj)
summary(expr_mtr)

# Create a unique sequence for breaks
breaks_input <- seq(min(expr_mtr), max(expr_mtr), length.out = 101)

# Define a color palette
colors <- NULL

####################################### heatmap ############################################################################

# Call the plotDeaHeatmap function
hm <- plotDeaHeatmap(
  object = spata_obj,
  across = "bayes_space",
  method_de = "wilcox",
  n_highest_lfc = 5,
  n_bcs = 100,
  breaks = breaks_input,
)

ggsave(plot = hm, file=file.path('/Users...../hm.pdf'),
       width = 5, height = 5)

####################################### dotplot ############################################################################

dp <- plotDeaDotPlot(
  object = spata_obj, 
  across = "bayes_space",
  n_highest_lfc = 5,
  by_group = TRUE, 
  scales = "free_y", 
  nrow = 1
)

################################## DEA Plot #########################################

DEAplot <- plotDeaDotPlot(
  object = spata_obj, 
  across = "bayes_space",
  color_by = "bayes_space",
  pt_clrp = "npg",
  size_by = "avg_log2FC",
  n_highest_lfc = 5, 
  by_group = FALSE
)

############################### END OF DEA #############################################################

############################### START OF LIGAND RECEPTOR DATABASE MAPPING ##############################
#INSTALL omnipath
install_github('saezlab/OmnipathR')

# laod
library(OmnipathR)
## We check some of the different interaction databases
get_interaction_resources()
## We check some of the different intercell categories
get_intercell_generic_categories()
## We import the intercell data into a dataframe
intercell <- import_omnipath_intercell(
  scope = "generic",
  aspect = "locational"
)

## We check the intercell annotations for the individual components of
## our previous complex. We filter our data to print it in a good format
DEA_mapping <- dplyr::filter(intercell, genesymbol %in% spata_obj_DEA$gene) %>%
  dplyr::distinct(genesymbol, parent, .keep_all = TRUE) %>%
  dplyr::select(category, genesymbol, parent) %>%
  dplyr::arrange(genesymbol)

#high conf 
icn <- import_intercell_network(high_confidence = TRUE)

## We check the intercell annotations for the individual components of
## our previous complex. We filter our data to print it in a good format
S2_DEA_mapping <- dplyr::filter(icn, icn$source_genesymbol %in% spata_obj_DEA$gene) %>%
  dplyr::distinct(icn$source_genesymbol, parent, .keep_all = TRUE) %>%
  dplyr::arrange(icn$source_genesymbol)

############################### START OF HYPERG ########################################################

#load library
library(hypeR)
library(openxlsx)

# convert cluster markers dataframe as list 
#Spata
#S1_DEA <- read.xlsx("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S1/S1_DEA.xlsx")
list_spata <- split(unlist(Spata_DEA$gene), Spata_DEA$bayes_space)
names(list_spata) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10")


# msig DB all
C1 <- msigdb_gsets("Homo sapiens", "C1", clean=TRUE)
C2_KEGG <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
C2_Biocarta <- msigdb_gsets("Homo sapiens", "C2", "CP:BIOCARTA", clean=TRUE)
C2_REACTOME <- msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean=TRUE)
C2_WIKIPATHWAYS <- msigdb_gsets("Homo sapiens", "C2", "CP:WIKIPATHWAYS", clean=TRUE)
HALLMARK <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
C5_GOBP <- msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
C5_GOCC <- msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
C5_GOMF <- msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
C5_HPO <- msigdb_gsets("Homo sapiens", "C5", "HPO", clean=TRUE)
C6 <- msigdb_gsets("Homo sapiens", "C6", clean=TRUE)
C7 <- msigdb_gsets("Homo sapiens", "C7", "IMMUNESIGDB", clean=TRUE)
C8 <- msigdb_gsets("Homo sapiens", "C8", clean=TRUE)

# now define list_cluster_markers
#list_cluster_markers <- list_S1
#list_cluster_markers <- list_S2
#list_cluster_markers <- list_S3
#list_cluster_markers <- list_S4
list_cluster_markers <- list_spata

# RUN HYPER G for cluster markers

C1_hyp_cluster <- hypeR(list_cluster_markers, C1, test="hypergeometric", background=7000)
C2_KEGG_hyp_cluster <- hypeR(list_cluster_markers, C2_KEGG, test="hypergeometric", background=7000)
C2_Biocarta_Hyp_cluster <- hypeR(list_cluster_markers, C2_Biocarta, test="hypergeometric", background=7000)
C2_REACTOME_Hyp_cluster <- hypeR(list_cluster_markers, C2_REACTOME, test="hypergeometric", background=7000)
C2_WIKI_Hyp_cluster <- hypeR(list_cluster_markers, C2_WIKIPATHWAYS, test="hypergeometric", background=7000)
HALLMARK_hyp_cluster <- hypeR(list_cluster_markers, HALLMARK, test="hypergeometric", background=7000)
GOBP_hyp_cluster <- hypeR(list_cluster_markers, C5_GOBP, test="hypergeometric", background=7000)
GOMF_hyp_cluster <- hypeR(list_cluster_markers, C5_GOMF, test="hypergeometric", background=7000)
GOCC_hyp_cluster <- hypeR(list_cluster_markers, C5_GOCC, test="hypergeometric", background=7000)
HPO_hyp_cluster <- hypeR(list_cluster_markers, C5_HPO, test="hypergeometric", background=7000)
C6_hyp_cluster <- hypeR(list_cluster_markers, C6, test="hypergeometric", background=7000)
C7_hyp_cluster <- hypeR(list_cluster_markers, C7, test="hypergeometric", background=7000)
C8_hyp_cluster <- hypeR(list_cluster_markers, C8, test="hypergeometric", background=7000)


# plots for cell and clsuter # change directory to save different plots

C1_plot <- hyp_dots(C1_hyp_cluster, merge=TRUE, pval=0.05, title="C1")
C2_KEGG_plot <- hyp_dots(C2_KEGG_hyp_cluster, merge=TRUE, pval=0.05, title="KEGG")
C2_BIOCARTA_plot <- hyp_dots(C2_Biocarta_Hyp_cluster, merge=TRUE, pval=0.05, title="Biocarta")
C2_REACTOME_plot <- hyp_dots(C2_REACTOME_Hyp_cluster, merge=TRUE, pval=0.05, title="REACTOME")
C2_wiki_plot <- hyp_dots(C2_WIKI_Hyp_cluster, merge=TRUE, pval=0.05, title="WIKIPATHWAY")
Hallmark_plot <- hyp_dots(HALLMARK_hyp_cluster, merge=TRUE, pval=0.05, title="HALLMARK")
GOBP_plot <- hyp_dots(GOBP_hyp_cluster, merge=TRUE, pval=0.05, title="GOBP")
GOMF_plot <- hyp_dots(GOMF_hyp_cluster, merge=TRUE, pval=0.05, title="GOMF")
GOCC_plot <- hyp_dots(GOCC_hyp_cluster, merge=TRUE, pval=0.05, title="GOCC")
HPO_plot <- hyp_dots(HPO_hyp_cluster, merge=TRUE, pval=0.05, title="HPO")
C6_plot <- hyp_dots(C6_hyp_cluster, merge=TRUE, pval=0.05, title="C6")
C7_plot <- hyp_dots(C7_hyp_cluster, merge=TRUE, pval=0.05, title="C7")
C8_plot <- hyp_dots(C8_hyp_cluster, merge=TRUE, pval=0.05, title="C8")


#set dir
setwd("/Users/....../HyperG")
#setwd("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S2/HyperG/Top500")
#setwd("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S3/HyperG/Top500")
#setwd("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S4/HyperG/Top500")

#save plots 
ggsave(C1_plot, filename = "C1_plot.pdf",width = 8, height = 5)
ggsave(C2_KEGG_plot, filename = "KEGG_plot.pdf",width = 8, height = 5)
ggsave(C2_BIOCARTA_plot, filename = "BIOCARTA_plot.pdf",width = 8, height = 5)
ggsave(C2_REACTOME_plot, filename = "REACTOME.pdf",width = 8, height = 5)
ggsave(C2_wiki_plot, filename = "WIKIPATHWAYS.pdf",width = 8, height = 5)
ggsave(Hallmark_plot, filename = "HALLMARK.pdf",width = 8, height = 5)
ggsave(GOBP_plot, filename = "GOBP.pdf",width = 8, height = 5)
ggsave(GOMF_plot, filename = "GOMF.pdf",width = 8, height = 5)
ggsave(GOCC_plot, filename = "GOCC.pdf",width = 8, height = 5)
ggsave(HPO_plot, filename = "HPO.pdf",width = 8, height = 5)
ggsave(C6_plot, filename = "C6.pdf",width = 8, height = 5)
ggsave(C7_plot, filename = "C7.pdf",width = 8, height = 5)
ggsave(C8_plot, filename = "C8.pdf",width = 8, height = 5)

#Save results
hyp_to_excel(C1_hyp_cluster, file_path="C1.xlsx")
hyp_to_excel(C2_KEGG_hyp_cluster, file_path="KEGG.xlsx")
hyp_to_excel(C2_Biocarta_Hyp_cluster, file_path="BIOCARTA.xlsx")
hyp_to_excel(C2_REACTOME_Hyp_cluster, file_path="REACTOME.xlsx")
hyp_to_excel(C2_WIKI_Hyp_cluster, file_path="WIKI.xlsx")
hyp_to_excel(HALLMARK_hyp_cluster, file_path="HALLMARK.xlsx")
hyp_to_excel(GOBP_hyp_cluster, file_path="GOBP.xlsx")
hyp_to_excel(GOMF_hyp_cluster, file_path="GOMF.xlsx")
hyp_to_excel(GOCC_hyp_cluster, file_path="GOCC.xlsx")
hyp_to_excel(HPO_hyp_cluster, file_path="HPO.xlsx")
hyp_to_excel(C6_hyp_cluster, file_path="C6.xlsx")
hyp_to_excel(C7_hyp_cluster, file_path="C7.xlsx")
hyp_to_excel(C8_hyp_cluster, file_path="C8.xlsx")

###############################################################################################################

######################################  GENE SET ENRICHMENT AND SURFACE PLOT ##################################

###############################################################################################################

# GSEA
# DEFINE TAM SIGNATURES
geneMat_TAM <- read.xlsx(file.path("/Users/......./Refined_GeneSet.xlsx"))
geneList_TAM <- lapply(1:ncol(geneMat_TAM), function(i) as.character(geneMat_TAM[, i]))
names(geneList_TAM) <- colnames(geneMat_TAM)
geneList_TAM <- lapply(geneList_TAM, function(i) unique(i[!is.na(i)]))
# geneList$ALL <- unique(unlist(geneList))

# ############################
# # Basic extracting functions
# 
# if(FALSE){
#   # the essential data.frame
#   getSpataDf(spata_obj)
#   
#   # dimensional reduction data
#   getUmapDf(spata_obj)
#   
#   # barcode spot coordinates
#   getCoordsDf(spata_obj)
# }
# 
# coords_df <- getCoordsDf(spata_obj)
# coords_df
# 
# joinWith(object = spata_obj, 
#          spata_df = coords_df,
#          features = "seurat_clusters", # cluster belonging
#          verbose = FALSE)
# 
# # output 
# joined_df


# GENES
#spata.genes.meta <- getGeneMetaData(spata_obj)
spata.genes <- getGenes(spata_obj)
geneList_TAM <- lapply(geneList_TAM, intersect, y = spata.genes)


### Adding TAM genesets
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Memory_effector_1', genes = geneList_TAM$Memory_effector_1, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'LA_TAM', genes = geneList_TAM$LA_TAM, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Neural_crest_like', genes = geneList_TAM$Neural_crest_like, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Transitory', genes = geneList_TAM$Transitory, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Transitory_Melanocytic', genes = geneList_TAM$Transitory_Melanocytic, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Melanocytic', genes = geneList_TAM$Melanocytic, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Exhausion_cell_cycle_CD8', genes = geneList_TAM$Exhausion_cell_cyle_CD8, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Exhausted_CD8', genes = geneList_TAM$Exhausted_CD8, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'IFN_TAMs', genes = geneList_TAM$IFN_TAMs, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Memory_effector_2', genes = geneList_TAM$Memory_effector_2, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Activated_CD8', genes = geneList_TAM$Activated_CD8, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class = 'mygs', name = 'Treg', genes = geneList_TAM$Treg, overwrite = TRUE)


##############
# PLOT SURFACE
# get genesets
getGeneSetOverview(spata_obj)
all_signatures <- getSignatureNames(spata_obj, assay_name = "gene")
str(all_signatures)
str(hallmark_signatures)
str(mygs)
#plots$S4_Image

genesetS2 <- c("mygs_IFN_TAMs",
               "mygs_LA_TAM",
               "mygs_Neural_crest_like",
               "mygs_Transitory",
               "mygs_Transitory_Melanocytic",
               "mygs_Melanocytic",
               "mygs_Exhausion_cell_cycle_CD8",
               "mygs_Exhausted_CD8",
               "mygs_Memory_effector_1",
               "mygs_Memory_effector_2",
               "mygs_Activated_CD8",
               "mygs_Treg")


# check which is the active matrix
getMatrix(spata_obj)
getMatrixNames(spata_obj)
activeMatrix(spata_obj)

getMatrixNames(spata_obj)
# results [1] "counts"       "LogNormalize"

# normalized data with other algorihms
spata_obj <- normalizeCounts(
  spata_obj,
  method = "CLR",
  overwrite = TRUE)

# if needed set the active matrix
spata_obj <- setActiveMatrix(spata_obj, "CLR")
spata_obj <- setActiveMatrix(spata_obj, "SCT")
spata_obj <- setActiveMatrix(spata_obj, "LogNormalize")

# surface plot 

# open application to obtain a list of plots
#plots <- plotSurfaceInteractive(object = spata_obj)

# ssGSEA
Spata_GS <- plotSurfaceComparison(object = spata_obj, 
                                  color_by = genesetS2,
                                  method_gs = "zscore",
                                  smooth = TRUE,
                                  pt_clrsp = "inferno",
                                  smooth_span = 0.2,
                                  pt_size = 2,
                                  display_image = FALSE
)

Spata_GS

#Retrieve genes from SPATA object
spata.genes <- getGenes(spata_obj)

#Retrieve IFN gene set
IFN <- geneList_TAM$IFN_TAMs

# LATAM
LATAM <- geneMat_TAM$LA_TAM

#Retrieve IFN gene set
CD8ac <- geneList_TAM$Activated_CD8

#Retrieve IFN gene set
CD8ex <- geneList_TAM$Exhausted_CD8

#Intersect IFN gene set with SPATA genes
refined_IFN <- intersect(IFN, spata.genes)

#Intersect IFN gene set with SPATA genes
refined_LATAM <- intersect(LATAM, spata.genes)

#Intersect IFN gene set with SPATA genes
refined_CD8ac <- intersect(CD8ac, spata.genes)

#Intersect IFN gene set with SPATA genes
refined_CD8ex <- intersect(CD8ex, spata.genes)

S3_GOI_S3 <- plotSurfaceComparison(object = spata_obj, 
                                   color_by = GOI_S3,
                                   method_gs = "mean",
                                   smooth = TRUE,
                                   pt_clrsp = "inferno",
                                   smooth_span = 0.2,
                                   pt_size = 2,
                                   display_image = FALSE
)

S3_GOI_S3



#################################################### CYTOSPACE #################################################################

################################# EXTRACT EXPRESSION MATRIX FROM SPATA OBJECT ######################################################

spata_obj <- S6

S6mat <- getCountMatrix(spata_obj)
S6cor <- getCoordsDf(spata_obj)
S6cor <- S6cor[, -5,-6]
S6cor <- S6cor[, -2]
S6mat[1:5, 1:5]

S6df <- as.data.frame(S6mat)

d <- S6df
names <- rownames(d)
rownames(d) <- NULL
data <- cbind(names,d)
colnames(data)[1] <- "V1"
# Rename columns in the existing data frame
colnames(S6cor)[colnames(S6cor) == "x"] <- "row"
colnames(S6cor)[colnames(S6cor) == "y"] <- "col"

# Define the file path
file_path <- '/Users/..../S6mat.txt'

# Write data to a tab-delimited text file without row names
write.table(data, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

write.table(S6cor, file = "/Users/..../S6Cor.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# RUN cytospace web tool and get the outputs

# install UCELL
BiocManager::install("UCell")

# # load library
# -----
library(openxlsx)
library(ggplot2)
library(UCell)
library(gridExtra)
library(reshape2)
library(tidyr)


# -----
# load geneset
# -----
geneset_OF_interest <- read.xlsx("/Users/...../geneset.xlsx")
genesetList <- lapply(1:ncol(geneset_OF_interest), function(i) {geneset_OF_interest[, i]})
names(genesetList) <- colnames(geneset_OF_interest)
genesetList <- lapply(genesetList, function(i) {unique(i[!is.na(i)])})

# -----
# load single cell dataset
# ------- 

melanoma_sc <- read.delim("/Users/..../melanoma_scRNA_GEP.txt", row.names = 1)

# load coordinates

#S1 S2 S3 S4
#melanoma_sc_coordinate <- read.csv("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/Cytoscape_Melanoma/Mel_spatial10x_deconvolution_basedon_Tirosh/results/Cytospace_Output_Ucell/S1/results/assigned_locations.csv")
#melanoma_sc_coordinate <- read.csv("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/Cytoscape_Melanoma/Mel_spatial10x_deconvolution_basedon_Tirosh/results/Cytospace_Output_Ucell/S2/results/assigned_locations.csv")
#melanoma_sc_coordinate <- read.csv("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/Cytoscape_Melanoma/Mel_spatial10x_deconvolution_basedon_Tirosh/results/Cytospace_Output_Ucell/S3/results/assigned_locations.csv")
melanoma_sc_coordinate <- read.csv('/Users/...../assigned_locations.csv')


# -----
# get the cells to work
# -----
cells_to_work <- intersect(colnames(melanoma_sc), melanoma_sc_coordinate$OriginalCID)

# -----
# subset single cell dataset
# -----
melanoma_sc <- melanoma_sc[, cells_to_work]


# -----
# deduplicate the coordinates df
# -----
melanoma_sc_coordinate <- dplyr::distinct(melanoma_sc_coordinate, OriginalCID, .keep_all = TRUE)

# -----
# calculate ucell
# -----
melanoma_ucell <- ScoreSignatures_UCell(melanoma_sc, features = genesetList)
melanoma_ucell <- as.data.frame(melanoma_ucell)


# -----
# prepare the df for ggplot
# -----
gg_df <- melanoma_sc_coordinate
# gg_df$LA_tam_ucell <- melanoma_ucell$LA_TAM.Genes_UCell[match(rownames(melanoma_ucell), melanoma_sc_coordinate$OriginalCID)]

# -----
# add ucell scores to the gg dataframe
# -----
melanoma_ucell$cell_name <- rownames(melanoma_ucell)
rownames(melanoma_ucell) <- NULL
gg_df <- merge(gg_df, melanoma_ucell, by.x = "OriginalCID", by.y = "cell_name")
colnames(gg_df) <- gsub("_UCell", "", colnames(gg_df))


# -----
# identify unique spots
# -----
if (FALSE) {
  gg_df$unique_spot <- paste0(gg_df$x, "_", gg_df$y)
  ggplot(gg_df, aes(X, Y, color = unique_spot)) +
    geom_point() +
    theme_bw() + 
    theme(legend.position = "none")  # Remove the legend
}

# -----
# calculate distance
# -----
coordinate_dist <- dist(gg_df[, c("row", "col")], method = "euclidean")
distance_matrix <- as.matrix(coordinate_dist)
min(distance_matrix[distance_matrix != 0])

gg_df$row <- gg_df$row + runif(nrow(gg_df), min = -2, max = 2)
gg_df$col <- gg_df$col + runif(nrow(gg_df), min = -2, max = 2)

# get out dir

outdir <- "/Users/..../results"


# # -----
# # plot the enrichment
# # -----
# p <- ggplot(gg_df, aes(x, y, color = CellType)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   scale_colour_brewer(palette = "Paired") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# p
# 

####
library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)

# Define the order of CellType for plotting
plot_order <- c("Macrophages", "Melanoma cells", "Fibroblasts","CD8 T cells","CD4 T cells", "B cells",  "Endothelial cells")

# Define custom colors for each CellType
custom_colors <- c("Melanoma cells" = "#F6D776",
                   "Endothelial cells" = "#508D69",
                   "Fibroblasts" = "#FFC5C5",
                   "Macrophages" = "#FF004D",
                   "CD8 T cells" = "#6DB9EF",
                   "CD4 T cells" = "#EEF5FF",
                   "B cells" = "#9ADE7B")

# Reorder gg_df based on plot_order
gg_df <- gg_df %>%
  mutate(CellType = factor(CellType, levels = plot_order))

#plot
p <- ggplot(gg_df, aes(row, col, color = CellType)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_colour_manual(values = custom_colors) +
  guides(alpha = FALSE) +  
  #scale_colour_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


print(p)

#save 
ggsave(plot = p, file=file.path(outdir,
                                "cell_type.pdf"),
       width = 6, height = 4)

##
# for (i in names(geneset_OF_interest)) {
#   
#   print(i)
#   p <- ggplot(gg_df, aes(x, y), order = TRUE) +
#     geom_point(aes(color = scale(gg_df[[i]]), order = gg_df[[i]])) +
#     scale_colour_gradient2(low = "lightblue", mid = "lightgrey",  high = "red", na.value = NA) +
#     theme_bw() + ggtitle(i) +
#     labs(color = "UCell score")
#   
#   ggsave(plot = p, file=file.path(outdir,
#                                   paste0(i, ".pdf")),
#          width = 8, height = 5)
#   
# }


#gg_df <- gg_df[, -24]

# try the alternative
library(ggplot2)

# Loop through each gene set
for (i in names(geneset_OF_interest)) {
  print(i)
  
  # Sort gg_df based on the expression of the current feature
  gg_df_sorted <- gg_df[order(gg_df[[i]]), ]
  
  # Create the ggplot object
  p <- ggplot(gg_df_sorted, aes(row, col)) +
    geom_point(aes(color = scale(gg_df_sorted[[i]])), 
               na.rm = TRUE) + # Set na.rm = TRUE to remove NA values
    scale_colour_gradient2(low = "lightblue", mid = "lightgrey",  high = "red", na.value = NA) +
    theme_bw() + 
    ggtitle(i) +
    labs(color = "UCell score")
  
  # Save the plot to a PDF file
  ggsave(plot = p, 
         file = file.path(outdir, paste0(i, ".pdf")),
         width = 8, 
         height = 5)
  
}

# -----
# plot the facet ggplot
# -----
new_gg_df <- gg_df
new_gg_df <- melt(new_gg_df, measure.vars = names(geneset_OF_interest))

p <- ggplot(new_gg_df, aes(x, y)) +
  geom_point(aes(color = value)) +
  scale_colour_gradient2(low = "#EEE2DE", high = "#B31312", na.value = NA) +
  theme_bw() + facet_wrap(~variable)

p

# #############
# library(ggplot2)
# library(reshape2)  # for melt function
# 
# # Melt the data
# new_gg_df <- melt(gg_df, measure.vars = names(geneset_OF_interest))
# 
# # Arrange the data frame to have higher values in the front
# new_gg_df <- new_gg_df[order(new_gg_df$value, decreasing = TRUE), ]
# 
# # Create a custom facet grid layout
# grid_layout <- expand.grid(variable = unique(new_gg_df$variable))
# 
# # Define the number of rows and columns for the grid
# n_col <- 3  # number of columns
# n_row <- ceiling(nrow(grid_layout) / n_col)  # number of rows
# 
# # Create plots and store them in a list
# plot_list <- lapply(seq_len(nrow(grid_layout)), function(i) {
#   current_variable <- grid_layout$variable[i]
#   
#   p <- ggplot(subset(new_gg_df, variable == current_variable), aes(x, y)) +
#     geom_point(aes(color = value), size = 0.5) +
#     scale_colour_gradient2(low = "#EEE2DE", high = "#B31312", na.value = NA) +
#     theme_bw() +
#     ggtitle(current_variable) +
#     scale_x_continuous(name = "X Axis", expand = c(0, 0)) + 
#     scale_y_continuous(name = "Y Axis", expand = c(0, 0))
#   
#   return(p)
# })
# 
# # Generate the facet grid plot with different scales
# final_plot <- wrap_plots(plotlist = plot_list, ncol = n_col, scales = "free")
# 
# # Save the final plot to a PDF file
# ggsave(plot = final_plot, 
#        file = file.path(outdir, "facet_grid_plots.pdf"),
#        width = 8, 
#        height = 10)

#######

############ ASSIGN CELL TYPES BASED ON UCELL #################################

# Select the columns of interest
columns_of_interest <- c("IFN_TAMs", "Inflam_TAM", "LA_TAM",
                         "Undifferentiated-Neural.crest-like", "Neural.crest-like", "Neural.crest-like-Transitory",
                         "Transitory", "Transitory-Melanocytic", "Melanocytic",
                         "Exhausion_cell_cyle", "Exhausion_HSP", "Exhausion",
                         "Memory_effector_1", "Early_activated_cells", "Memory_effector_2",
                         "AXL_Melnoma", "MITF_Melanoma")

# Find the column name with the highest value for each row
gg_df$CellTypeFinal <- apply(gg_df[columns_of_interest], 1, function(x) {
  columns_of_interest[which.max(x)]
})

# View the updated gg_df with the new column "CellTypeFinal"
head(gg_df)

# NOW PLOT 

# Define the order of CellType for plotting
plot_order <- c("LA_TAM", "MITF_Melanoma", "Transitory-Melanocytic", "IFN_TAMs", "Inflam_TAM",
                "Undifferentiated-Neural.crest-like", "Neural.crest-like", "Neural.crest-like-Transitory",
                "Transitory", "Melanocytic",
                "Exhausion_cell_cyle", "Exhausion_HSP", "Exhausion",
                "Memory_effector_1", "Early_activated_cells", "Memory_effector_2",
                "AXL_Melnoma")

# Define custom colors for each CellType
custom_colors <- c("MITF_Melanoma" = "#F6D776",
                   "Transitory-Melanocytic" = "#FE7A36",
                   "Undifferentiated-Neural.crest-like" = "#FFF8C9",
                   "Endothelial cells" = "#508D69",
                   "Fibroblasts" = "#FFC5C5",
                   "LA_TAM" = "#FF004D",
                   "CD8 T cells" = "#6DB9EF",
                   "CD4 T cells" = "#EEF5FF",
                   "B cells" = "#9ADE7B")

# Reorder gg_df based on plot_order
gg_df <- gg_df %>%
  mutate(CellTypeFinal = factor(CellTypeFinal, levels = plot_order))

#plot
p <- ggplot(gg_df, aes(x, y, color = CellTypeFinal)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_colour_manual(values = custom_colors) +
  guides(alpha = FALSE) +  
  #scale_colour_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


print(p)

#save 
ggsave(plot = p, file=file.path(outdir,
                                "cell_type.pdf"),
       width = 8, height = 5)


# violin plot 

# plot results
plotViolinplot(
  object = spata_obj, 
  across = "bayes_space",
  variables = TGF, 
  clrp = "jama"
)


