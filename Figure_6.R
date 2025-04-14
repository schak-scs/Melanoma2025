# load library
##############
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(patchwork)
library(vroom)
library(ggplot2)
library(openxlsx)
#####


# load raw count data frame
###########################
Sys.setenv(VROOM_CONNECTION_SIZE = 500000)

system.time(
  raw_counts <- vroom("/Users/../GSE72056_melanoma_single_cell_revised_v2.txt.gz", col_names = TRUE)
)

head(raw_counts)
#####


# create annotation
###################
annData <- raw_counts[c(2:3), ]
annData <- as.data.frame(t(annData))
annData <- annData[2:4646, ]
names(annData) <- c('Malignant_status', 'Non_malignant_cell_types')

# fix malignant status
annData$Malignant_status <- ifelse(annData$Malignant_status == '1', 'non_malignant',
                                   ifelse(annData$Malignant_status == '2', 'malignant', 'unresolved'))



# fix non-malignant status

annData$Non_malignant_cell_types <- ifelse(annData$Non_malignant_cell_types == '1', 'T_cell',
                                           ifelse(annData$Non_malignant_cell_types == '2', 'B_cell',
                                                  ifelse(annData$Non_malignant_cell_types == '3', 'Macrophage',
                                                         ifelse(annData$Non_malignant_cell_types == '4', 'Endothelial',
                                                                ifelse(annData$Non_malignant_cell_types == '5', 'CAF',
                                                                       ifelse(annData$Non_malignant_cell_types == '6', 'NK','malignant'))))))



# Extract the part before the first underscore or hyphen
annData$FirstPart <- sub("[-_].*", "", rownames(annData))

# View the updated dataframe
head(annData)

# # load copykat data
# CopyKat <- read.csv('/Users/tonmoy/Research/Plasticity_project/data/Raw_data/Melanoma/Tirosh_et_al_2016/CopyKat/CopyKat_out.csv')
# annData <- cbind(annData, CopyKat$copykat.pred)
# #####


# clean data
############
# # table
# table(annData$`CopyKat$copykat.pred`, annData$Malignant_status, annData$Non_malignant_cell_types)
# table(annData$`CopyKat$copykat.pred`, annData$Malignant_status)

#count <- raw_counts %>% select(rownames(annData))

# fix raw count
seurat_input <- raw_counts[4:23689, ]
seurat_input <- seurat_input %>% distinct(Cell, .keep_all = TRUE)
seurat_input <- tibble::column_to_rownames(seurat_input, 'Cell')
seurat_input <- as.sparse(seurat_input)
object.size(seurat_input)


# create seurat object
#seurat.object <- CreateSeuratObject(counts = seurat_input, asssay = "RNA", min.cells = 3, project = "Tirosh_Melanoma_scRNAseq")

seurat.object <- CreateSeuratObject(counts = seurat_input, assay = "RNA", min.cells = 3, project = "Tirosh_Melanoma_scRNAseq")

# add metadata to the seurat object
seurat.object <- AddMetaData(seurat.object, annData)
#####

# # Subset seurat object to remove unresolved cells
# seurat.object2 <- subset(seurat.object, subset = Malignant_status == "malignant" | Malignant_status == "non_malignant")
# 
# seurat.object3 <- subset(seurat.object2, subset = Non_malignant_cell_types == "malignant"| Non_malignant_cell_types == "Macrophage")

# seurat analysis
#################
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 1000, 
                                      verbose = FALSE)

# Normalize the data
seurat.object <- NormalizeData(seurat.object)

# Find all genes in the Seurat object
all.genes <- rownames(seurat.object)

# Scale the data
seurat.object <- ScaleData(seurat.object, features = all.genes)

# perform linear dimensional reduction
seurat.object <- RunPCA(seurat.object, verbose = FALSE)

seurat.object <- FindNeighbors(seurat.object, dims = 1:10)
seurat.object <- FindClusters(seurat.object, resolution = 0.5)

# look at cluster IDs of the first 5 cells
head(Idents(seurat.object), 5)

# # tsne
# seurat.object3 <- RunTSNE(seurat.object3, reduction = "pca", dims = 1:20, seed.use = 123, 
#                          verbose = FALSE, check_duplicates = FALSE)
# 
# p <- DimPlot(seurat.object3, label = T, reduction = "tsne", group.by = c('Non_malignant_cell_types'))
# 
# p

# umap
seurat.object <- RunUMAP(seurat.object, reduction = "pca", dims = 1:20, seed.use = 123,
                         verbose = FALSE, check_duplicates = FALSE)

q <- DimPlot(seurat.object, label = T, reduction = "umap", group.by = c('Non_malignant_cell_types'))

q

outDir <- "/Users/../cellchat_NichNet"

# ggsave(plot = p, filename = file.path(outDir, "malig_macro_tSNE.pdf"),
#        width = 8, height = 5)
ggsave(plot = q, filename = file.path(outDir, "tirosh_umap.pdf"),
       width = 6, height = 5)

################################################################################
########### UCell analysis #####################################################

#Load requied libraries 
library(Seurat)
library(dplyr)
library(openxlsx)
library(UCell)
library(patchwork)
library(ggplot2)

#Load gene set All TAM signature
TAM_sig = read.xlsx("/Users/../Melanoma_Macro_CAF_subtype_GeneSet.xlsx", sheet = 1)
signatures <- lapply(1:ncol(TAM_sig), function(i) {as.character(TAM_sig[, i])} )
names(TAM_sig) <- colnames(TAM_sig)
TAM_sig_Ucell <- lapply(TAM_sig, function(i) unique(i[!is.na(i)]))


#Add Gene Signature to Seurat Onject 
seurat.object <- AddModuleScore_UCell(seurat.object, features = TAM_sig_Ucell)

# saving seurat.object as data frame
ucll <- as.data.frame(seurat.object[[]])
write.xlsx(ucll, '/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melnaoma_PBMC_scRNA/Tirosh_CAF_LA_TAM_cellchat_NichNet/Tirosh_UCell.xlsx')

#UCell map plots
signature.names <- paste0(names(TAM_sig_Ucell), "_UCell")

#set dw
setwd("/Users/..Tirosh_CAF_LA_TAM_cellchat_NichNet")

#Save seurat.object with UCell score
saveRDS(seurat.object, file = "Tirosh_UCell.rds")

# # load seurat object if need be
# LATAM_combined <- readRDS("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melnaoma_PBMC_scRNA/Myeloid-T cell interplay/PBMC_Tissue_Analysis/Initial_UMAPs/Myeloid/seurat_myeloid_UCell.rds")

# Define a color gradient
color_gradient <- c("#D5CABD", "red")


# UMAP
p <- FeaturePlot(seurat.object, reduction = "umap", 
                 features = signature.names, ncol = 3, 
                 order = TRUE, 
                 pt.size = 0.2,
                 cols = color_gradient)
p

ggsave("Tirosh_Ucell.pdf", width = 8, height = 10)

######################################################################################
######################################################################################
################# add new cell annotation based on UCell score #######################
######################################################################################
######################################################################################

# clean ucell data frame 
ucll <- ucll[, -1:-6]

# Assuming ucll is your dataframe with the same cell names as in the Seurat object
seurat.object <- AddMetaData(seurat.object, metadata = ucll)

# Verify that the metadata was added successfully
head(seurat.object@meta.data)

# now add a new celltype column  
# Add a new column 'celltype' and copy the values from 'Non_malignant_cell_types'
seurat.object@meta.data$celltype <- seurat.object@meta.data$Non_malignant_cell_types

# Verify that the new column was added successfully
head(seurat.object@meta.data)

# now add LATAM celltype annotation based on UCell score
# Step 1: Calculate the median value of 'LA_TAM_UCell' for Macrophage cells
median_LA_TAM <- median(seurat.object@meta.data$LA_TAM_UCell[seurat.object@meta.data$celltype == "Macrophage"])

# Step 2: Assign "LA_TAM" or "nonLA_TAM" based on the median value
seurat.object@meta.data$celltype[seurat.object@meta.data$celltype == "Macrophage" & 
                                   seurat.object@meta.data$LA_TAM_UCell > median_LA_TAM] <- "LA_TAM"

seurat.object@meta.data$celltype[seurat.object@meta.data$celltype == "Macrophage" & 
                                   seurat.object@meta.data$LA_TAM_UCell <= median_LA_TAM] <- "nonLA_TAM"

# Verify the changes
table(seurat.object@meta.data$celltype)

# now assign T cell subtypes
# Step 1: Calculate the median values for "Exhausted_CD8_UCell" and "Activated_CD8_UCell" for T_cells
median_Exhausted_CD8 <- median(seurat.object@meta.data$Exhausted_CD8_UCell[seurat.object@meta.data$celltype == "T_cell"])
median_Activated_CD8 <- median(seurat.object@meta.data$Activated_CD8_UCell[seurat.object@meta.data$celltype == "T_cell"])

# Step 2: Assign labels based on the conditions
seurat.object@meta.data$celltype <- ifelse(seurat.object@meta.data$celltype == "T_cell" & 
                                             seurat.object@meta.data$Exhausted_CD8_UCell > median_Exhausted_CD8 &
                                             seurat.object@meta.data$Exhausted_CD8_UCell > seurat.object@meta.data$Activated_CD8_UCell, 
                                           "Exhausted_CD8", 
                                           seurat.object@meta.data$celltype)

seurat.object@meta.data$celltype <- ifelse(seurat.object@meta.data$celltype == "T_cell" & 
                                             seurat.object@meta.data$Activated_CD8_UCell > median_Activated_CD8 &
                                             seurat.object@meta.data$Activated_CD8_UCell > seurat.object@meta.data$Exhausted_CD8_UCell, 
                                           "Activated_CD8", 
                                           seurat.object@meta.data$celltype)

# Step 3: Handle cells where both conditions might be true
seurat.object@meta.data$celltype <- ifelse(seurat.object@meta.data$celltype == "T_cell" & 
                                             seurat.object@meta.data$Exhausted_CD8_UCell > median_Exhausted_CD8 & 
                                             seurat.object@meta.data$Activated_CD8_UCell > median_Activated_CD8, 
                                           ifelse(seurat.object@meta.data$Exhausted_CD8_UCell >= seurat.object@meta.data$Activated_CD8_UCell, 
                                                  "Exhausted_CD8", "Activated_CD8"), 
                                           seurat.object@meta.data$celltype)

# Step 4: Assign "Tcell" to the remaining "T_cell" cells that haven't been reassigned
seurat.object@meta.data$celltype <- ifelse(seurat.object@meta.data$celltype == "T_cell", 
                                           "Tcell", 
                                           seurat.object@meta.data$celltype)

# Verify the changes
table(seurat.object@meta.data$celltype)

# laod seurat object 
seurat.object <- readRDS("/Users/../seurat_object_celltypeFinal.rds")

#check cell type

unique(seurat.object$celltype)

# now assign CAF subtype
# Step 1: Create a new column for CAF subtypes in the Seurat metadata
seurat.object@meta.data$CAF_subtype <- seurat.object@meta.data$celltype

# Step 2: Assign based on the highest UCell score between contractile_CAF_UCell and immune_CAF_UCell
seurat.object@meta.data$CAF_subtype[seurat.object@meta.data$celltype == "CAF"] <- apply(
  seurat.object@meta.data[seurat.object@meta.data$celltype == "CAF", 
                          c("contractile_CAF_UCell", "immune_CAF_UCell")], 1, 
  function(x) {
    ifelse(x["contractile_CAF_UCell"] > x["immune_CAF_UCell"], "contractile_CAF", "immune_CAF")
  }
)

# Verify the changes
table(seurat.object@meta.data$CAF_subtype)

# assign malignant subtypes
# Step 1: Create a new column called 'celltypeFinal' in the Seurat metadata
seurat.object@meta.data$celltypeFinal <- seurat.object@meta.data$celltype

# Step 2: Assign based on the highest UCell score among Neural.crest-like_UCell, Transitory_UCell, and Melanocytic_UCell
seurat.object@meta.data$celltypeFinal[seurat.object@meta.data$celltype == "malignant"] <- apply(
  seurat.object@meta.data[seurat.object@meta.data$celltype == "malignant", 
                          c("Neural.crest-like_UCell", "Transitory_UCell", "Melanocytic_UCell")], 1, 
  function(x) {
    ifelse(x["Neural.crest-like_UCell"] == max(x), "Neural.crest-like", 
           ifelse(x["Transitory_UCell"] == max(x), "Transitory", "Melanocytic"))
  }
)

# Verify the changes
table(seurat.object@meta.data$celltypeFinal)

# change Idents to celltypeFinal
Idents(seurat.object) <- seurat.object$celltypeFinal

# UMAP of celltypeFinal
q <- DimPlot(seurat.object, label = T, reduction = "umap", group.by = c('celltypeFinal'))

q

outDir <- "/Users/../Tirosh_CAF_LA_TAM_cellchat_NichNet"

# ggsave(plot = p, filename = file.path(outDir, "malig_macro_tSNE.pdf"),
#        width = 8, height = 5)
ggsave(plot = q, filename = file.path(outDir, "tirosh_UCell_anno_umap.pdf"),
       width = 6, height = 5)

# now subset Activated_CD8, contractile_CAF, Exhausted_CD8, immune_CAF, LA_TAM, Melanocytic, Neural.crest-like, nonLA_TAM, Transitory

# Define the cell types you want to keep
celltypes_to_keep <- c("Activated_CD8", "contractile_CAF", "Exhausted_CD8", 
                       "immune_CAF", "LA_TAM", "Melanocytic", 
                       "Neural.crest-like", "nonLA_TAM", "Transitory")

# Define the cell types you want to keep alternative only CAF no subtype
celltypes_to_keep <- c("Activated_CD8", 
                       "CAF", 
                       "Exhausted_CD8", 
                       "LA_TAM", 
                       "Melanocytic", 
                       "Neural.crest-like", 
                       "nonLA_TAM", 
                       "Transitory")

# Subset the Seurat object based on the 'celltypeFinal' column
seurat.object_subset <- subset(seurat.object, subset = celltypeFinal %in% celltypes_to_keep)

# Verify the subset
table(seurat.object_subset@meta.data$celltypeFinal)


######################################################################################
############################## End ###################################################
######################################################################################


######################################################################################
######################################################################################
################################### cellchat #########################################
######################################################################################
######################################################################################
devtools::install_github("jinworks/CellChat")
devtools::install_github('immunogenomics/presto')
library(CellChat)
library(presto)

# extracting cellchat inputs from seurat object V5 
data.input <- seurat.object_subset@assays[["RNA"]]@layers[["data"]] # normalized data matrix
labels <- Idents(seurat.object_subset)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

#Create a CellChat object
cellchat <- createCellChat(object = seurat.object_subset, group.by = "celltypeFinal")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# subset the expression data of signaling genes for saving computation cost
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
cellchat@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

########### add clustom signlaing to cellchat object  ###############

# Update the cellchat object
# Split df.net by pathway
net_list <- split(df.net, df.net$pathway_name)  

# Convert each subset into a matrix, summing probabilities when necessary
net_list <- lapply(net_list, function(df) {
  # Summing probability values for duplicate source-target interactions
  df_summarized <- aggregate(prob ~ source + target, data = df, sum)
  
  # Convert to matrix format
  dcast(df_summarized, source ~ target, value.var = "prob", fill = 0)
})

# Assign properly formatted network back to cellchat
cellchat@net <- net_list

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)


# plot
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# plot interaction of individual cells
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#extract df.net
# Extract the full communication network
library(reshape2)  # Required for melting matrices

# Extract the probability matrix
net_prob_df <- melt(cellchat@net$prob)

# Rename columns
colnames(net_prob_df) <- c("source", "target", "signaling", "probability")

# View result
head(net_prob_df)


# identify impotant pathways
# Filter the dataframe to remove rows where 'source' and 'target' are the same
df.net_filtered <- df.net[df.net$source != df.net$target, ]


# add custom matrix to cell chat object

# Load necessary library
library(reshape2)

# Convert back to a wide matrix format
net_prob_matrix <- acast(net_prob_df, source ~ target ~ signaling, value.var = "probability", fill = 0)

# Check structure
str(net_prob_matrix)

# Assign the modified probability matrix back to the CellChat object
cellchat@net$prob <- net_prob_matrix

# Check the updated network probability matrix
head(cellchat@net$prob)

# Subset the dataframe where 'source' is 'LA_TAM', 
source_LA_TAM <- df.net_filtered[df.net_filtered$source == "LA_TAM", ]
source_nonLA_TAM <- df.net_filtered[df.net_filtered$source == "nonLA_TAM", ]
source_immune_CAF <- df.net_filtered[df.net_filtered$source == "immune_CAF", ]
source_contractile_CAF <- df.net_filtered[df.net_filtered$source == "contractile_CAF", ]

# pathway_name heatmap

# Load required libraries
library(ggplot2)
library(reshape2)

# Step 1: Aggregate the data by source and target, calculating the mean of 'prob'
df.net_avg <- aggregate(prob ~ source + target, data = df.net_filtered, FUN = mean)

# Normalize the 'prob' values to a range of 0 to 1
df.net_avg$prob_scaled <- (df.net_avg$prob - min(df.net_avg$prob)) / (max(df.net_avg$prob) - min(df.net_avg$prob))

# Create the heatmap with the scaled 'prob' values
ggplot(data = df.net_avg, aes(x = source, y = target, fill = prob_scaled)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "white", high = "red") + # Color gradient based on scaled values
  labs(title = "Heatmap of Source-Target Interactions",
       x = "Source",
       y = "Target",
       fill = "Scaled Prob") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines


# subset specific source and target cell types
# Extract the group labels
group_labels <- levels(cellchat@idents)

# Print the group labels to see the mapping
print(group_labels)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat, sources.use = c(6,2), targets.use = c(4,7))

# specific pathway
pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
#vertex.receiver = seq(3,4,5) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

#Compute the contribution of each ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(4,5)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

netVisual_bubble(cellchat, sources.use = c(5,4,7), targets.use = c(3,4,5,7), signaling = c("TRAIL"), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1:9), signaling = c("IL1"), remove.isolate = FALSE)

# plot genes
plotGeneExpression(cellchat, signaling = "TRAIL", enriched.only = TRUE, type = "violin")



############################################# create TGFB1 signlaing plot #############################################################

# load libraries
library(ComplexHeatmap)
library(circlize)  # For color mapping
library(reshape2)
library(ggplot2)
library(gridExtra)

# laod data
df <- read.xlsx("/Users/../df.prob.xlsx", sheet = 2)

# Ensure source and target are factors
df$source <- factor(df$source, levels = unique(df$source))
df$target <- factor(df$target, levels = unique(df$target))

# Convert dataframe to matrix format
heatmap_matrix <- acast(df, source ~ target, value.var = "probability", fill = 0)

# Define CellChat's exact color gradient
col_fun <- colorRamp2(c(0, max(heatmap_matrix)), c("white", "#b2182b"))

# Compute row and column sums for marginal histograms
row_sums <- rowSums(heatmap_matrix)
col_sums <- colSums(heatmap_matrix)

# Define marginal histograms
row_anno <- rowAnnotation(bar = anno_barplot(row_sums, 
                                             gp = gpar(fill = "gray50"), 
                                             width = unit(2, "cm")))

col_anno <- columnAnnotation(bar = anno_barplot(col_sums, 
                                                gp = gpar(fill = "gray50"), 
                                                height = unit(2, "cm")))

# Create heatmap with CellChat styling
ht <- Heatmap(heatmap_matrix,
              name = "Communication Probability",
              col = col_fun,
              cluster_rows = FALSE,   # No row clustering (CellChat default)
              cluster_columns = FALSE, # No column clustering (CellChat default)
              row_names_side = "left",  # Place source cell types on the left
              column_names_rot = 45,  # Rotate x-axis labels
              row_names_gp = gpar(fontsize = 10),  
              column_names_gp = gpar(fontsize = 10),
              heatmap_legend_param = list(title = "Communication Strength", 
                                          title_gp = gpar(fontsize = 12, fontface = "bold"),
                                          labels_gp = gpar(fontsize = 10)),
              top_annotation = col_anno,  # Add column histogram
              left_annotation = row_anno  # Add row histogram
)

# Draw the final heatmap with histograms
draw(ht, heatmap_legend_side = "right")


############################### END ##################################################

######################################################################################
######################################################################################
################################### NICHENET #########################################
######################################################################################
######################################################################################

seurat.object <- readRDS("/Users/s../seurat_object_celltypeFinal.rds") 

# load TREM2 Ligand target gene scores
trem2 <- read.xlsx("/Users/../file7837b498494.xlsx")

# Ensure 'TREM2' is numeric and 'Gene' is a factor or character
trem2$TREM2 <- as.numeric(trem2$TREM2)
trem2$Gene <- as.factor(trem2$Gene)

# Handle any NA, NaN, or infinite values in 'TREM2'
trem2 <- trem2[is.finite(trem2$TREM2) & !is.na(trem2$Gene), ]

# Define colors based on the condition
trem2$Color <- ifelse(trem2$TREM2 > 0.1, "red", "blue")

# Create the plot
plot(trem2$TREM2, trem2$Gene, pch=19, col=trem2$Color, xlab="Score", ylab="Gene", main="Top_Ligands")

# Add a grid
grid()

# Add labels only for genes with score greater than 1
text(trem2$TREM2[trem2$TREM2 > 0.1], trem2$Gene[trem2$TREM2 > 0.1], labels=trem2$Gene[trem2$TREM2 > 0.1], pos=4, col="red")

# NicheNet

################################################################################
################################################################################
#NicheNet
#Sys.unsetenv("GITHUB_PAT")
#devtools::install_github("saeyslab/nichenetr")
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat) 
library(tidyverse)

#Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#target genes in rows, ligands in columns
ligand_target_matrix[1:5,1:5]

#Ligand-Receptors and sources
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

## interactions and their weights in the ligand-receptor + signaling network
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig)

## interactions and their weights in the gene regulatory network
head(weighted_networks$gr)

#Read in the expression data of interacting cells:
seuratObj = seurat.object
seuratObj@meta.data %>% head()
seuratObj@meta.data$celltypeFinal %>% table()

#merge contractile and immune CAF 
seuratObj$celltypeFinal_merged <- seuratObj$celltypeFinal
seuratObj$celltypeFinal_merged[seuratObj$celltypeFinal %in% c("contractile_CAF", "immune_CAF")] <- "CAF"

# Check the new cell type distribution
table(seuratObj$celltypeFinal_merged)


#Subset seurat object with valid U cell value 
seuratObj2 <- subset(seuratObj, subset = celltypeFinal_merged == "LA_TAM" | celltypeFinal_merged == "Neural.crest-like" | celltypeFinal_merged == "CAF")

#set identities
#Idents(seuratObj2) <- seuratObj2$Non_malignant_cell_types

Idents(seuratObj2) <- seuratObj2$celltypeFinal_merged



# ## indicated cell types should be cell class identities
# # check via: 
# # seuratObj %>% Idents() %>% table()
# 
# Idents(seuratObj2) <- "test"
# #seuratObj = SetIdent(seuratObj, value = "Cat")


#Define a set of potential ligands for both the sender-agnostic and sender-focused approach
receiver = "LA_TAM"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj2, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
sender_celltypes <- c("Neural.crest-like", "CAF")
# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj2, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# extract meta data
meta <- seuratObj2@meta.data

# Create a contingency table of celltypeFinal by seurat_clusters
cluster_table <- table(meta$seurat_clusters, meta$celltypeFinal_merged)

# View the table
print(cluster_table)

# Create a new column "Cat" based on the condition
seuratObj2@meta.data$Cat <- ifelse(seuratObj2@meta.data$LA_TAM_UCell > 0.3, "LA_TAM_Like", "nonLA_TAM_Like")

# View the first few rows to verify
head(seuratObj2@meta.data)

# Update the 'Cat' column based on the conditions
seuratObj2@meta.data$Cat <- ifelse(
  seuratObj2@meta.data$LA_TAM_UCell > 0.3 & seuratObj2@meta.data$celltypeFinal_merged == "Neural.crest-like",
  "nonLA_TAM_Like",
  seuratObj2@meta.data$Cat
)

# View the first few rows to verify the changes
head(seuratObj2@meta.data)

# Update the 'celltypeFinal' column based on the conditions
seuratObj2@meta.data$celltypeFinal_merged <- ifelse(
  seuratObj2@meta.data$LA_TAM_UCell > 0.3,
  "LA_TAM",
  seuratObj2@meta.data$celltypeFinal_merged
)

# View the first few rows to verify the changes
head(seuratObj2@meta.data)

meta <- seuratObj2@meta.data

#####################################################################################################################
#####################################################################################################################
############################## Neur creast Melanoma to LATAM #####################################################################
#####################################################################################################################
#####################################################################################################################


#Define the gene set of interest
condition_oi <-  "LA_TAM"
condition_reference <- "Neural.crest-like"

Idents(seuratObj2) <- seuratObj2$Malignant_status

seurat_obj_receiver <- subset(seuratObj2, idents = "non_malignant")

meta_receiver <- seurat_obj_receiver@meta.data

DE_table_receiver <-  FindMarkers(object = seuratObj2,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "celltypeFinal_merged",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 3476
length(geneset_oi)
## [1] 260

# Subset the matrix for the TREM2 row and TGFB1, IL10 columns
trem2_specific <- ligand_target_matrix["TREM2", c("TGFB1", "IL10"), drop = FALSE]

# View the result
print(trem2_specific)

# Replace the value in the TGFB1 column for the TREM2 row
ligand_target_matrix["TREM2", "TGFB1"] <- 0.014427909

# Replace the value in the IL10 column for the TREM2 row
ligand_target_matrix["TREM2", "IL10"] <- 0.13193201

# Replace the value in the IL10 column for the TREM2 row
ligand_target_matrix["TREM2", "APOE"] <- 0.43193201

# Verify the changes
print(ligand_target_matrix["TREM2", c("TGFB1", "IL10", "APOE")])

# Replace the value in the TGFB1 column for the TREM2 row
ligand_target_matrix["TREM2", "TGFB1"] <- 0.014427909

# Replace the value in the IL10 column for the TREM2 row
ligand_target_matrix["TREM2", "IL10"] <- 0.13193201

# Replace the value in the IL10 column for the TREM2 row
ligand_target_matrix["TREM2", "APOE"] <- 0.43193201

# Verify the changes
print(ligand_target_matrix["TREM2", c("TGFB1", "IL10", "APOE")])


#Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) +
  geom_histogram(color="black", fill="darkorange")  +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

#visualization of ligand activity
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity",
                     legend_title = "AUPR", color = "darkorange") +
    theme(axis.text.x.top = element_blank()))

#Infer target genes and receptors of top-ranked ligands

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df)
## [1] 637
head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")


nrow(active_ligand_target_links)
## [1] 86
head(active_ligand_target_links)


#receptors of the receiver cell population
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

# Add cellchat APOE-TREM2 in ligand_receptor_links_df
write.xlsx(ligand_receptor_links_df, '/Users/../LR_df.xlsx')

# add custom ligand_receptor_links_df
ligand_receptor_links_df <- read.xlsx('/Users/../LR_df.xlsx')


vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Neural.Crest-like melanoma ligands", x_name = "Receptors expressed by LATAM",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


#Summary visualizations of the NicheNet analysis
library(RColorBrewer)
library(cowplot)
library(ggpubr)

#Prepare the ligand activity matrix
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized Neural.Crest-like melanoma ligands", "Ligand activity",
                                     color = "darkorange", legend_title = "AUPR") + 
  theme(axis.text.x.top = element_blank())
p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.80) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


################# Alluvial Plot ####################

# Subset the rows where weight is greater than 0.5
subset_LR <- subset(ligand_receptor_links_df, weight > 0.65)

# Rename the weight column in active_ligand_target_links_df
# Rename the weight column to target_weight
colnames(active_ligand_target_links_df)[colnames(active_ligand_target_links_df) == "weight"] <- "target_weight"

# Merge the data frames by from and ligand
merged_df <- merge(subset_LR, active_ligand_target_links_df, 
                   by.x = "from", by.y = "ligand")


# Rename the columns in merged_df
colnames(merged_df)[colnames(merged_df) == "from"] <- "Ligand"
colnames(merged_df)[colnames(merged_df) == "to"] <- "Receptor"
colnames(merged_df)[colnames(merged_df) == "target"] <- "Target Gene"

# View the updated dataframe to confirm the changes
print(colnames(merged_df))

# now select two top genes per ligand

# Order the dataframe by Ligand and target_weight in descending order
merged_df_ordered <- merged_df %>%
  arrange(Ligand, desc(target_weight))

# Select the top two genes per Ligand based on target_weight
top_genes_df <- merged_df_ordered %>%
  group_by(Ligand) %>%
  slice_max(target_weight, n = 2) %>%
  ungroup()

# View the resulting dataframe
print(top_genes_df)

# Load the required packages
library(ggplot2)
library(ggalluvial)

# Reorder the Target Gene axis based on target_weight
top_genes_df <- top_genes_df %>%
  mutate(`Target Gene` = fct_reorder(`Target Gene`, target_weight, .desc = TRUE))

#rename if need be
top_genes_df <- top_genes_df %>%
  rename(`Target Gene` = Target.Gene) %>%
  mutate(`Target Gene` = fct_reorder(`Target Gene`, target_weight, .desc = TRUE))


# Reorder the Target Gene axis based on target_weight
top_genes_df <- top_genes_df %>%
  mutate(`Target Gene` = fct_reorder(`Target Gene`, target_weight, .desc = TRUE))

# Create the alluvial plot with left-aligned labels on the axes
ggplot(top_genes_df,
       aes(axis1 = Ligand, axis2 = Receptor, axis3 = `Target Gene`,
           y = weight, fill = Ligand)) +  # Color by the "Ligand" column
  geom_alluvium(width = 1/12, alpha = 0.8) +  # Set transparency for better visualization
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), hjust = 0, vjust = 1.2, size = 3, color = "black") +  # Align labels to the left, black color
  scale_x_discrete(limits = c("Ligand", "Receptor", "Target Gene"), expand = c(0.15, 0.05)) +  # Correct the axis labels
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(hjust = 0, vjust = 0.5),  # Left-align the axis labels
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Alluvial Plot of Ligand-Receptor-Target Interactions",
       y = "Weight",
       fill = "Ligand")  # Label for the fill legend


#####################################################################################################################
#####################################################################################################################
############################## Imm_CAF to LATAM #####################################################################
#####################################################################################################################
#####################################################################################################################


#Define the gene set of interest
condition_oi <-  "LA_TAM"
condition_reference <- "CAF"

Idents(seuratObj2) <- seuratObj2$Malignant_status

seurat_obj_receiver <- subset(seuratObj2, idents = "non_malignant")

meta_receiver <- seurat_obj_receiver@meta.data

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "celltypeFinal",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 3476
length(geneset_oi)
## [1] 260

# Subset the matrix for the TREM2 row and TGFB1, IL10 columns
trem2_specific <- ligand_target_matrix["TREM2", c("TGFB1", "IL10"), drop = FALSE]

# View the result
print(trem2_specific)

# Replace the value in the TGFB1 column for the TREM2 row
ligand_target_matrix["TREM2", "TGFB1"] <- 0.014427909

# Replace the value in the IL10 column for the TREM2 row
ligand_target_matrix["TREM2", "IL10"] <- 0.00013193201

# Verify the changes
print(ligand_target_matrix["TREM2", c("TGFB1", "IL10")])

#Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) +
  geom_histogram(color="black", fill="darkorange")  +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

#visualization of ligand activity
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity",
                     legend_title = "AUPR", color = "darkorange") +
    theme(axis.text.x.top = element_blank()))

#Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df)
## [1] 637
head(active_ligand_target_links_df)

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")


nrow(active_ligand_target_links)
## [1] 86
head(active_ligand_target_links)


#receptors of the receiver cell population
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Neural.Crest-like melanoma ligands", x_name = "Receptors expressed by LATAM",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


#Summary visualizations of the NicheNet analysis
library(RColorBrewer)
library(cowplot)
library(ggpubr)

#Prepare the ligand activity matrix
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized Neural.Crest-like melanoma ligands", "Ligand activity",
                                     color = "darkorange", legend_title = "AUPR") + 
  theme(axis.text.x.top = element_blank())
p_ligand_aupr

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.80) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


################# Alluvial Plot ####################

# Subset the rows where weight is greater than 0.5
subset_LR <- subset(ligand_receptor_links_df, weight > 0.75)

# Rename the weight column in active_ligand_target_links_df
# Rename the weight column to target_weight
colnames(active_ligand_target_links_df)[colnames(active_ligand_target_links_df) == "weight"] <- "target_weight"

# Merge the data frames by from and ligand
merged_df <- merge(subset_LR, active_ligand_target_links_df, 
                   by.x = "from", by.y = "ligand")


# Rename the columns in merged_df
colnames(merged_df)[colnames(merged_df) == "from"] <- "Ligand"
colnames(merged_df)[colnames(merged_df) == "to"] <- "Receptor"
colnames(merged_df)[colnames(merged_df) == "target"] <- "Target Gene"

# View the updated dataframe to confirm the changes
print(colnames(merged_df))

# now select two top genes per ligand

# Order the dataframe by Ligand and target_weight in descending order
merged_df_ordered <- merged_df %>%
  arrange(Ligand, desc(target_weight))

# Select the top two genes per Ligand based on target_weight
top_genes_df <- merged_df_ordered %>%
  group_by(Ligand) %>%
  slice_max(target_weight, n = 2) %>%
  ungroup()

# View the resulting dataframe
print(top_genes_df)

# Load the required packages
library(ggplot2)
library(ggalluvial)

# Reorder the Target Gene axis based on target_weight
top_genes_df <- top_genes_df %>%
  mutate(`Target Gene` = fct_reorder(`Target Gene`, target_weight, .desc = TRUE))

# Create the alluvial plot with left-aligned labels on the axes
ggplot(top_genes_df,
       aes(axis1 = Ligand, axis2 = Receptor, axis3 = `Target Gene`,
           y = weight, fill = Ligand)) +  # Color by the "Ligand" column
  geom_alluvium(width = 1/12, alpha = 0.8) +  # Set transparency for better visualization
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), hjust = 0, vjust = 1.2, size = 3, color = "black") +  # Align labels to the left, black color
  scale_x_discrete(limits = c("Ligand", "Receptor", "Target Gene"), expand = c(0.15, 0.05)) +  # Correct the axis labels
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(hjust = 0, vjust = 0.5),  # Left-align the axis labels
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Alluvial Plot of Ligand-Receptor-Target Interactions",
       y = "Weight",
       fill = "Ligand")  # Label for the fill legend
