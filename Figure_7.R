#Load requied libraries 
library(Seurat)
library(dplyr)
library(openxlsx)
library(UCell)
library(patchwork)
library(ggplot2)
library(magrittr)
library(ArchR)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(viridisLite)
library(rstatix)
library(ggsignif)
library(Nebulosa)
library(BiocFileCache)

# Load rds file Seurat object (Provided by TD and CK)
#feldman <- readRDS("/Users/sajibchakraborty/Documents/EkoEstMed/Immune Sub Type/Melanoma/scRNA-seq/feldman_seurat.rds")
# Load annotation file 
ann <- read.delim("/Users/../GSE146613_metadata.tsv")


#create 1st seurat object APOE2
if (FALSE) {
  data_dir <- '/Users/s../APOE variant/APOE2'
  list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
  data <- Read10X(data.dir = data_dir)
  seurat_object = CreateSeuratObject(counts = data)
}

APOE2 <- seurat_object

#create 2nd seurat object APOE4
if (FALSE) {
  data_dir <- '/Users/../APOE variant/APOE4'
  list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
  data <- Read10X(data.dir = data_dir)
  APOE4 = CreateSeuratObject(counts = data)
}

# create two ann data for two seurat object 
ann_APOE2 <- subset(ann, subset = genotype == "E2")
ann_APOE4 <- subset(ann, subset = genotype == "E4")

#rename ann dataframes CEll ID barcodes: replace 2 with 1 at the end, and replace _ with -
library(stringr)
ann_APOE4$cell_ID <- str_replace(ann_APOE4$cell_ID, "2", "1") # replace 2 with 1
ann_APOE4$cell_ID <- str_replace(ann_APOE4$cell_ID, "_", "-") # replace _ with - for apoe4
ann_APOE2$cell_ID <- str_replace(ann_APOE2$cell_ID, "_", "-") # replace _ with - for apoe2

# subset row of APOE seurat object based on ann_APOE cell Identity

# ADD ROWNAMES TO ANN DATA OTHERWISE IT WILL NOT WORK
rownames(ann_APOE2) <- ann_APOE2$cell_ID
rownames(ann_APOE4) <- ann_APOE4$cell_ID

# # Add cell annotation to seurat objects APOE2 and APOE4  
APOE2 <- AddMetaData(APOE2, ann_APOE2, col.name = NULL)
APOE4 <- AddMetaData(APOE4, ann_APOE4, col.name = NULL)

#merge two seurat objects APOE2+APOE4
seurat_object <- merge(APOE2, y = APOE4, add.cell.ids = c("APOE2", "APOE4"), project = "APOE_VAriant")

# again subset based on E2 and E4 genotype # THe prbolem is seurat object contains cells without metadata

#match cell ID in ann data and Seurat object

# 
# # Subset seurat object by taking the 
# macro <- subset(feldman, subset = Cluster.number == "3")

# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- NormalizeData(seurat_object)

#Identification of highly variable features (feature selection)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

#Perform linear dimensional reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))


# Examine and visualize PCA results a few different ways
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
seurat_cluster <- RunUMAP(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "pca")

#PCA heatmap
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
JackStrawPlot(seurat_object, dims = 1:15)


#Cluster the cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

#UMAP only
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
UMAP <- DimPlot(seurat_object, reduction = "umap")

#Genotype UMAP
gentype <- DimPlot(seurat_object, reduction = "umap", group.by= "genotype")

# Pre and Post
cellType <- DimPlot(seurat_object, reduction = "umap", group.by= "celltype")



################################################# END OF INITIAL ANALYSIS #############################################################

################################################# LAOD SAVED SEURAT OBJECT  #############################################################

# LAOD 
seurat_object <- readRDS("/Users/../Seurat_Object_APOE.rds")

# remove cells with celltype

# Check if 'celltype' metadata column exists in the Seurat object
if ("celltype" %in% colnames(seurat_object@meta.data)) {
  # Subset the Seurat object to exclude cells with NA or empty 'celltype'
  seurat_object <- subset(seurat_object, subset = !(is.na(celltype) | celltype == ""))
} else {
  print("The 'celltype' metadata column does not exist in the Seurat object.")
}




# test the laoded seurat object 
#seurat_object <- RunUMAP(seurat_object, dims = 1:10)

UMAP <- DimPlot(seurat_object, reduction = "umap")
UMAP

#ggsave(cellType, filename ="/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma_Spatial_Trans and APOE variant/APOE variant/Celltype.pdf", width = 8, height = 5)

#Genotype UMAP
gentype <- DimPlot(seurat_object, reduction = "umap", group.by= "genotype")
gentype

# Test celltype
cellType <- DimPlot(seurat_object, reduction = "umap", group.by= "celltype")
cellType


################################################# END OF SEURAT OBJECT LOADING AND TESTING #############################################################

################################################# START OF UCELL FOR LA-TAM ANNOTATION  #############################################################

# ################################################################################
# ################################## UCELL #######################################
# 
# #Load requied libraries 
library(Seurat)
library(dplyr)
library(openxlsx)
library(UCell)
library(patchwork)
library(ggplot2)
# 
# #Load gene set All TAM signature
TAM_sig = read.xlsx("/Users/..t/LM22_LA_TAM_Mouse.xlsx", sheet = 4)
signatures <- lapply(1:ncol(TAM_sig), function(i) {as.character(TAM_sig[, i])} )
names(TAM_sig) <- colnames(TAM_sig)
TAM_sig_Ucell <- lapply(TAM_sig, function(i) unique(i[!is.na(i)]))
# 
# #Add Gene Signature to Seurat Onject 
seurat_object <- AddModuleScore_UCell(seurat_object, features = TAM_sig_Ucell)
# 
# #Save seurat.object with UCell score
# saveRDS(macro, file = "seurat_APOE_UCell_SC.rds")
# 
# # saving seurat.object as data frame
ucll <- as.data.frame(seurat_object[[]])
# 
# #UCell map plots
signature.names <- paste0(names(TAM_sig_Ucell), "_UCell")


# set idents of seurat object 
Idents(object = seurat_object) <- "celltype"

# 
# # UMAP
p <- FeaturePlot(seurat_object, reduction = "umap",
                 features = signature.names, ncol = 3,
                 #group.by = celltype,
                 split.by = genotype,
                 order = F,
                 pt.size = 0.3,
                 cols = c("darkred", "orange", "yellow", "white"))
p

#define color patelettes
my_palette <- ArchRPalettes$beach

# UMAP Feature plot 

p1 <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')

p2 <- (FeaturePlot(seurat_object, features = c("LA_TAM_UCell", "IFN_TAMs_UCell"), split.by = "genotype", pt.size = 0.2))


p1 | p2

# violin plot
dev.off()
VlnPlot(object = seurat_object, features = c("LA_TAM_UCell", "IFN_TAMs_UCell"), group.by = "celltype", pt.size = 0, split.by = "genotype")

# dot plot
DotPlot(seurat_object, features = c("LA_TAM_UCell", "IFN_TAMs_UCell"), cols = c("blue", "lightgrey"), split.by = "genotype") + RotatedAxis()

# plot density
# but before I need to subset CKD and Normal
# Subset cells based on disease column in metadata
E2 <- subset(seurat_object, subset = genotype == "E2")
E4 <- subset(seurat_object, subset = genotype == "E4")

# Subset the Seurat object: CKD
p1 <- plot_density(E2, features = c("LA_TAM_UCell", "IFN_TAMs_UCell"), size = 0.5)

p2 <- plot_density(E4, features = c("LA_TAM_UCell", "IFN_TAMs_UCell"), size = 0.5)

p1 | p2

DotPlot(E2, features = c("Apoe", "Trem2", "Cd163", "Ccl2", "Ccl4", "Cd68"), cols = "RdYlBu") + RotatedAxis()
DotPlot(E4, features = c("Apoe", "Trem2", "Cd163", "Ccl2", "Ccl4", "Cd68"), cols = "RdYlBu") + RotatedAxis()

################################################# END OF UCELL  #############################################################

################################################# START OF CELL CHAT  #############################################################

# install
devtools::install_github("sqjin/CellChat")
install.packages('NMF')
install.packages("circlize")

# laod
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Subset cells based on genotype column  in metadata
E2 <- subset(seurat_object, subset = genotype == "E2")
E4 <- subset(seurat_object, subset = genotype == "E4")


#set Idents of E2 and E4
#Idents(object = seurat_object) = "celltype"
Idents(object = E2) <- "celltype"
Idents(object = E4) <- "celltype"

################################# CREATING CELLCHAT DATABASE ################################################

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")) # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

################################### END OF CELLCHAT DATABASE #################################################

# # Extract the CellChat input files from a Seurat V3 object
# # For seurat_object
# data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
# labels <- Idents(seurat_object)
# meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
# # column name of meta should labels
# names(meta) <- "labels"

################################### CREATING CELLCHATE2 ONJECT AND RUNNING CELLCHAT #############################################

# Extract the CellChat input files from a Seurat V3 object
# For E2
data.inputE2 <- GetAssayData(E2, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(E2)
metaE2 <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
# column name of meta should labels
names(metaE2) <- "labels"

# Create a CellChat object using data matrix as input

# For E2
cellchatE2 <- createCellChat(object = data.inputE2, meta = metaE2, group.by = "labels")

#Add cell information into meta slot of the object
#for E2
cellchatE2 <- addMeta(cellchatE2, meta = metaE2, meta.name = "labels")
cellchatE2 <- setIdent(cellchatE2, ident.use = "labels") # set "labels" as default cell identity
levels(cellchatE2@idents) # show factor levels of the cell labels
groupSizeE2 <- as.numeric(table(cellchatE2@idents)) # number of cells in each cell group

# add database to cellchat object 

# set the used database in the object
cellchatE2@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
#for E2
cellchatE2 <- subsetData(cellchatE2) # This step is necessary even if using the whole database


# identify overexpression of genes and interactions
cellchatE2 <- identifyOverExpressedGenes(cellchatE2)
cellchatE2 <- identifyOverExpressedInteractions(cellchatE2)

# # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.mouse)

#Compute the communication probability and infer cellular communication network
cellchatE2 <- computeCommunProb(cellchatE2)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchatE2 <- filterCommunication(cellchatE2, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
df.netE2 <- subsetCommunication(cellchatE2)
df.netE2 <- subsetCommunication(cellchatE2, sources.use = c(2,11), targets.use = c(16,19,21,23,24))
#df.net <- subsetCommunication(cellchat, signaling = c("IL10"))


#Infer the cell-cell communication at a signaling pathway level
cellchatE2 <- computeCommunProbPathway(cellchatE2)

#Calculate the aggregated cell-cell communication network
cellchatE2 <- aggregateNet(cellchatE2)

groupSize <- as.numeric(table(cellchatE2@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatE2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatE2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf("cellchat_E2.pdf", width = 20, height = 20)
mat <- cellchatE2@net$weight
par(mfrow = c(6,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("CCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchatE2, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchatE2, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchatE2, signaling = pathways.show, layout = "chord")


################################ END OF CELLCHATE2 OBJECT #######################################

################################ START OF CELLCHATE4 OBJECT #######################################

#For E4
data.inputE4 <- GetAssayData(E4, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(E4)
metaE4 <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
# column name of meta should labels
names(metaE4) <- "labels"


# Create a CellChat object using data matrix as input
# For E4
cellchatE4 <- createCellChat(object = data.inputE4, meta = metaE4, group.by = "labels")


#for E4
cellchatE4 <- addMeta(cellchatE4, meta = metaE4, meta.name = "labels")
cellchatE4 <- setIdent(cellchatE4, ident.use = "labels") # set "labels" as default cell identity
levels(cellchatE4@idents) # show factor levels of the cell labels
groupSizeE4 <- as.numeric(table(cellchatE4@idents)) # number of cells in each cell group


# set the used database in the object
cellchatE4@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
#for E2
cellchatE4 <- subsetData(cellchatE4) # This step is necessary even if using the whole database

############################ CELL CHAT ANALYSIS of E4 #######################################

# identify overexpression of genes and interactions
cellchatE4 <- identifyOverExpressedGenes(cellchatE4)
cellchatE4 <- identifyOverExpressedInteractions(cellchatE4)

# # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.mouse)

#Compute the communication probability and infer cellular communication network
cellchatE4 <- computeCommunProb(cellchatE4)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchatE4 <- filterCommunication(cellchatE4, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
df.netE4 <- subsetCommunication(cellchatE4)
df.netE4 <- subsetCommunication(cellchatE4, sources.use = c(2,11), targets.use = c(16,19,21,23,24))
#df.net <- subsetCommunication(cellchat, signaling = c("IL10"))

#Infer the cell-cell communication at a signaling pathway level
cellchatE4 <- computeCommunProbPathway(cellchatE4)

#Calculate the aggregated cell-cell communication network
cellchatE4 <- aggregateNet(cellchatE4)

groupSize <- as.numeric(table(cellchatE4@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchatE4@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchatE4@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf("cellchat_E4.pdf", width = 20, height = 20)
mat <- cellchatE4@net$weight
par(mfrow = c(6,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

##################################### Marge of Cellchat Objects #################################################

# merge cellchat objects
object.list <- list(E2 = cellchatE2, E4 = cellchatE4)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object

pdf("E2_E4_heatmap.pdf", width = 10, height = 5)
gg1 + gg2
dev.off()

# interaction weights
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#for selected cells

#Differential number of interactions or interaction strength among different cell types
group.cellType <- c(rep("Mø cluster 1", 4), rep("Mø cluster 2", 4), rep("Stromal cells", 4), rep("CD8+ T cells", 4), rep("CD4+ T cells", 4))
group.cellType <- factor(group.cellType, levels = c("Mø cluster 1", "Mø cluster 2", "Stromal cells", "CD8+ T cells", "CD4+ T cells"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# visualize
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#differential number of interactions or interaction strength
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)



