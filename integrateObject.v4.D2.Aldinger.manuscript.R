#  Integrating the D2 cerebellar organoid object with the Aldinger 
#  downsampled object.


library(dplyr)
library(Seurat)

# Versions
R.Version()$version.string
#[1] "R version 4.0.3 (2020-10-10)"

packageVersion('dplyr')
# [1] ‘1.0.2’

packageVersion('Seurat')
# [1] '4.3.0'

packageVersion('SeuratObject')
# [1] ‘4.1.3’


setwd("/project/tuannguy_229/202107 GEO processed data")

# Load Cerebellar organoid data
load(file='object.Cerebellar-only.D2.final.Robj') 
object
# An object of class Seurat 
# 67744 features across 15302 samples within 3 assays 
# Active assay: integrated (19732 features, 3000 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap


# Rename meta of integrated_snn_res...
object[["organoid_res.0.2"]] <- object@meta.data$integrated_snn_res.0.2
object[["organoid_res.0.8"]] <- object@meta.data$integrated_snn_res.0.8
object[["organoid_res.1.2"]] <- object@meta.data$integrated_snn_res.1.2

object@meta.data[13:16] <- NULL
head(object@meta.data)

# Add age column
object[['age']] <- rep(x = '2-month', times=length(object@meta.data$orig.ident))


# Load Aldinger's downsampled data
load(file='object.seurat.v4.fetal.downsample500.Aldinger.Robj')
down
# An object of class Seurat 
# 90171 features across 9790 samples within 3 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 2 other assays present: RNA, SCT
# 3 dimensional reductions calculated: pca, umap, tsne

# Remove some extra meta
down@meta.data[c(5:7,13:15,17:31)] <- NULL


# Set to SCT assay 
DefaultAssay(object) <- 'SCT'
DefaultAssay(down) <- 'SCT'

###############################################
# Integrate organoid and fetal Seurat objects #
###############################################
sample.list <- list(object, down)

# Both are already SCT normalization with CC.Difference score 


# Get all the unique gene names
all_features <- lapply(sample.list, row.names) %>% Reduce(intersect, .) 
length(all_features)
# [1] 16184

# Using same CC.Difference that was generated before 
sample.list <- lapply(X = sample.list, FUN = SCTransform, vars.to.regress = 'CC.Difference')


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)

sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT", anchor.features = features)

object.integrated.cca <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = all_features)
object.integrated.cca
# An object of class Seurat 
# 109876 features across 25092 samples within 3 assays 
# Active assay: integrated (15570 features, 2969 variable features)
# 2 other assays present: RNA, SCT


#####################################
# Run for PCA and UMAP and Clusters #
#####################################
object.integrated.cca <- RunPCA(object.integrated.cca, verbose = T)

object.integrated.cca <- RunUMAP(object.integrated.cca, reduction = "pca", dims = 1:30)

DefaultAssay(object.integrated.cca) <- "integrated"

object.integrated.cca <- FindNeighbors(object.integrated.cca, reduction = "pca", dims = 1:30)

object.integrated.cca <- FindClusters(object.integrated.cca, resolution = c(0.2,0.8,1.2))


DimPlot(object.integrated.cca)

################################################
# SAVE Integrated D2 organoid and Fetal object # 
################################################
save(object.integrated.cca, file='object.integrated.Cerebellar.D2.Aldinger.Robj') 






