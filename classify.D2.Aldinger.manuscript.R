#  Use the object generated from theintegration with Aldinger data
#  to label the clusters.

library(dplyr)
library(Seurat)


R.Version()$version.string
# [1] "R version 4.0.3 (2020-10-10)"

packageVersion('dplyr')
# [1] ‘1.0.2’

packageVersion('Seurat')
# [1] ‘4.3.0’

setwd("/project/tuannguy_229/202107 GEO processed data")

# Load integrated fetal and organoid data set
load(file='object.integrated.Cerebellar.D2.Aldinger.Robj') 
# object.integrated.cca
# An object of class Seurat 
# 109876 features across 25092 samples within 3 assays 
# Active assay: integrated (15570 features, 2969 variable features)
# 2 other assays present: RNA, SCT
# 2 dimensional reductions calculated: pca, umap

# Use the unbiased res.0.8 number 
Idents(object.integrated.cca) <- 'integrated_snn_res.0.8'

levels(Idents(object.integrated.cca))
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25"

################################################################
# With Alexander, we assigned these one-to-one classification. #
################################################################

Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(0))) <- 'Progenitors'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(1))) <- 'GCP'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(6))) <- 'Glia'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(2))) <- 'PC'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(4))) <- 'GC'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(3))) <- 'immature-iCN/\nimmature-PC'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(22))) <- 'iCN'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(5))) <- 'BG'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(8))) <- 'PIP'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(10,9,21))) <- 'Choroid'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(13))) <- 'Roof-Plate'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(11))) <- 'VZ'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(23))) <- 'Committed-OPC'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(15))) <- 'Brain-Stem'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(14,19))) <- 'Div-Choroid'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(16))) <- 'Meninges'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(18))) <- 'RL'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(17))) <- 'eCN/Unibrush'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(20))) <- 'OPC'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(7))) <- 'Endo'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(12))) <- 'Div-VZ'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(24))) <- 'Micro'
Idents(object = object.integrated.cca, cells = WhichCells(object.integrated.cca, idents = c(25))) <- 'MLI'

object.integrated.cca[["final.clusters"]] <- Idents(object.integrated.cca)
#DimPlot(object.integrated.cca, group.by = 'final.clusters', pt.size=1, label=T)


##########
## SAVE ##
##########
save(object.integrated.cca, file = "/Users/quadratolab2/Documents/Single Cell Data/10X USC/202107/object.integrated.cca.manuscipt.jan2024.Robj")


