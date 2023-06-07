### random seed -------------------------------------------------------------
set.seed(123)

### load library ------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(harmony)

### save image --------------------------------------------------------------
save.image(file = "project_image.RData")

### load data ---------------------------------------------------------------

### Data are available at:
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150169

Ctrl1 <- ReadMtx(
  mtx="../GSM4525522_Ctrl-1_matrix.mtx.gz",
  cells = "../GSM4525522_Ctrl-1_barcodes.tsv.gz",
  features = "../GSM4525522_Ctrl-1_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
Ctrl1 <- CreateSeuratObject(counts = Ctrl1, min.cells = 10, min.features = 200)
Ctrl1[["Condition"]] <-  c('Ctrl')
Ctrl1[["Batch"]] <-  c('Batch-1')
Ctrl1[["sample_name"]] <-  c('Ctrl-1')
Ctrl1[["percent_mt"]] <- PercentageFeatureSet(Ctrl1, pattern = "^mt-")
Ctrl1[["percent_ribo"]] <- PercentageFeatureSet(Ctrl1, pattern = "^Rp[ls]")


Ctrl2 <- ReadMtx(
  mtx="../GSM4525523_Ctrl-2_matrix.mtx.gz",
  cells = "../GSM4525523_Ctrl-2_barcodes.tsv.gz",
  features = "../GSM4525523_Ctrl-2_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
Ctrl2 <- CreateSeuratObject(counts = Ctrl2, min.cells = 10, min.features = 200)
Ctrl2[["Condition"]] <-  c('Ctrl')
Ctrl2[["Batch"]] <-  c('Batch-2')
Ctrl2[["sample_name"]] <-  c('Ctrl-2')
Ctrl2[["percent_mt"]] <- PercentageFeatureSet(Ctrl2, pattern = "^mt-")
Ctrl2[["percent_ribo"]] <- PercentageFeatureSet(Ctrl2, pattern = "^Rp[ls]")


Ctrl3 <- ReadMtx(
  mtx="../GSM4525524_Ctrl-3_matrix.mtx.gz",
  cells = "../GSM4525524_Ctrl-3_barcodes.tsv.gz",
  features = "../GSM4525524_Ctrl-3_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
Ctrl3 <- CreateSeuratObject(counts = Ctrl3, min.cells = 10, min.features = 200)
Ctrl3[["Condition"]] <-  c('Ctrl')
Ctrl3[["Batch"]] <-  c('Batch-2')
Ctrl3[["sample_name"]] <-  c('Ctrl-3')
Ctrl3[["percent_mt"]] <- PercentageFeatureSet(Ctrl3, pattern = "^mt-")
Ctrl3[["percent_ribo"]] <- PercentageFeatureSet(Ctrl3, pattern = "^Rp[ls]")


D0_1 <- ReadMtx(
  mtx="../GSM4525525_D0-1_matrix.mtx.gz",
  cells = "../GSM4525525_D0-1_barcodes.tsv.gz",
  features = "../GSM4525525_D0-1_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D0_1 <- CreateSeuratObject(counts = D0_1, min.cells = 10, min.features = 200)
D0_1[["Condition"]] <-  c('D0')
D0_1[["Batch"]] <-  c('Batch-1')
D0_1[["sample_name"]] <-  c('D0-1')
D0_1[["percent_mt"]] <- PercentageFeatureSet(D0_1, pattern = "^mt-")
D0_1[["percent_ribo"]] <- PercentageFeatureSet(D0_1, pattern = "^Rp[ls]")


D0_2 <- ReadMtx(
  mtx="../GSM4525526_D0-2_matrix.mtx.gz",
  cells = "../GSM4525526_D0-2_barcodes.tsv.gz",
  features = "../GSM4525526_D0-2_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D0_2 <- CreateSeuratObject(counts = D0_2, min.cells = 10, min.features = 200)
D0_2[["Condition"]] <-  c('D0')
D0_2[["Batch"]] <-  c('Batch-2')
D0_2[["sample_name"]] <-  c('D0-2')
D0_2[["percent_mt"]] <- PercentageFeatureSet(D0_2, pattern = "^mt-")
D0_2[["percent_ribo"]] <- PercentageFeatureSet(D0_2, pattern = "^Rp[ls]")


D0_3 <- ReadMtx(
  mtx="../GSM4525527_D0-3_matrix.mtx.gz",
  cells = "../GSM4525527_D0-3_barcodes.tsv.gz",
  features = "../GSM4525527_D0-3_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D0_3 <- CreateSeuratObject(counts = D0_3, min.cells = 10, min.features = 200)
D0_3[["Condition"]] <-  c('D0')
D0_3[["Batch"]] <-  c('Batch-2')
D0_3[["sample_name"]] <-  c('D0-3')
D0_3[["percent_mt"]] <- PercentageFeatureSet(D0_3, pattern = "^mt-")
D0_3[["percent_ribo"]] <- PercentageFeatureSet(D0_3, pattern = "^Rp[ls]")


### data integration --------------------------------------------------------
seurat_all <- merge(Ctrl1, y = c(Ctrl2, Ctrl3, D0_1, D0_2, D0_3))
VlnPlot(seurat_all, "percent_mt")
seurat_filter <- subset(seurat_all, percent_mt < 10)
VlnPlot(seurat_filter, "percent_mt")
seurat_filter <- seurat_filter %>% 
  NormalizeData(.,scale.factor = 1e6) %>% 
  ScaleData(.) %>%
  FindVariableFeatures(.) %>% 
  RunPCA(.)

### harmony integration
Idents(seurat_filter) <- "sample_name"
seurat_filter <- seurat_filter %>% 
  RunHarmony(c("sample_name","Batch"), plot_convergence = TRUE)
seurat_filter <- RunUMAP(seurat_filter, reduction = "harmony", dims = 1:30)
Idents(seurat_filter) <- "Condition"
DimPlot(seurat_filter, split.by = "Condition")
DimPlot(seurat_filter)
FeaturePlot(seurat_filter, c("Ptprc", "Itgam", "Cx3cr1", "Ly6g"))
FeaturePlot(seurat_filter, c("Tmem119","Mrc1","Ccr2", "P2ry12"))
seurat_filter <- FindNeighbors(seurat_filter,reduction = "harmony", dims = 1:30)
seurat_filter <- FindClusters(seurat_filter, resolution = 0.1)
DimPlot(seurat_filter, label = T)


### Analyze microglia and related populations -------------------------------
seurat_0123 <- subset(seurat_filter, idents = c(0,1,2,3))
seurat_0123 <- seurat_0123 %>%
  NormalizeData(.,scale.factor = 1e6) %>% 
  ScaleData(.) %>%
  FindVariableFeatures(.) %>% 
  RunPCA(.)
seurat_0123 <- seurat_0123 %>% 
  RunHarmony(c("sample_name","Batch"), plot_convergence = TRUE)
seurat_0123 <- RunUMAP(seurat_0123, reduction = "harmony", dims = 1:40)
DimPlot(seurat_0123, split.by = "Condition")
DimPlot(seurat_0123)
FeaturePlot(seurat_0123, c("Ptprc", "Itgam", "Cx3cr1", "Ly6g"))
FeaturePlot(seurat_0123, c("Tmem119", "P2ry12", "Mertk", "Hexb", "Sall1", "Mrc1", "Ms4a7", "Pf4", "Ccr2"))
seurat_0123 <- FindNeighbors(seurat_0123,reduction = "harmony", dims = 1:40)
seurat_0123 <- FindClusters(seurat_0123, resolution = 0.1)
DimPlot(seurat_0123, label = T) 
DimPlot(seurat_0123, label = T, split.by = "Condition")
FeaturePlot(seurat_0123, c("Lgals3", "Tmem119"), split.by = "Condition")


### Subset microglia 0124 ---------------------------------------------------
seurat_microglia <- subset(seurat_0123, 
                           idents = c(0,1,2,4))
seurat_microglia <- seurat_microglia %>%
  NormalizeData(.,scale.factor = 1e6) %>% 
  ScaleData(.) %>%
  FindVariableFeatures(.) %>% 
  RunPCA(.)
seurat_microglia <- seurat_microglia %>% 
  RunHarmony(c("sample_name","Batch"), plot_convergence = TRUE)
seurat_microglia <- RunUMAP(seurat_microglia, reduction = "harmony", dims = 1:40)
DimPlot(seurat_microglia, split.by = "Condition")
DimPlot(seurat_microglia)
FeaturePlot(seurat_microglia, c("Ptprc", "Itgam", "Cx3cr1", "Ly6g"))
FeaturePlot(seurat_microglia, c("Tmem119", "P2ry12", "Mertk", "Hexb", "Sall1", "Mrc1", "Ms4a7", "Pf4", "Ccr2"))
seurat_microglia <- FindNeighbors(seurat_microglia,reduction = "harmony", dims = 1:40)
seurat_microglia <- FindClusters(seurat_microglia, resolution = 0.1)
DimPlot(seurat_microglia, label = T)
DimPlot(seurat_microglia, label = T, split.by = "Condition")
FeaturePlot(seurat_microglia, c("Lgals3", "Tmem119"), split.by = "Condition")


### subset Lgals3+ microglia ------------------------------------------------
seurat_Lgals3s <- subset(seurat_microglia, Lgals3 > 1)
Idents(seurat_Lgals3s) <- "Condition"
DimPlot(seurat_Lgals3s, split.by = "Condition")
FeaturePlot(seurat_Lgals3s, "Tmem119",split.by = "Condition", cols = c("grey", "red"), max.cutoff = 1)
VlnPlot(seurat_Lgals3s, "Tmem119", split.by = "Condition")



### Get LogCPM of Tmem119 ---------------------------------------------------
Tmem119_expression <- as.data.frame(seurat_Lgals3s@assays$RNA@data["Tmem119",])

Tmem119_expression <- cbind(rownames(Tmem119_expression), Tmem119_expression)
colnames(Tmem119_expression) <- c("CellID","Tmem119")

Tmem119_expression_name <- as.data.frame(cbind(seurat_Lgals3s@meta.data$Condition, seurat_Lgals3s@assays$RNA@data@Dimnames[[2]]))
colnames(Tmem119_expression_name) <- c("Condition","CellID")

write.csv(left_join(Tmem119_expression, Tmem119_expression_name,"CellID"), "Lgals3+_Microglia_Tmem119_Expression.csv")
