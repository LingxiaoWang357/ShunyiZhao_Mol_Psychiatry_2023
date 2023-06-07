### Load Packages -----------------------------------------------------------
library(Seurat)
library(eoffice)


### Preprocessing Tabula Muris data -----------------------------------------

### Seurat Objects for the Tabula Muris database are available at:
### https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263

# Update old Seurat objects
new_SeuratObj_list <- list()

organ_name <- vector()

obj_list <- list.files(path = "./Old_seurat_obj/")
for (i in 1:20) {
  print(obj_list[i])
  load(obj_list[i])
  print(tiss@version)
  new_SeuratObj_list[[i]] <- UpdateSeuratObject(tiss)
  organ_name[i] <- new_SeuratObj_list[[i]]@project.name
  print(new_SeuratObj_list[[i]]@version)
}

names(new_SeuratObj_list) <- organ_name

# Save updated Seurat Objects
saveRDS(new_SeuratObj_list, file = "updated_SeuratObj_all.RDS")


### Data processing ---------------------------------------------------------

# Merge all data into one Seurat Object

all_organ_Seurat <- merge(new_SeuratObj_list[[1]], y = new_SeuratObj_list[[2]], project = "Shunyi_Tabula")

for (i in 3:20) {
  all_organ_Seurat <- merge(all_organ_Seurat, y = new_SeuratObj_list[[i]], project = "Shunyi_Tabula")
}

unique(all_organ_Seurat@meta.data$orig.ident)

# Save merged Seurat Object
saveRDS(all_organ_Seurat, file = "all_organ_Seurat.RDS")

# Normalize data into Log(CPM) and store in assay$RNA@data
all_organ_Seurat <- NormalizeData(all_organ_Seurat, normalization.method = "LogNormalize", scale.factor = 1e6) #Log(CPM)

# Identify highly variable features
all_organ_Seurat <- FindVariableFeatures(all_organ_Seurat, selection.method = "vst", nfeatures = 2000)

# Scale data
all.genes <- rownames(all_organ_Seurat)

all_organ_Seurat <- ScaleData(all_organ_Seurat, features = all.genes)

# Linear dimensional reduction
all_organ_Seurat <- RunPCA(all_organ_Seurat, features = VariableFeatures(object = all_organ_Seurat))

# tSNE reduction
all_organ_Seurat <- RunTSNE(all_organ_Seurat, dims = 1:50)

Idents(all_organ_Seurat) <- "orig.ident"

p1 <- DimPlot(all_organ_Seurat, reduction = "tsne",pt.size = 0.5)

p2 <- FeaturePlot(all_organ_Seurat, "Cx3cr1", slot = "data", reduction = "tsne", cols = c("grey","red"), max.cutoff = 1,pt.size = 0.5)

# Subset Cx3cr1>0 cells
all_organ_Seurat_Cx3cr1 <- subset(x = all_organ_Seurat, subset = Cx3cr1 > 0)

Idents(all_organ_Seurat_Cx3cr1) <- "orig.ident"

unique(Idents(all_organ_Seurat_Cx3cr1))

all_organ_Seurat_Cx3cr1 <- RunTSNE(all_organ_Seurat_Cx3cr1, dims = 1:10)

DimPlot(all_organ_Seurat_Cx3cr1, reduction = "tsne")

saveRDS(all_organ_Seurat_Cx3cr1, file = "all_organ_Seurat_Cx3cr1.RDS")

all_organ_Seurat_Cx3cr1 <- RunUMAP(all_organ_Seurat_Cx3cr1, dims = 1:10)

p3 <- DimPlot(all_organ_Seurat_Cx3cr1, reduction = "umap")

### Visualization -----------------------------------------------------------

### Plot liver
Idents(new_SeuratObj_list$Liver) <- "cell_ontology_class"
Idents(new_SeuratObj_list$Liver)
#Normalize data into Log(CPM) and store in assay$RNA@data
new_SeuratObj_list$Liver <- NormalizeData(new_SeuratObj_list$Liver, normalization.method = "LogNormalize", scale.factor = 1e6)
#Identify highly variable features
new_SeuratObj_list$Liver <- FindVariableFeatures(new_SeuratObj_list$Liver, selection.method = "vst", nfeatures = 2000)
#Scale data
all.genes <- rownames(new_SeuratObj_list$Liver)
new_SeuratObj_list$Liver <- ScaleData(new_SeuratObj_list$Liver, features = all.genes)
#Perform linear dimensional reduction
new_SeuratObj_list$Liver <- RunPCA(new_SeuratObj_list$Liver, features = VariableFeatures(object = new_SeuratObj_list$Liver))
#Perfrom non-linear dimentional reduction
new_SeuratObj_list$Liver <- RunTSNE(new_SeuratObj_list$Liver, dims = 1:10)
p4 <- DimPlot(new_SeuratObj_list$Liver, reduction = "tsne")
p5 <- FeaturePlot(new_SeuratObj_list$Liver, "Cx3cr1", slot = "data", reduction = "tsne", cols = c("grey","red"), max.cutoff = 1)
new_SeuratObj_list$Liver <- RunUMAP(new_SeuratObj_list$Liver, dims = 1:10)
DimPlot(new_SeuratObj_list$Liver, reduction = "umap")

### Plot lung
Idents(new_SeuratObj_list$Lung) <- "cell_ontology_class"
Idents(new_SeuratObj_list$Lung)
#Normalize data into Log(CPM) and store in assay$RNA@data
new_SeuratObj_list$Lung <- NormalizeData(new_SeuratObj_list$Lung, normalization.method = "LogNormalize", scale.factor = 1e6)
#Identify highly variable features
new_SeuratObj_list$Lung <- FindVariableFeatures(new_SeuratObj_list$Lung, selection.method = "vst", nfeatures = 2000)
#Scale data
all.genes <- rownames(new_SeuratObj_list$Lung)
new_SeuratObj_list$Lung <- ScaleData(new_SeuratObj_list$Lung, features = all.genes)
#Perform linear dimensional reduction
new_SeuratObj_list$Lung <- RunPCA(new_SeuratObj_list$Lung, features = VariableFeatures(object = new_SeuratObj_list$Lung))
#Perfrom non-linear dimentional reduction
new_SeuratObj_list$Lung <- RunTSNE(new_SeuratObj_list$Lung, dims = 1:10)
p6 <- DimPlot(new_SeuratObj_list$Lung, reduction = "tsne")
p7 <- FeaturePlot(new_SeuratObj_list$Lung, "Cx3cr1", slot = "data", reduction = "tsne", cols = c("grey","red"), max.cutoff = 1)
new_SeuratObj_list$Lung <- RunUMAP(new_SeuratObj_list$Lung, dims = 1:10)
DimPlot(new_SeuratObj_list$Lung, reduction = "umap")


### Plot Kidney
Idents(new_SeuratObj_list$Kidney) <- "cell_ontology_class"
Idents(new_SeuratObj_list$Kidney)
#Normalize data into Log(CPM) and store in assay$RNA@data
new_SeuratObj_list$Kidney <- NormalizeData(new_SeuratObj_list$Kidney, normalization.method = "LogNormalize", scale.factor = 1e6)
#Identify highly variable features
new_SeuratObj_list$Kidney <- FindVariableFeatures(new_SeuratObj_list$Kidney, selection.method = "vst", nfeatures = 2000)
#Scale data
all.genes <- rownames(new_SeuratObj_list$Kidney)
new_SeuratObj_list$Kidney <- ScaleData(new_SeuratObj_list$Kidney, features = all.genes)
#Perform linear dimensional reduction
new_SeuratObj_list$Kidney <- RunPCA(new_SeuratObj_list$Kidney, features = VariableFeatures(object = new_SeuratObj_list$Kidney))
#Perfrom non-linear dimentional reduction
new_SeuratObj_list$Kidney <- RunTSNE(new_SeuratObj_list$Kidney, dims = 1:10)
p8 <- DimPlot(new_SeuratObj_list$Kidney, reduction = "tsne")
p9 <- FeaturePlot(new_SeuratObj_list$Kidney, "Cx3cr1", slot = "data", reduction = "tsne", cols = c("grey","red"), max.cutoff = 1)
new_SeuratObj_list$Kidney <- RunUMAP(new_SeuratObj_list$Kidney, dims = 1:10)
DimPlot(new_SeuratObj_list$Kidney, reduction = "umap")


### Plot spleen
Idents(new_SeuratObj_list$Spleen) <- "cell_ontology_class"
Idents(new_SeuratObj_list$Spleen)
#Normalize data into Log(CPM) and store in assay$RNA@data
new_SeuratObj_list$Spleen <- NormalizeData(new_SeuratObj_list$Spleen, normalization.method = "LogNormalize", scale.factor = 1e6)
#Identify highly variable features
new_SeuratObj_list$Spleen <- FindVariableFeatures(new_SeuratObj_list$Spleen, selection.method = "vst", nfeatures = 2000)
#Scale data
all.genes <- rownames(new_SeuratObj_list$Spleen)
new_SeuratObj_list$Spleen <- ScaleData(new_SeuratObj_list$Spleen, features = all.genes)
#Perform linear dimensional reduction
new_SeuratObj_list$Spleen <- RunPCA(new_SeuratObj_list$Spleen, features = VariableFeatures(object = new_SeuratObj_list$Spleen))
#Perfrom non-linear dimentional reduction
new_SeuratObj_list$Spleen <- RunTSNE(new_SeuratObj_list$Spleen, dims = 1:10)
p10 <- DimPlot(new_SeuratObj_list$Spleen, reduction = "tsne")
p11 <- FeaturePlot(new_SeuratObj_list$Spleen, "Cx3cr1", slot = "data", reduction = "tsne", cols = c("grey","red"), max.cutoff = 1)
new_SeuratObj_list$Spleen <- RunUMAP(new_SeuratObj_list$Spleen, dims = 1:10)
DimPlot(new_SeuratObj_list$Spleen, reduction = "umap")
# Single cell heatmap of feature expression
DoHeatmap(new_SeuratObj_list$Spleen, features = "Cx3cr1", )


### Plot bone marrow
Idents(new_SeuratObj_list$Marrow) <- "cell_ontology_class"
Idents(new_SeuratObj_list$Marrow)
#Normalize data into Log(CPM) and store in assay$RNA@data
new_SeuratObj_list$Marrow <- NormalizeData(new_SeuratObj_list$Marrow, normalization.method = "LogNormalize", scale.factor = 1e6)
#Identify highly variable features
new_SeuratObj_list$Marrow <- FindVariableFeatures(new_SeuratObj_list$Marrow, selection.method = "vst", nfeatures = 2000)
#Scale data
all.genes <- rownames(new_SeuratObj_list$Marrow)
new_SeuratObj_list$Marrow <- ScaleData(new_SeuratObj_list$Marrow, features = all.genes)
#Perform linear dimensional reduction
new_SeuratObj_list$Marrow <- RunPCA(new_SeuratObj_list$Marrow, features = VariableFeatures(object = new_SeuratObj_list$Marrow))
#Perform non-linear dimensional reduction
new_SeuratObj_list$Marrow <- RunTSNE(new_SeuratObj_list$Marrow, dims = 1:20)
p12 <- DimPlot(new_SeuratObj_list$Marrow, reduction = "tsne", pt.size = 1)
p13 <- FeaturePlot(new_SeuratObj_list$Marrow, "Cx3cr1", slot = "data", reduction = "tsne", cols = c("grey","red"), max.cutoff = 1, pt.size = 1)


### Export plots
p1+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p2+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p3+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p4+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p5+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p6+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p7+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p8+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p9+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p10+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p11+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p12+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)
p13+NoLegend()+NoAxes()+FontSize(
  x.text = 0,
  y.text = 0,
  x.title = 0,
  y.title = 0,
  main = 0)


### Get LogCPM value of Cx3cr1 from all tissue----------------------------------

for (p in 1:20) {
  cell_ontology <- unique(Idents(new_SeuratObj_list[[p]])) #Get tissue name
  for (q in 1:length(cell_ontology)){
    if (is.na(cell_ontology[q]) == F) {
      mtx <- as.matrix(GetAssayData(object = subset(new_SeuratObj_list[[p]], 
                                                    idents = cell_ontology[q],),
                                    slot = "data"))["Cx3cr1",]
      a <- paste0(paste(names(new_SeuratObj_list[p]), cell_ontology[q], "Cx3cr1_LogCPM", sep = "_"),".csv")
      write.csv(mtx, a)
    }
  }
}

