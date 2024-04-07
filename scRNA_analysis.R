setwd("C:/Users/Jacky/Desktop/PDAC")
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(harmony)
library(DoubletFinder)
library(glmGamPoi)
library(Rcpp)
library(openxlsx)
library(tibble)
library(here)
library(phylogram)
library(ggrastr)
library(ggbeeswarm)
library(beeswarm)
library(vipor)
library(Cairo)
library(data.table)
library(SCP)
library(BiocParallel)

seu_obj <- readRDS("./scRNA/RDS/Total_Harmony.rds")

sample <- unique(seu_obj@meta.data$orig.ident)

Idents(seu_obj) <- "orig.ident"
for (sample_name in sample) {
  extracted_sample <- subset(seu_obj, idents = sample_name)  # 提取样本
  # 将提取的样本对象保存到与样本名称相同的对象名中
  assign(sample_name, extracted_sample)
}

seu_obj <- merge(P20181121,y = c(PACA.P20181128,PACA.P20190225,PACA.P20190306,P1,P2,P3,P4,P5,PDAC_WY_01,PDSC_WY_02,PDAC_WY_03,P03,P04,P05,P07,P10,P12,P13,P14,P15,P19,P20,P22,P23,P26,Case1,Case2,Case3,MDA1,MDA2,AdjNorm_TISSUE_1,AdjNorm_TISSUE_2,BL5.h5,BL10.h5,BL11.h5,BL12.h5,TIPC249,TIPC262,TIPC282,TIPC301,TIPC309,TIPC416,TIPC418,TIPC432))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu_obj <- CellCycleScoring(seu_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

seu_obj1 <- readRDS("Harmony_seurat.rds")
all_objects <- ls()
objects_to_remove <- grep("^CD8", all_objects, value = TRUE)
rm(list = objects_to_remove)

selected_metadata_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
new_metadata <- seu_obj@meta.data[, selected_metadata_columns]
seu_obj@meta.data <- new_metadata

cell_counts <- table(seu_obj@meta.data$orig.ident)
samples_to_keep <- names(cell_counts[cell_counts >= 500])
seu_obj <- subset(seu_obj, subset = orig.ident %in% samples_to_keep)
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")

seu_obj[['SCT']] <- NULL
seu_obj[['Seurat']] <- NULL
#saveRDS(seu_obj,"merge.rds")

table(seu_obj@meta.data$orig.ident)

seu_obj <- Standard_SCP(srt = seu_obj)
Idents(seu_obj) <- "orig.ident"
CellDimPlot(srt = seu_obj, group.by = "orig.ident")

CellDimPlot(
  srt = seu_obj, group.by = "Harmonyclusters", reduction = "HarmonyUMAP2D",
  title = "Seurat", theme_use = "theme_blank",label = T
)

CellDimPlot(
  srt = seu_obj, group.by = "orig.ident", reduction = "HarmonyUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)

###3. ϸ??ע?? ʹ??2022??NC scTypeע??????
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

#Cell type assignment??ע??ǰ??׼????????reference?Ĺ??? #source("./gene_sets_prepare.R") ???ؼ???
source("C:/Users/Jacky/Desktop/PDAC/scRNA/annotation/gene_sets_prepare.R")
# load cell type annotation function
source("C:/Users/Jacky/Desktop/PDAC/scRNA/annotation/sctype_score_.R")

# DB file??ע?Ͳο????󣬿??԰??ø?ʽ?????Լ???ע?ͻ??򼯣????иĽ?
db_ = "C:/Users/Jacky/Desktop/PDAC/scRNA/annotation/ScTypeDB_full.xlsx"
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets ????ע???ź?reference?Ľ?��
gs_list = gene_sets_prepare(db_, tissue)

#Finally, let's assign cell types to each cluster: ?Ծ?????ÿ??cluster????ע??
es.max = sctype_score(scRNAseqData = seu_obj[["RNA"]]@counts, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seu_obj@meta.data$Harmonyclusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu_obj@meta.data[seu_obj@meta.data$Harmonyclusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu_obj@meta.data$Harmonyclusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown" ????��ע????Ϣ????Ϊunknown
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

seu_obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu_obj@meta.data$customclassif[seu_obj@meta.data$Harmonyclusters == j] = as.character(cl_type$type[1])
}

table(seu_obj@meta.data$Harmonyclusters.y)

selected_metadata_columns <- c("orig.ident", "nCount_RNA", "nFeature_RNA","pMT","Harmonyclusters")
new_metadata <- seu_obj@meta.data[, selected_metadata_columns]
seu_obj@meta.data <- new_metadata

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
#metatable <- read_excel("./metadata/metadata.xlsx")
#row_names <- rownames(seu_obj@meta.data)
#seu_obj@meta.data$cell_id <- row_names  
#metadata <- FetchData(seu_obj, "cell_id")
#metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
#metadata$Harmonyclusters <- metadata$cell_id #??ȡÿ??cell??????ID
#metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
#rownames(metadata) <- metadata$cell_id
#seu_obj <- AddMetaData(seu_obj, metadata = metadata)

#Idents(seu_obj) <- seu_obj@meta.data$Harmonyclusters.y
#clusters_to_keep <- c( "2" ,"3" , "5" , "6" , "7" , "8" , "9" , "10", "12" , "13" , "14", "15")
#seu_obj <- subset(seu_obj, subset = Harmonyclusters.y %in% clusters_to_keep)
#seu_obj <- Integration_SCP(srtMerge = seu_obj, batch = "orig.ident", integration_method = "Harmony")

Idents(seu_obj) <- seu_obj@meta.data$Harmonyclusters
seu_obj.markers <- FindAllMarkers(seu_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seu_obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(data.frame(ID=rownames(seu_obj.markers),seu_obj.markers), file = "cluster_marker_genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

Tm_mainmarkers <- c( "NKTR", "CD1D")
for (i in seq_along(Tm_mainmarkers)) {
  FeaturePlot(seu_obj, features = Tm_mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10),raster=FALSE)
  ggsave2(paste0("FeaturePlot_mainmarkers_", Tm_mainmarkers[i], ".svg"), path = "./Tm_results", width = 20, height = 30, units = "cm")
}

cell_counts <- table(seu_obj@meta.data$Harmonyclusters)
samples_to_keep <- names(cell_counts[cell_counts >= 500])
seu_obj <- subset(seu_obj, subset = Harmonyclusters %in% samples_to_keep)

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/celltype.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "Harmonyclusters")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "Harmonyclusters") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

CellDimPlot(
  srt = seu_obj, group.by = "celltype", reduction = "HarmonyUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)

metadata <- seu_obj@meta.data
write.table(metadata, file = "metadata_fastMNN.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
Tm_mainmarkers <- c( "TSC22D1", "PRDM1", "TRPS1", "ZBED2")
for (i in seq_along(Tm_mainmarkers)) {
  FeaturePlot(seu_obj, features = Tm_mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10),raster=FALSE)
  ggsave2(paste0("FeaturePlot_mainmarkers_", Tm_mainmarkers[i], ".svg"), path = "./Tm_results", width = 20, height = 30, units = "cm")
}

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/Neo_T.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "cell_id")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

Idents(seu_obj) <- seu_obj@meta.data$Predicted
DimPlot(seu_obj, cells.highlight=WhichCells(seu_obj, idents="1"),
        cols.highlight = "#62D6F6" , pt.size = 0.00001) #ָ??????ϸ??????ɫ

DimPlot(seu_obj, cells.highlight=WhichCells(seu_obj, idents="1"))

library("Seurat")
library("Matrix")
library("readxl")
library("dplyr")
library(easypackages)

libraries("Seurat", "Matrix", "readxl","RColorBrewer",'Rmagic',
          'patchwork','dplyr','viridis','ggplot2','pals','SeuratWrappers')

library(cols4all)
mycol <- c4a('vivid',12)

clust.cp.graded = unique(c(stepped3(16),stepped(20),stepped2(20)))
CLcp=clust.cp.graded[c(11,9,25,3,1,18,26,27,28,8,6,5,19)]

clust.cp=CLcp
clust.cp = c('#A1D99B','#31A354')

Idents(seu_obj) <- seu_obj@meta.data$Predicted
ids.cluster.library.AllData = as.data.frame(table(Idents(seu_obj), seu_obj@meta.data$orig.ident))
colnames(ids.cluster.library.AllData) = c('ID','Library','CellCount')

ids.cluster.library.AllData$ID <- factor(ids.cluster.library.AllData$ID, levels = c(0, 1))  #调换顺序

Fig1D=
  ggplot(ids.cluster.library.AllData, aes(fill=ID, y= CellCount,
                                          x=Library)) +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),
                        x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp)+
  theme(axis.text.x = element_text(#face="bold", color="#993333", 
    size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= (CellCount),
                         x=as.integer(Library)),
            position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                        x=(Library)),
           position="fill", stat="identity", width = 0.5)+
  ggtitle("Distribution of cell types in time and space")

print(Fig1D)

Fig1D <- ggplot(ids.cluster.library.AllData, aes(fill=ID, y=CellCount, x=Library)) +  
  geom_bar(position="fill", stat="identity", width = 0.5) +  
  scale_fill_manual(values = clust.cp) +  
  theme(axis.text.x = element_text(size=8, angle=-45, hjust=0, vjust = 0.5),  
        panel.background = element_rect(fill = "white"),  
        plot.background = element_rect(fill = "white")) +  
  geom_area(mapping = aes(y=CellCount, x=as.integer(Library)),  
            position="fill", stat="identity", alpha=0.2, size=.5, colour="white") +  
  ggtitle("Distribution of cell types in time and space")  

print(Fig1D)

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/CD4_annotation.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "cell_id")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

saveRDS(seu_obj,"45sample_merge+CD4_annotation.rds")
seu_obj <- readRDS("./scRNA/RDS/45sample_merge+CD4_filtered.rds")

Idents(seu_obj) <- "celltypes"

sample2 <- c("CD8", "EOMES+ NK-like", "GZMK+ Tem", "IL7R+ Tm", "Proliferating", "TC17", "Temra", "terminal Tex", "TXK+ NK-like", "ZNF683+CXCR6+ Trm")

for (sample_name in sample2) {
  extracted_sample <- subset(seu_obj, idents = sample_name)  # 提取样本
  # 将提取的样本对象保存到与样本名称相同的对象名中
  assign(sample_name, extracted_sample)
}

seu_obj <- merge(CD8,y = c(EOMES+ NK-like,GZMK+ Tem,IL7R+ Tm,Proliferating,TC17,Temra,terminalTex,TXK+ NK-like,ZNF683+CXCR6+ Trm))

table(seu_obj@meta.data$celltypes)

endothelial_mainmarkers <- c("CXCL13","CTLA4","PDCD1","LAG3")

endothelial_mainmarkers <- c("ZNF683")
for (i in seq_along(endothelial_mainmarkers)) {
  FeaturePlot(seu_obj, features = endothelial_mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10),raster=FALSE)
  ggsave2(paste0("FeaturePlot_mainmarkers_", endothelial_mainmarkers[i], ".svg"), path = "./", width = 20, height = 30, units = "cm")
}

CellDimPlot(
  srt = seu_obj, group.by = "Harmonyclusters", reduction = "HarmonyUMAP2D",
  title = "Seurat", theme_use = "theme_blank", label = T
)

Idents(seu_obj) <- "pMT"

pMT_upper <- 15
seu_obj <- subset(seu_obj, subset = pMT < pMT_upper)

seu_obj <- readRDS("./scRNA/RDS/45sample_merge+CD4_filtered_Harmony.rds")

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/fastMNN_annotation.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "fastMNNclusters")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "fastMNNclusters") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

saveRDS(seu_obj,"45sample_merge+fastMNNclusters_annotation.rds")
seu_obj <- readRDS("./scRNA/RDS/45sample_merge+fastMNNclusters_annotation.rds")

CellDimPlot(
  srt = seu_obj, group.by = "celltype", reduction = "fastMNNUMAP2D",
  title = "Seurat", theme_use = "theme_blank", label = T
)

cols_to_remove <- c("Harmony_SNN_res.0.6", "Harmonyclusters", "celltypes")  
seu_obj@meta.data <- seu_obj@meta.data[, !(names(seu_obj@meta.data) %in% cols_to_remove)]

CellDimPlot(
  srt = seu_obj, group.by = "orig.ident", reduction = "fastMNNUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)

seu_obj@meta.data$orig.ident

CellDimPlot(
  srt = seu_obj, group.by = "fastMNNclusters", reduction = "fastMNNUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/Virus_T.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "cell_id")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

Idents(seu_obj) <- seu_obj@meta.data$virus
DimPlot(seu_obj, cells.highlight=WhichCells(seu_obj, idents="1"),
        cols.highlight = "#28C580" , pt.size = 0.00001) #ָ??????ϸ??????ɫ

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/Virus_T.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "cell_id")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

seu_obj@meta.data <- select(seu_obj@meta.data, -Label)

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")
metatable <- read_excel("./metadata/Neo_T.xlsx")
row_names <- rownames(seu_obj@meta.data)
seu_obj@meta.data$cell_id <- row_names  
metadata <- FetchData(seu_obj, "cell_id")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)
















