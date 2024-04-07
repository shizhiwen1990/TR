setwd("C:/Users/Jacky/Desktop/DEGs")

seu_obj <- readRDS("./RDS/annotation.rds")
#在Seurat对象的metadata中删除指定的列
seu_obj@meta.data <- select(seu_obj@meta.data, -cell_id)
seu_obj@meta.data <- select(seu_obj@meta.data, -HPV_Specific)
seu_obj@meta.data <- select(seu_obj@meta.data, -sample_id)

metatable <- read_excel("./metadata/Neo_metadata.xlsx")
metadata <- FetchData(seu_obj, "cell.name")
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x = metadata, y = metatable, by = "cell.name")
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)
View(seu_obj@meta.data)

CellDimPlot(
  srt = seu_obj, group.by = c("Phase"), reduction = "SeuratUMAP2D",
  title = "Seurat", theme_use = "theme_blank",label = TRUE
)

CellDimPlot(
  srt = seu_obj, group.by = c("T_cell_label"), reduction = "SeuratUMAP2D",
  title = "Seurat", theme_use = "theme_blank"
)

metatable <- read_excel("./metadata/patients_metadata.xlsx")
metadata <- FetchData(seu_obj, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
rownames(metadata) <- metadata$cell_id
seu_obj <- AddMetaData(seu_obj, metadata = metadata)
View(seu_obj@meta.data)

meatadata <- seu_obj@meta.data
write.table(data.frame(ID=rownames(meatadata),meatadata), file = "metadata.txt", sep = "\t", col.names = T, row.names = F, quote = F)

saveRDS(seu_obj,"annotation.rds")

setwd("C:/Users/Jacky/Desktop/HPV_scRNA")
library(scRNAtoolVis)
library(ggplot2)
library(tidyverse)
library(ggrepel)

setwd("C:/Users/Jacky/Desktop/PDAC/scRNA/DEGs")
cell_types <- unique(seu_obj@meta.data$celltype)  

for (cell_type in cell_types) {  
  # 根据celltype标识细胞  
  Idents(seu_obj) <- seu_obj@meta.data$celltype  
  
  # 提取特定细胞类型的子集  
  cell_subset <- subset(seu_obj, idents = cell_type)  
  
  # 根据T_cell_label重新标识细胞子集  
  Idents(cell_subset) <- cell_subset@meta.data$Label  
  
  # 查找差异表达基因  
  DEGs <- FindMarkers(cell_subset, ident.1 = "1", ident.2 = "0", min.pct = 0.25)  
  
  # 将差异表达基因保存到文件  
  output_file <- paste0(cell_type, "_DEGs.txt")  
  write.table(data.frame(ID = rownames(DEGs), DEGs), file = output_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)  
}

file_list <- list.files(pattern = "\\.txt$")

# 读取所有文件并合并
combined_data <- lapply(file_list, function(file) {
  # 读取文件
  data <- read.table(file, header = TRUE, sep = "\t")  # 假设文件是用制表符分隔的
  
  # 添加文件名作为新列
  data$cluster <- file
  
  # 返回修改后的数据框
  return(data)
}) %>% bind_rows()

names(combined_data)[names(combined_data) == "ID"] <- "gene" #修改ID名称为gene
write.table(combined_data, file = "combined_data.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

df1 <- read.csv("combined_data.csv", header = TRUE, encoding = 'UTF-8')

mycol <- c( "steelblue","#E64B357F","#4DBBD57F","#F98400","#DF6FA0","#00A0877F","#3C54887F","#F39B7F7F","#8491B47F","#E2D200", 'pink')

#添加其它参数修改标签,这里调整一下字体和文字大小:
jjVolcano(diffData = df1,
          tile.col = corrplot::COL2('RdBu', 15)[2:12],
          size  = 3.5,
          fontface = 'italic')

#修改点颜色:
jjVolcano(diffData = df1) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol)
















































