CD8_2023STM <- Standard_SCP(srt = CD8_2023STM)
Idents(CD8_2023STM) <- "orig.ident"
CellDimPlot(srt = CD8_2023STM, group.by = "orig.ident")

sample1 <- unique(seu_obj@meta.data$orig.ident)
samples <- c("BL5.h5","BL10.h5","BL11.h5","BL12.h5","TIPC249","TIPC262","TIPC282","TIPC301","TIPC309","TIPC416","TIPC418","TIPC432")

for (sample in samples) {
  # 加载样本数据
  # 假设你的样本数据是以某种方式加载的，这里使用的是假设函数loadData()
  sample_data <- get(sample)  # 请替换为实际的数据加载方式
  
  # 应用相同的处理流程
  sample_data <- PercentageFeatureSet(sample_data, pattern = "^MT-", col.name = "pMT")
  sample_data <- SCTransform(sample_data, vars.to.regress = c("pMT", "nCount_RNA"), verbose = TRUE, variable.features.n = 6000)
  
  SCT <- sample_data@assays$SCT@scale.data
  SCT_df <- as.data.frame(SCT)
  SCT_df$ID <- rownames(SCT_df)
  SCT_df <- SCT_df[, c("ID", setdiff(names(SCT_df), "ID"))]
  
  # 将结果写入文件，每个样本一个文件
  file_name <- paste0(sample, "_SCT.txt")
  write.table(SCT_df, file = file_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}


TIPC418 <- PercentageFeatureSet(TIPC418, pattern = "^MT-", col.name = "pMT")
TIPC418 <- SCTransform(TIPC418, vars.to.regress = c("pMT","nCount_RNA"), verbose = TRUE, variable.features.n = 6000)

SCT <- TIPC418@assays$SCT@scale.data
SCT_df <- as.data.frame(SCT)
SCT_df$ID <- rownames(SCT_df)
SCT_df <- SCT_df[, c("ID", setdiff(names(SCT_df), "ID"))]
write.table(SCT_df, file = "TIPC418_SCT.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

TIPC432 <- PercentageFeatureSet(TIPC432, pattern = "^MT-", col.name = "pMT")
TIPC432 <- SCTransform(TIPC432, vars.to.regress = c("pMT","nCount_RNA"), verbose = TRUE, variable.features.n = 6000)

SCT <- TIPC432@assays$SCT@scale.data
SCT_df <- as.data.frame(SCT)
SCT_df$ID <- rownames(SCT_df)
SCT_df <- SCT_df[, c("ID", setdiff(names(SCT_df), "ID"))]
write.table(SCT_df, file = "TIPC432_SCT.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

TIPC416 <- PercentageFeatureSet(TIPC416, pattern = "^MT-", col.name = "pMT")
TIPC416 <- SCTransform(TIPC416, vars.to.regress = c("pMT","nCount_RNA"), verbose = TRUE, variable.features.n = 6000)

SCT <- TIPC416@assays$SCT@scale.data
SCT_df <- as.data.frame(SCT)
SCT_df$ID <- rownames(SCT_df)
SCT_df <- SCT_df[, c("ID", setdiff(names(SCT_df), "ID"))]
write.table(SCT_df, file = "TIPC416_SCT.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

table(seu_obj@meta.data$orig.ident)

row_names <- rownames(CD8_2023STM@meta.data)
CD8_2023STM@meta.data$cell_id <- row_names  
metadata <- CD8_2023STM@meta.data
write.table(metadata, file = "STM_metadata.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

setwd("C:/Users/Jacky/Desktop/SCT")
genes_file <- "filtered_gene_list.txt"
genes_data <- read.table(genes_file, header = TRUE, stringsAsFactors = FALSE)
genes_of_interest <- genes_data[, 1]  # 假设基因名是第一列

for (sample in samples) {
  # 读取每个样本的SCT文件
  file_name <- paste0(sample, "_SCT.txt")
  data <- read.table(file_name, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 筛选特定的行
  filtered_data <- data[data$ID %in% genes_of_interest,]
  
  # 对筛选后的数据进行标准化（z-score）
  # 忽略第一列（ID列），因为它不包含数值数据
  filtered_data[,-1] <- apply(filtered_data[,-1], 2, function(x) (x - mean(x)) / sd(x))
  
  # 保存标准化后的数据到新文件
  output_file_name <- paste0(sample, "_Filtered_SCT.txt")
  write.table(filtered_data, file = output_file_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

merged_data <- NULL

for (sample in samples) {
  # 读取每个样本的过滤后的SCT文件
  file_name <- paste0(sample, "_Filtered_SCT.txt")
  data <- read.table(file_name, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 如果是第一个样本，直接赋值给merged_data
  # 否则，合并到现有的数据框架
  if (is.null(merged_data)) {
    merged_data <- data
  } else {
    merged_data <- merge(merged_data, data, by = "ID", all = TRUE)
  }
}

data <- t(merged_data)
write.table(data.frame(ID=rownames(data),data), file = "merge.txt", sep = "\t", col.names = F, row.names = T, quote = F)

cell_counts <- table(seu_obj@meta.data$orig.ident)
samples_to_keep <- names(cell_counts[cell_counts >= 1000])
sample2 <- c("BL11.h5", "Case2", "Case3", "MDA1", "MDA2", "P03", "P10", "P14", "P15", "P19", "P20", "P20181121", "P22", "P26", "P4", "PACA.P20181128", "PACA.P20190306", "PDAC_WY_03", "PDSC_WY_02", "TIPC249", "TIPC262", "TIPC301", "TIPC309", "TIPC416", "TIPC418", "TIPC432", "TIPC282", "PDAC_WY_01")
seu_obj <- subset(seu_obj, subset = orig.ident %in% data)

sample2 <- c("P20181121", "PACA.P20181128", "PACA.P20190225", "PACA.P20190306", "P1", "P2", "P3", "P4", "P5", "PDAC_WY_01", "PDSC_WY_02", "PDAC_WY_03", "P03", "P04", "P05", "P07", "P10", "P12", "P13", "P14", "P15", "P19", "P20", "P22", "P23", "P26", "Case1", "Case2", "Case3", "MDA1", "MDA2", "AdjNorm_TISSUE_1", "AdjNorm_TISSUE_2", "BL5.h5", "BL10.h5", "BL11.h5", "BL12.h5")


for (sample_name in sample2) {
  extracted_sample <- subset(seu_obj, idents = sample_name)  # 提取样本
  # 将提取的样本对象保存到与样本名称相同的对象名中
  assign(sample_name, extracted_sample)
}

for (sample in sample2) {
  # 加载样本数据
  # 假设你的样本数据是以某种方式加载的，这里使用的是假设函数loadData()
  sample_data <- get(sample)  # 请替换为实际的数据加载方式
  
  # 应用相同的处理流程
  sample_data <- PercentageFeatureSet(sample_data, pattern = "^MT-", col.name = "pMT")
  sample_data <- SCTransform(sample_data, vars.to.regress = c("pMT", "nCount_RNA"), verbose = TRUE, variable.features.n = 20000)
  
  SCT <- sample_data@assays$SCT@scale.data
  SCT_df <- as.data.frame(SCT)
  SCT_df$ID <- rownames(SCT_df)
  SCT_df <- SCT_df[, c("ID", setdiff(names(SCT_df), "ID"))]
  
  # 将结果写入文件，每个样本一个文件
  file_name <- paste0(sample, "_SCT.txt")
  write.table(SCT_df, file = file_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

sample2 <- c("Case2", "Case3", "MDA1", "MDA2", "P03", "P14", "P15", "P19", "P20181121", "PACA.P20181128", "PACA.P20190306", "PDAC_WY_03", "PDSC_WY_02", "TIPC249", "TIPC262", "TIPC301", "TIPC309", "TIPC416", "TIPC418", "TIPC432", "TIPC282", "PDAC_WY_01")

for (sample in sample2) {
  # 读取每个样本的SCT文件
  file_name <- paste0(sample, "_SCT.txt")
  data <- read.table(file_name, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 筛选特定的行
  filtered_data <- data[data$ID %in% genes_of_interest,]
  
  # 对筛选后的数据进行标准化（z-score）
  # 忽略第一列（ID列），因为它不包含数值数据
  filtered_data[,-1] <- apply(filtered_data[,-1], 2, function(x) (x - mean(x)) / sd(x))
  
  # 保存标准化后的数据到新文件
  output_file_name <- paste0(sample, "_Filtered_SCT.txt")
  write.table(filtered_data, file = output_file_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

merged_data <- NULL

for (sample in sample2) {
  # 读取每个样本的过滤后的SCT文件
  file_name <- paste0(sample, "_Filtered_SCT.txt")
  data <- read.table(file_name, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 如果是第一个样本，直接赋值给merged_data
  # 否则，合并到现有的数据框架
  if (is.null(merged_data)) {
    merged_data <- data
  } else {
    merged_data <- merge(merged_data, data, by = "ID", all = TRUE)
  }
}

data1 <- t(merged_data)
write.table(data.frame(ID=rownames(data1),data1), file = "merge.txt", sep = "\t", col.names = F, row.names = T, quote = F)


























