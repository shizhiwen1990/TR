setwd("C:/Users/Jacky/Desktop/PDAC/scRNA/TF/")
#加载分析包
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
#可视化相关包，多加载点没毛病
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(limma)
library(SCENIC)
library(AUCell)
library(data.table)

sample2 <- c("P20181121", "P20181121", "P20181121", "P20181121")

Idents(seu_obj) <- "orig.ident"
for (sample_name in sample2) {
  extracted_sample <- subset(seu_obj, idents = sample_name)  # 提取样本
  # 将提取的样本对象保存到与样本名称相同的对象名中
  assign(sample_name, extracted_sample)
}

sce_SCENIC <- open_loom("./P20181121/sce_SCENIC.loom")

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)#将上一步矩阵文件转化为list
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

cellinfo <- P20181121@meta.data[,c('celltype','Label',"nFeature_RNA","nCount_RNA")]#细胞meta信息
colnames(cellinfo)=c('celltype', 'group','nGene' ,'nUMI')

cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

cols_to_keep <- rownames(cellTypes)
sub_regulonAUC <- sub_regulonAUC[, cols_to_keep]

#如果行名和regulonAUC不一致的话，可以用该句改行名
#rownames <- rownames(cellTypes)
#new_rownames <- gsub("P20181121_", "", rownames)
#rownames(cellTypes) <- new_rownames

#rownames <- rownames(cellinfo)
#new_rownames <- gsub("P20181121_", "", rownames)
#rownames(cellinfo) <- new_rownames

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),#从aucellresults获取AUC矩阵
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)#去除含有NA的行

anaAUC <- sub_regulonAUC
anaAUC <- anaAUC[onlyNonDuplicatedExtended(rownames(anaAUC)),]
ana.cellinfo<-P20181121@meta.data#细胞信息

#非批量导出
#Tcell.cellinfo <- subset(ana.cellinfo, ana.cellinfo$celltype=='GZMK+ Tem')
#targets<-data.table(FileName=rownames(Tcell.cellinfo),Target=Tcell.cellinfo$Label)#提取组合分组

#接下来其实就是常规的差异分析了
#lev<-unique(targets$Target)
#f <- factor(targets$Target, levels=lev) 
#design <- model.matrix(~0+f) #样本矩阵
#colnames(design) <- c("Group1","Group0") #design矩阵重命名

#eset=getAUC(anaAUC) #获取所有细胞TF的AUC值
#eset=eset[,targets$FileName]#我们只选取要分析的T细胞
#eset<-t(scale(t(eset))) #scale
#dim(eset)

#cont.wt <- makeContrasts("Group0",levels=design) #GM vs BM
#fit <- lmFit(eset, design)
#fit2 <- contrasts.fit(fit, cont.wt) 
#fit2 <- eBayes(fit2) 
#DEGs=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)#差异分析结果

#DEGs = subset(DEGs, select=c("adj.P.Val","P.Value","logFC"))
#colnames(DEGs)=c("FDR","P.Value","logFC")#选择关键的三列即可
#write.table(data.frame(ID=rownames(DEGs),DEGs), file = "GZMK+Tem_DEGs.txt", sep = "\t", col.names = T, row.names = F, quote = F)

table(P20181121@meta.data$celltype)
subtypes <- c("GZMK+early_Tem", "terminal_Tex", "Proliferating")

for (subtype in subtypes) {
  # 从 ana.cellinfo 中筛选出指定亚型的 T 细胞
  subtype_cellinfo <- subset(ana.cellinfo, ana.cellinfo$celltype == subtype)
  
  # 提取组合分组信息
  targets <- data.table(FileName = rownames(subtype_cellinfo), Target = subtype_cellinfo$Label)
  
  # 常规的差异分析
  lev <- unique(targets$Target)
  f <- factor(targets$Target, levels = lev) 
  design <- model.matrix(~0 + f)
  colnames(design) <- c("Group1", "Group0")
  
  # 获取指定亚型 T 细胞的 AUC 值
  eset <- getAUC(anaAUC)
  eset <- eset[, targets$FileName]
  eset <- t(scale(t(eset)))
  
  # 差异分析
  cont.wt <- makeContrasts("Group0", levels = design)
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, cont.wt) 
  fit2 <- eBayes(fit2) 
  DEGs <- topTable(fit2, adjust = "BH", sort.by = "logFC", n = Inf)
  
  # 选择关键的三列
  DEGs <- subset(DEGs, select = c("adj.P.Val", "P.Value", "logFC"))
  colnames(DEGs) <- c("FDR", "P.Value", "logFC")
  
  # 将结果写入文件
  write.table(data.frame(ID = rownames(DEGs), DEGs), 
              file = paste0(subtype, "_DEGs.txt"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

sce_SCENIC <- open_loom("./PACA.P20181121/sce_SCENIC.loom")
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)#将上一步矩阵文件转化为list
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

cellinfo <- P20181121@meta.data[,c('celltype','Label',"nFeature_RNA","nCount_RNA")]#细胞meta信息
colnames(cellinfo)=c('celltype', 'group','nGene' ,'nUMI')

cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

cols_to_keep <- rownames(cellTypes)
sub_regulonAUC <- sub_regulonAUC[, cols_to_keep]

table(cellTypes)
subtypes <- c("GZMK+early_Tem", "terminal_Tex")

for (subtype in subtypes) {
  # 从 ana.cellinfo 中筛选出指定亚型的 T 细胞
  subtype_cellinfo <- subset(ana.cellinfo, ana.cellinfo$celltype == subtype)
  
  # 提取组合分组信息
  targets <- data.table(FileName = rownames(subtype_cellinfo), Target = subtype_cellinfo$Label)
  
  # 常规的差异分析
  lev <- unique(targets$Target)
  f <- factor(targets$Target, levels = lev) 
  design <- model.matrix(~0 + f)
  colnames(design) <- c("Group1", "Group0")
  
  # 获取指定亚型 T 细胞的 AUC 值
  eset <- getAUC(anaAUC)
  eset <- eset[, targets$FileName]
  eset <- t(scale(t(eset)))
  
  # 差异分析
  cont.wt <- makeContrasts("Group0", levels = design)
  fit <- lmFit(eset, design)
  fit2 <- contrasts.fit(fit, cont.wt) 
  fit2 <- eBayes(fit2) 
  DEGs <- topTable(fit2, adjust = "BH", sort.by = "logFC", n = Inf)
  
  # 选择关键的三列
  DEGs <- subset(DEGs, select = c("adj.P.Val", "P.Value", "logFC"))
  colnames(DEGs) <- c("FDR", "P.Value", "logFC")
  
  # 将结果写入文件
  write.table(data.frame(ID = rownames(DEGs), DEGs), 
              file = paste0(subtype, "_DEGs.txt"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}





























