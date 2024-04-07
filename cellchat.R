setwd("C:/Users/Jacky/Desktop/PDAC/scRNA/细胞通讯")
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
library(CellChat)

sample2 <- c("P15","P19")
NC1 <- readRDS("./P15+P19.rds")

NC1@meta.data <- select(NC1@meta.data, -customclassif)

metatable <- read_excel("./metadata/Neo_T.xlsx")
row_names <- rownames(NC1@meta.data)
NC1@meta.data$cell_id <- row_names  
metadata <- FetchData(NC1, "cell_id")
metadata$cell_id <- rownames(metadata) #??ȡÿ??cell??UMI_ID
metadata <- left_join(x = metadata, y = metatable, by = "cell_id") #?ϲ??ٴ???Ϣ
rownames(metadata) <- metadata$cell_id
NC1 <- AddMetaData(NC1, metadata = metadata)

CellDimPlot(
  srt = NC1, group.by = "fastMNNclusters", reduction = "fastMNNUMAP2D",
  title = "Seurat", theme_use = "theme_blank",label = T
)

CellDimPlot(
  srt = NC1, group.by = "Label", reduction = "fastMNNUMAP2D",
  title = "Seurat", theme_use = "theme_blank",label = T
)

NC1@meta.data$fastMNNclusters

Idents(NC1) <- "orig.ident"
for (sample_name in sample2) {
  extracted_sample <- subset(NC1, idents = sample_name)  # 提取样本
  # 将提取的样本对象保存到与样本名称相同的对象名中
  assign(sample_name, extracted_sample)
}

table(P15@meta.data$Label)

data.input <- P15@assays$RNA$data
meta = P15@meta.data[,c("Label","cell_id")]
unique(meta$Label)

###创建CellChat 对象###
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Label")
###将细胞信息添加到对象的meta slot中, 如果在创建cellchat对象时未添加细胞meta信息
#cellchat <- addMeta(cellchat, meta = meta)
#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
#levels(cellchat@idents) # show factor levels of the cell labels
#groupSize <- as.numeric(table(cellchat@idents)) 

#设置配体受体交互数据库
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis,挑选特定的细胞通讯分子
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB

###预处理用于细胞通信分析的表达数据
cellchat <- subsetData(cellchat) ## subset the expression data of signaling，首先取一下要计算的基因子集，因为CellChatDB.human中包含4万多基因
#future::plan("multisession", workers = 4) #设置使用的电脑核心数，并行计算节省时间，windows系统需将'multisession'改为'multiprocess'
#future::plan("multiprocess", workers = 8)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#推断相互作用网络，预处理完表达值后，cellchat会计算相互作用的可能性和相互作用网络：
# trim参数用来过滤掉只在少量细胞表达的基因，
# 0.1表示如果某基因在一种细胞类型中表达比例少于10%，
# 则该基因在这种细胞中平均表达量归0，也就是该配体/受体不在这种细胞中表达
# population.size设置为TRUE表示需要考虑每种细胞类型中细胞数量的多少，对于分选的细胞需要设置为FALSE
cellchat <- computeCommunProb(cellchat, trim = 0.1, population.size = FALSE) 
# min.cells用来过滤掉细胞数目较少的细胞类型
cellchat <- filterCommunication(cellchat, min.cells = 30)
# 在信号通路水平推断相互作用，到这里细胞通讯的推断已经完成，所有数据都存储在了cellchat对象中，接下来就是可视化了。
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@LR$LRsig %>% View()  # 查看受体配体库，用于pathway_name和interaction_name参数

###总的相互作用网络:
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
)

#每种细胞分别作为source的相互作用网络:
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#所有显示重要通信的信号通路均可通过cellchat@netP$pathways获取。
cellchat@netP$pathways
pathways.show <- c("CD99")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CXCL")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("ADGRE5")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("NECTIN")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CD46")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CD80")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("TIGIT")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CD86")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("IFN-II")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("VCAM")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("LCK")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("PVR")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CD80")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CD137")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("SN")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("CEACAM")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("PD-L1")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("PDL2")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("MHC-I")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

pathways.show <- c("ICAM")
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
pathways.show <- c("CD80")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("TS_CD8", 6), rep("non_TS_CD8", 6), rep("Proliferating_epithelial", 6)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

###计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
pathways.show <- c("CD80", "CXCL", "ADGRE5", "NECTIN", "CD46", "TIGIT", "CD86", "IFN-II", "PD-L1", "PDL2", "CEACAM")
netAnalysis_contribution(cellchat, signaling = pathways.show)

#我们还可以可视化由单个配体受体对调节的细胞-细胞通信。我们提供一个函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因。
pairLR.CLEC <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CLEC[3,] # show one ligand-receptor pair

# Hierarchy plot
pathways.show <- c("TIGIT")
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#自动保存所有推断网络的模块以进行快速探索
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

###可视化由多个配体受体或信号通路调节的细胞通信
#气泡图
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:11), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c("TS_CD8","non_TS_CD8"), targets.use = c(1:32), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c("TS_CD8"), targets.use = c(1:17), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c("non_TS_CD8"), targets.use = c(1:17), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 10, targets.use = c(1:11), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use =  c("TS_CD8","non_TS_CD8"), targets.use = c(1:32), signaling = c("TIGIT","CXCL","CD80"), remove.isolate = FALSE)
netVisual_chord_gene(cellchat, sources.use = c("TS_CD8","non_TS_CD8"), targets.use = c(1:32), lab.cex = 0.5,legend.pos.y = 30)
# show all the interactions received by Neo
netVisual_chord_gene(cellchat, sources.use = c(1:11), targets.use = 8, legend.pos.x = 30)
netVisual_chord_gene(cellchat, sources.use = c("TS_CD8","non_TS_CD8"), targets.use = c(1:32), legend.pos.x = 30)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1:11), targets.use = 8, signaling = c("CCL","CLEC"),legend.pos.x = 8)
netVisual_chord_gene(cellchat, sources.use = c("TS_CD8","non_TS_CD8"), targets.use = 8, signaling = c("CCL","CLEC"),legend.pos.x = 8)

###使用小提琴/点图绘制信号基因表达分布
plotGeneExpression(cellchat, signaling = "CD80")
plotGeneExpression(cellchat, signaling = "MHC-I", enriched.only = FALSE)

plotGeneExpression(cellchat, signaling = "MHC-I", enriched.only = TRUE)
plotGeneExpression(cellchat, signaling = c("CD80","MHC-I","TIGIT","CXCL","ICAM"), enriched.only = FALSE)
library(ggplot2)
your_gene <- "ACKR1"
subtypes <- c("Endothelial", "TS_CD8")

table(P15@meta.data$Label)

###第四部分：细胞通信网络系统分析
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("ICAM")

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 4, font.size = 10)

#在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("TIGIT", "CD80","CXCL"))
gg1 + gg2
#识别对某些细胞组的传出或传入信号贡献最大的信号
#我们还可以回答以下问题：哪些信号对某些细胞组的传出或传入信号贡献最大。
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 20, height = 20, font.size = 6)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 20, height = 20, font.size = 6)
ht1 + ht2
ht3 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CLEC","CEACAM","CXCL","ADGRE5","TIGIT","CD99","LCK","CD96"), pattern = "outgoing") #特定信号通路的展示
ht4 <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I","NECTIN","CD86","VCAM","LCK","CD99","PRV","CD80","CD137","SN","PD-L1"), pattern = "incoming") #特定信号通路的展示
ht3 + ht4

###确定全局通信模式，探索多个细胞类型和信号通路如何协调在一起
#识别和可视化分泌细胞的传出通信模式
library(NMF)
library(ggalluvial)
#运行selectK推断模式的数量
selectK(cellchat, pattern = "outgoing") #当传出模式数为 3 时，Cophenetic 和Silhouette值都开始突然下降
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")
#识别和可视化目标细胞的传入通信模式
selectK(cellchat, pattern = "incoming")
nPatterns = 5 #当传入模式的数量为 5 时，Cophenetic 值开始下降。
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchat_P15.rds")















































