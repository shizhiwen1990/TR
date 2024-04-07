

#=======================================================================================
#                                   pySCENIC
#=======================================================================================
官方教程：
[pySCENIC — pySCENIC latest documentation](https://pyscenic.readthedocs.io/en/latest/)

#----------------------以下内容在Linux环境中操作----------------------------------------

#安装下载及环境设置
# 安装一个conda，为什么安装他可以理解为Rstuido之于R,后期在环境设置、软件安装上很方便。
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash安装，按照指引，都选yes，这样一些依赖的python包都安装了。
bash Miniconda3-latest-Linux-x86_64.sh 
#  激活conda环境
source ~/.bashrc

#设置镜像
conda config --add channels r 
conda config --add channels conda-forge 
conda config --add channels bioconda
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/
conda config --set show_channel_urls yes

#安装pyscenic，并创建分析环境
conda create -n pyscenic python=3.9#创建一个pyscenic 的python环境，pyscenic要求python版本3.6及以上，目前python出到3.9了，我用3.9
conda activate pyscenic #激活pyscenic 环境
#安装依赖包
conda install -y numpy
conda install -y -c anaconda cytoolz
conda install -y scanpy
#安装pyscenic
pip install pyscenic
#文件下载
#TF注释
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
#转录组因子列表
下载地址
https://github.com/aertslab/pySCENIC/tree/master/resources
#reference
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
#文件准备
#pyscenic的输入文件是行为基因名，列为细胞ID的矩阵，所以在seurat对象中导出矩阵的时候需要转置一下，可以用标准化矩阵，也可以用counts矩阵，影响不大！
#表达矩阵、meta----R中进行
write.csv(t(as.matrix(human_data@assays$RNA@counts)),file = "sce_exp.csv")
# cellInfo <- human_data@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
# colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
# head(cellInfo)
# write.csv(cellInfo, file = "cellInfo.csv")

#转化为loom文件，Linux下的python脚本
#编辑脚本
vim trans.py
#输入以下内容
import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("sce_exp.csv");#R中导出的表达矩阵
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sce.loom",x.X.transpose(),row_attrs,col_attrs)

#保存并退出
#运行trans.py
python trans.py

ls
#这样在文件夹中会出现sce.loom文件，就是接下来输入pyscenic的文件。


#=======================================================================================
#                                   pySCENIC数据分析
#=======================================================================================

#分析第一步：GRN---运行完得到sce.adj.csv文件
##这一步的目的
#推断转录因子与提供的表达矩阵基因的共表达模块，基于grnboost2，R中时GENIE3
pyscenic grn --num_workers 10 \
--sparse \
--method grnboost2 \
--output sce.adj.csv \
sce.loom \
hs_hgnc_tfs.txt

#分析第二步：RcisTarget---运行完得到sce.regulons.csv文件
#这一步的目的
#进行TF-motif富集分析，识别直接靶标
#得到转录因子(TF)与其对应的直接作用的靶点,称为regulon(每一个regulon是1个TF和其调控的靶基因)
pyscenic ctx --num_workers 10 \
--output sce.regulons.csv \
--expression_mtx_fname sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
sce.adj.csv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather


#分析第三步：AUCell---运行完得到sce_SCENIC.loom文件，即分析结果
#这一步的目的
#使用AUCell对每个细胞的每个regulon活性进行评分。
pyscenic aucell --num_workers 3 \
--output sce_SCENIC.loom \
sce.loom \
sce.regulons.csv




#############################################################################
#                                在R中可视化 
#############################################################################

setwd("/home/shpc_100828/Pyscenic/")
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

##读取pyscenic第三步分析的文件sce_SCENIC.sce_SCENIC
sce_SCENIC <- open_loom("sce_SCENIC.loom")
# exprMat <- get_dgem(sce_SCENIC)#从sce_SCENIC文件提取表达矩阵
# exprMat_log <- log2(exprMat+1) # log处理
#这里的表达矩阵其实就是我们在pyscenic分析第一步的输入矩阵，可见这些文件都是在一起的

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
#提取第二步分析的regulons,column.attr.name填Regulons，具体按照实际情况提示选择
#也就是我们所输入的基因和找到的转录因子组成的表达文件
regulons <- regulonsToGeneLists(regulons_incidMat)#将上一步矩阵文件转化为list
class(regulons)


#提取pyscenic第三步分析中AUC结果
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

#以上就是一些主要文件了、够后续分析和可视化
#####################################################################################
##==============================加载seurat对象、RSS分析=======================================
#在可视化之前，我们再做一个分析，计算RSS值，计算regulon特异性评分
human_data <- readRDS("~/Pyscenic/human_data.rds")
cellinfo <- human_data@meta.data[,c('celltype','group',"nFeature_RNA","nCount_RNA")]#细胞meta信息
colnames(cellinfo)=c('celltype', 'group','nGene' ,'nUMI')
######计算细胞特异性TF
#在实际数据分析应用中，我认为比较靠谱的应用在于，细胞分了亚群，例如macrophage，有不同的特征
#我们可以查看不同亚群特异性的TF，有助于了解亚群的功能！！！！
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),#从aucellresults获取AUC矩阵
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)#去除含有NA的行
#可视化细胞特异性TF
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 3,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')

rssPlot

#提取数据，可以自己可视化dotplot，或者热图
rss_data <- rssPlot$plot$data
devtools::install_github("XiaoLuo-boy/ggheatmap")
library(ggheatmap)
library(reshape2)
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
col_ann <- data.frame(group= c(rep("Neutrophil",1),
                               rep("Macrophage",1),
                               rep("mDC",1),
                               rep("T cell",1),
                               rep("Mast",1)))#列注释
rownames(col_ann) <- colnames(rss_data)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC")
names(groupcol) <- c("Neutrophil","Macrophage","mDC", "T cell","Mast")
col <- list(group=groupcol)

text_columns <- sample(colnames(rss_data),0)#不显示列名

p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p

######既然可以分析细胞中的特异性TF，那么也可以分析样本中特异性的TF
#我们可以提取某一个细胞的表达矩阵，做pyscenic，后期再进行RSS分析的时候
#可以使用样本，就可以看出哪些TF是这个样本细胞中特异性的TF，这样生物学意义也就有了
##############################################################################

##==============================TF_AUC与seurat结合===========================
#普通展示
next_regulonAUC <- regulonAUC[,match(colnames(human_data),colnames(regulonAUC))]
dim(next_regulonAUC)#将AUC结果于seurat对象结合

regulon_AUC <- regulonAUC@NAMES
human_data@meta.data = cbind(human_data@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))

#自己选定感兴趣的或者比较重要的转录因子，这里我是随机的
TF_plot <- c("ZNF561(+)","FOXP3(+)","YY1(+)","HOXB2(+)",
             "TBX21(+)","TCF12(+)","STAT2(+)","SOX21(+)",
             "RBBP5(+)","NR2F6(+)","NELFE(+)","MAFG(+)")

DotPlot(human_data, features = TF_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

####

DotPlot(human_data, features = TF_plot, group.by = 'group')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  theme(legend.direction = "horizontal", 
        legend.position = "bottom")+
  labs(x=NULL,y=NULL)

###
FeaturePlot(human_data, features ="FOXP3(+)")


#############################################################################

##==============================TF平均表达活性===========================

cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])
#去除extend的TF
# sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
# dim(sub_regulonAUC)

#计算平均表达
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

#scale处理\类似于热图数据的标准化
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 


regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

#热图
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byGroup_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6),
                                   show_row_names = F)) 
hm


################################################------------------------
#------rank可视化rss-----------------------------------------------------

B_rss <- as.data.frame(rss)#rss特异性TF结果
#需要作图的细胞类型
celltype <- c("Neutrophil","Macrophage","mDC","T cell","Mast")
rssRanklist <- list()

for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), color="black", size=3, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', 
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[i]] <- p
}


library(cowplot)

plot_grid(rssRanklist[[1]],rssRanklist[[2]],rssRanklist[[3]],
          rssRanklist[[4]],rssRanklist[[5]],ncol=3)



#=======================================================================================
#                              pySCENIC的差异分析及其他思路
#=======================================================================================
#加载R包
library(limma)
library(SCENIC)
library(AUCell)
library(data.table)


anaAUC <- regulonAUC
anaAUC <- anaAUC[onlyNonDuplicatedExtended(rownames(anaAUC)),]
ana.cellinfo<-human_data@meta.data#细胞信息

#选取需要分析的细胞，这里我们以T cell为例子
Tcell.cellinfo <- subset(ana.cellinfo, ana.cellinfo$celltype=='T cell')
targets<-data.table(FileName=rownames(Tcell.cellinfo),Target=Tcell.cellinfo$group)#提取组合分组

#接下来其实就是常规的差异分析了
lev<-unique(targets$Target)
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #样本矩阵
colnames(design) <- lev #design矩阵重命名

eset=getAUC(anaAUC) #获取所有细胞TF的AUC值
eset=eset[,targets$FileName]#我们只选取要分析的T细胞
eset<-t(scale(t(eset))) #scale

#比较分析limma
cont.wt <- makeContrasts("GM-BM",levels=design) #GM vs BM
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)#差异分析结果

tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")#选择关键的三列即可
write.table(tT,file="Tcell_TF_GM_BM.csv",sep="\t",quote=F)#保存结果

#接下来就可以按照不同的阈值筛选显著的TF了，其实这样的差异分析对于结果的解释和
#应用上更有说服力，这样我们可以直观的看到某些转录因子在不同细胞状态下的差异
#可视化的话可以用热图、火山图等展示。

logFoldChange<-0.5
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]#显著筛选
diffSig=na.omit(diffSig)
rownames(diffSig)<-gsub("\\(","",rownames(diffSig))
rownames(diffSig)<-gsub("\\)","",rownames(diffSig))
rownames(diffSig)<-gsub("\\+","",rownames(diffSig))
write.csv(diffSig, file = "diffSig_TF.csv")



#我们用热图展示以下差异结果
#先计算每个样本细胞的平均表达值
regulonActivity_byCellType <- sapply(split(rownames(Tcell.cellinfo), Tcell.cellinfo$orig.ident),
                                     function(cells) rowMeans(getAUC(anaAUC)[,cells]))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType
#这里有一个基础知识，很多文章展示TF不显示+号，用gsub函数就可以去除
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\(","",rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\)","",rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\+","",rownames(regulonActivity_byCellType_Scaled))

#选择差异显著的TF进行可视化
selTF <- c("YY1","JUND","BCLAF1","SAP30","NR3C1","TAF7","ZNF148",
           "UBTF","FOS","ELF1","MXI1", "ZFX","ZMIZ1","TWIST1",
           "CEBPD","HOXD9","ETV4","TFDP1","BRF1","ZNF423","ZNF420",
           "HDAC2","GLIS1")

sel_regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[selTF,]

annotation_col = data.frame(group = c("BM","BM","BM","BM", "GM", "GM"))
rownames(annotation_col)<-colnames(regulonActivity_byCellType_Scaled)
ann_colors =  list(group = c("BM" = "#619cff", "GM" = "#f8766d"))

pheatmap::pheatmap(sel_regulonActivity_byCellType_Scaled,
                   color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100),
                   breaks=seq(-2, 2, length.out = 100),
                   border_color=NA,
                   cluster_rows = F,
                   cluster_cols=F,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   scale="row",
                   main = "T cell TF")

#=======================================================================================
#                                   pySCENIC 转录因子与靶基因
#=======================================================================================
#我们还可以可视化感兴趣的TF调控的靶基因
#TF与靶基因的关系在pyscenic分析得到的第二个文件

sce_regulons <- read.csv("sce.regulons.csv")
sce_regulons <- sce_regulons[-2, ]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ]
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

#举例子我这里关注FOXP3和TGIF1这两个TF

FOXP3 <- subset(sce_regulons, TF=='FOXP3')
FOXP3 <- FOXP3[which(FOXP3$AUC>0.1),]
FOXP3 <- FOXP3[, c("TF","TargetGenes")]
FOXP3$TargetGenes <-gsub("\\[","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\]","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\(","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\)","",FOXP3$TargetGenes)
FOXP3$TargetGenes <-gsub("\\'","",FOXP3$TargetGenes)
library(stringr)
split_FOXP3<-str_split(FOXP3$TargetGenes,",")
FOXP31 <- as.data.frame(split_FOXP3[[1]])
FOXP32<- as.data.frame(split_FOXP3[[2]])
FOXP33<-as.data.frame(split_FOXP3[[3]])
FOXP32<-as.data.frame(split_FOXP3[[4]])
FOXP35<- as.data.frame(split_FOXP3[[5]])

names(FOXP31) <- 'TF'
names(FOXP32) <- 'TF'
names(FOXP33) <- 'TF'
names(FOXP34) <- 'TF'
names(FOXP35) <- 'TF'

FOXP3 <- rbind(FOXP31,FOXP32,FOXP33,FOXP34,FOXP35)

FOXP3_target <- FOXP3[seq(1,nrow(FOXP3),2), ]
FOXP3_score <- FOXP3[seq(0,nrow(FOXP3),2), ]

FOXP3_gene <- data.frame(FOXP3_target,FOXP3_score)
FOXP3_gene <- FOXP3_gene[!duplicated(FOXP3_gene$FOXP3_target), ]
FOXP3_gene$gene <- 'FOXP3'
colnames(FOXP3_gene) <- c("target","score",'tf')

#同理得到TGIF1及其靶基因，此处省略1万字
#two thousand years later
#TGIF1_gene

TF_target <- rbind(FOXP3_gene,TGIF1_gene)
TF_target$score <- as.numeric(TF_target$score)
###接下来就是网络图了

#节点数据
paths <- c("FOXP3", "TGIF1")#列重命名
nodelist <- list()
for (i in 1:length(paths)){
  node <- subset(TF_target, tf == paths[i])#提取数据
  
  nodes <- data.frame(name = unique(union(node$tf, node$target)))#整理为datafram
  nodes$value <- c(sum(node$score)/10, node$score)#加上values
  
  nodelist[[i]] <- nodes
}  #提取每个大节点数据


nodes <- rbind(nodelist[[1]],nodelist[[2]])#将三个节点文件合并
nodes$cluster <- c(rep("FOXP3",1),rep("FOXP3_gene",20),
                   rep("TGIF1",1),rep("TGIF1_gene",6))#分组，为了后续颜色设置

edges <- TF_target[c("tf","target","score")]#边缘文件
edges$class <- edges$tf

library(ggraph)
library(tidygraph)
layout_cir <- tbl_graph(nodes = nodes, edges = edges)#构建ggraph作图文件
#作图
ggraph(layout_cir,layout='linear',circular = TRUE) +#选择circle
  geom_node_point(aes(size=value,colour = cluster))+#节点，大小用我们赋的值表示，颜色用分组
  geom_node_text(aes(x = 1.03 * x,
                     y = 1.03 * y,
                     label=name,
                     color=cluster,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                 hjust='outward') +#文字设置。x，y是为了调整位置。angle是为了调整角度，以后其他所有网络图angle都用此公式，文字会向外发散排列
  geom_edge_arc(aes(colour=class))+#连线为曲线
  theme_void()+#theme主题
  theme(legend.position = "none")+
  scale_colour_manual(values =c('#407972',
                                '#961E28',
                                '#D46724',
                                '#0f8096'))+#节点颜色
  scale_edge_colour_manual(values = c('#961E28',
                                      '#D46724',
                                      '#0f8096'))+#连线颜色
  scale_size_continuous(range = c(2,8))+#点的大小范围设置
  coord_cartesian(xlim=c(-1.5,1.5),ylim = c(-1.5,1.5))#设置坐标位置，防止图溢出作图边界显示不全




