setwd("C:/Users/Jacky/Desktop/PDAC/scTCR")
library(scRepertoire)
library(gridExtra)
library(ggplot2)
library(plotly)

S1 <- read.csv("./matrix/all_contig_annotations_S01_Neo.csv")
S2 <- read.csv("./matrix/all_contig_annotations_S01_nonNeo.csv")
S3 <- read.csv("./matrix/all_contig_annotations_S02_Neo.csv")
S4 <- read.csv("./matrix/all_contig_annotations_S02_nonNeo.csv")
S5 <- read.csv("./matrix/all_contig_annotations_S03_Neo.csv")
S6 <- read.csv("./matrix/all_contig_annotations_S03_nonNeo.csv")
S7 <- read.csv("./matrix/all_contig_annotations_S04.csv")
S8 <- read.csv("./matrix/all_contig_annotations_S05.csv")
S9 <- read.csv("./matrix/all_contig_annotations_S06.csv")
S10 <- read.csv("./matrix/all_contig_annotations_S07.csv")
S11 <- read.csv("./matrix/all_contig_annotations_S08_Neo.csv")
S12 <- read.csv("./matrix/all_contig_annotations_S08_nonNeo.csv")
S13 <- read.csv("./matrix/all_contig_annotations_S09_Neo.csv")
S14 <- read.csv("./matrix/all_contig_annotations_S09_nonNeo.csv")
S15 <- read.csv("./matrix/all_contig_annotations_S10_Neo.csv")
S16 <- read.csv("./matrix/all_contig_annotations_S10_nonNeo.csv")
S17 <- read.csv("./matrix/all_contig_annotations_S11_Neo.csv")
S18 <- read.csv("./matrix/all_contig_annotations_S11_nonNeo.csv")

contig_list <- list(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16)

combined <- combineTCR(contig_list, 
                       samples = c("WY_01_Neo", "WY_01_Non","WY_02_Neo","WY_02_Non","WY_03_Neo", "WY_03_Non", "BL5", "BL10", "BL11", "BL12", "TIPC249_Neo", "TIPC249_Non", "TIPC262_Neo", "TIPC262_Non","TIPC282_Neo","TIPC282_Non"), 
                       cells ="T-AB", removeNA = T, removeMulti = F)

###可以使用addVariable函数添加其他信息，如tissue ，age ，分期等，这个很实用，后续可以使用这些分组进行 clone 的各类比较
combined <- addVariable(combined, name = "Tissue",  
                       variables = c("PDAC", "PDAC", "PDSC","PDSC", "PDAC", "PDAC", "CP","CP", "CP","CP","PDAC", "PDAC","PDAC", "PDAC", "PDAC","PDAC"))

combined <- addVariable(combined, name = "Specific",  
                        variables = c("Neo", "Non", "Neo","Non", "Neo","Non", "Non","Non", "Non","Non","Neo","Non","Neo","Non", "Neo","Non"))

###2，VCJ分析以及可视化
#使用quantContig 探索每个样本的unique clone信息
#百分比
p1 <- quantContig(combined, cloneCall="gene+nt", scale = T)+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))#X坐标加点斜体，出图好看
#绝对数量
p2 <- quantContig(combined, cloneCall="gene+nt", scale = F)+ 
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))

注：这里cloneCall函数有四个参数
(1) “gene” ：使用包含 TCR/Ig 的 VDJC 基因
(2) “nt”：使用 CDR3 区域的核苷酸序列
(3) “aa” ：使用 CDR3 区域的氨基酸序列
(4) “gene+nt” 使用包含 TCR/Ig + CDR3 区域的核苷酸序列的 VDJC 基因
###柱形图具体的数值可以添加 exportTable = T函数导出
quantContig_output <- quantContig(combined, cloneCall="gene+nt",    
                                  scale = T, exportTable = T)
write.table(data.frame(ID=rownames(quantContig_output),quantContig_output), file = "unique_clonotype.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#2.2 Clonotype Abundance克隆型丰度
abundanceContig(combined, cloneCall = "gene", scale = F)
#2.3 Length of Clonotypes克隆型长度,使用lengtheContig函数查看 CDR3 序列的长度分布。
p11 <- lengthContig(combined, cloneCall="aa", chain = "both", group.by = "Tissue")#左图
p12 <- lengthContig(combined, cloneCall="aa", chain = "both", group.by = "Specific")#左图
p12 <- p12 + scale_fill_manual(values=c("#e64b35", "#4dbbd5"))

p22 <- lengthContig(combined, cloneCall="aa", chain = "TRB") #右图

注意：chain有两种参数
1. “both” for combined chain visualization
2. “TRA”, “TRB”, “TRD”, “TRG”, “IGH” or “IGL” to select single chain

#Compare Clonotypes比较克隆型,使用compareClonotypes函数查看样本之间的克隆型的比例和动态变化。
WY_02 <- compareClonotypes(combined, numbers = 10, 
                  samples = c("WY_02_Neo","WY_02_Non"), #samples可以选择任意样本，但是如果top number clonotype sequences之间没有共享的话，则没有连线。 
                  #samples = c("t2_Normal", "t2_PBMC","t3_PBMC","t3_Normal","t3_Center"),  
                  cloneCall="aa",  #cloneCall可以选择“gene+nt”
                  graph = "alluvial")

HPV51 <- compareClonotypes(combined, numbers = 10, 
                          samples = c("HPV51_KSA_T","HPV51_Non_T","HPV51_QVD_T","HPV51_Non_m","HPV51_QVD_m"), #samples可以选择任意样本，但是如果top number clonotype sequences之间没有共享的话，则没有连线。 
                          #samples = c("t2_Normal", "t2_PBMC","t3_PBMC","t3_Normal","t3_Center"),  
                          cloneCall="aa",  #cloneCall可以选择“gene+nt”
                          graph = "alluvial")

HPV34 <- compareClonotypes(combined, numbers = 10, 
                           samples = c("HPV34_Non_T","HPV34_QVD_T","HPV34_QVD_m","HPV34_KSA_m","HPV34_Non_m"), #samples可以选择任意样本，但是如果top number clonotype sequences之间没有共享的话，则没有连线。 
                           #samples = c("t2_Normal", "t2_PBMC","t3_PBMC","t3_Normal","t3_Center"),  
                           cloneCall="aa",  #cloneCall可以选择“gene+nt”
                           graph = "alluvial")

HPV15 <- compareClonotypes(combined, numbers = 10, 
                           samples = c("HPV15_QVD_T","HPV15_KSA_T"), #samples可以选择任意样本，但是如果top number clonotype sequences之间没有共享的话，则没有连线。 
                           #samples = c("t2_Normal", "t2_PBMC","t3_PBMC","t3_Normal","t3_Center"),  
                           cloneCall="aa",  #cloneCall可以选择“gene+nt”
                           graph = "alluvial")

saveRDS(combined,"PDAC_scTCR.rds")
combined <- readRDS("PDAC_scTCR.rds")

###Clonal Space Homeostasis克隆空间稳态,通过clonalHomeostasis函数查看， 各特定比例的克隆型（cloneTypes：Rare，，，Hyperexpanded）所占据该样本的相对比例
library(cols4all)
mycol <- c4a('vivid',12)

clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04,
                                 Small = 0.001,
                                 Medium = 0.01,
                                 Large = 0.1,
                                 Hyperexpanded = 1))

clonalHomeostasis(combined, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04,
                                 Small = 0.001,
                                 Medium = 0.01,
                                 Large = 0.1,
                                 Hyperexpanded = 1)) +
  scale_fill_manual(name = "Clonotype Group", 
                    values = mycol)


#Clonal Proportion克隆比例,通过clonalProportion函数查看克隆型的比例，按克隆型的出现频率将其进行排名，1:10表示每个样本中的前10个克隆型。
clonalProportion(combined, cloneCall = "nt")

#Overlap Analysis 重叠分析, 使用clonalOverlap函数分析样本之间的相似性，使用 clonesizeDistribution函数对样本进行聚类。
clonalOverlap(combined, cloneCall = "gene+nt",   
              method = "overlap")

clonesizeDistribution(combined, 
                      cloneCall = "gene+nt", 
                      method="ward.D2")

###Diversity Analysis多样性分析
clonalDiversity(combined, 
                cloneCall = "gene+nt",  
                n.boots = 100)

clonalDiversity(combined, 
                cloneCall = "gene+nt",  
                group = "Specific",
                n.boots = 100)

combined <- readRDS("scTCR.rds")

WY_01 <- combined$WY_01
WY_02 <- combined$WY_02
WY_03 <- combined$WY_03
BL5 <- combined$BL5
BL10 <- combined$BL10
BL11 <- combined$BL11
BL12 <- combined$BL12
TIPC249 <- combined$TIPC249
TIPC262 <- combined$TIPC262
TIPC282 <- combined$TIPC282
TIPC301_413_416 <- combined$TIPC301_413_416

var_names <- c("WY_01", "WY_02", "WY_03", "BL5", "BL10", "BL11", "BL12", "TIPC249", "TIPC262", "TIPC282", "TIPC301_413_416")  

for (var_name in var_names) {  
  # 使用as.symbol将字符串转换为符号，然后用eval求值  
  data <- eval(as.symbol(var_name))  
  
  # 生成文件名，假设文件名和变量名一致  
  file_name <- paste0(var_name, "_T.txt")  
  
  # 保存文件，这里使用了write.table，你也可以选择write.csv等函数  
  write.table(data.frame(ID = rownames(data), data),   
              file = file_name,   
              sep = "\t",   
              col.names = TRUE,   
              row.names = FALSE,   
              quote = FALSE)  
}


write.table(data.frame(ID=rownames(WY_01),WY_01), file = "WY_01_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV7_Non_m),HPV7_Non_m), file = "HPV7_Non_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV7_Non_T),HPV7_Non_T), file = "HPV7_Non_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV7_QVD_m),HPV7_QVD_m), file = "HPV7_QVD_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)

write.table(data.frame(ID=rownames(HPV34_KSA_m),HPV34_KSA_m), file = "HPV34_KSA_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV34_Non_m),HPV34_Non_m), file = "HPV34_Non_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV34_Non_T),HPV34_Non_T), file = "HPV34_Non_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV34_QVD_m),HPV34_QVD_m), file = "HPV34_QVD_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV34_QVD_T),HPV34_QVD_T), file = "HPV34_QVD_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)


write.table(data.frame(ID=rownames(HPV51_KSA_T),HPV51_KSA_T), file = "HPV51_KSA_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV51_Non_m),HPV51_Non_m), file = "HPV51_Non_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV51_Non_T),HPV51_Non_T), file = "HPV51_Non_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV51_QVD_m),HPV51_QVD_m), file = "HPV51_QVD_m.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(ID=rownames(HPV51_QVD_T),HPV51_QVD_T), file = "HPV51_QVD_T.txt", sep = "\t", col.names = T, row.names = F, quote = F)


























































