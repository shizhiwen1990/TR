setwd("C:/Users/Jacky/Desktop/test")
#Let us filter combine the data of 32
rm(list=ls())
library(latticeExtra)
#patientInfo<-read.csv("./annotation/180418LGLExpList.csv", row.names = 1)
patientInfo<-read.csv("./metadata/sample_information_1.csv", row.names = 1)
patientInfoAll<-paste(patientInfo[,1],patientInfo[,2],patientInfo[,3], sep="_")
names(patientInfoAll)<-rownames(patientInfo)
setwd("C:/Users/Jacky/Desktop/PDAC/scTCR/matrix")

for(sample in paste("S0", 1:16, sep="")){
  dataFiltered1A1B<-read.csv(paste("all_contig_annotations_", sample, ".csv", sep=""), header=T)
  dataPlot<-as.matrix(table(dataFiltered1A1B[, c("v_gene", "j_gene")]))
  dataPlot3<-as.data.frame(matrix(1, nrow=dim(dataPlot)[1]*dim(dataPlot)[2], ncol=3))
  colnames(dataPlot3)<-c("x", "y", "z")
  for(ii in 1:dim(dataPlot)[1]){
    for(jj in 1:dim(dataPlot)[2]){
      dataPlot3[(ii-1)*dim(dataPlot)[2]+jj, 1]<-rownames(dataPlot)[ii]
      dataPlot3[(ii-1)*dim(dataPlot)[2]+jj, 2]<-as.integer(dataPlot[ii,jj])
      dataPlot3[(ii-1)*dim(dataPlot)[2]+jj, 3]<-colnames(dataPlot)[jj]
    }
  }
  dataPlot3$x<-as.factor(dataPlot3$x)
  dataPlot3$y<-as.integer(dataPlot3$y)
  dataPlot3$z<-as.factor(dataPlot3$z)
  #png(paste("result/", sample,  patientInfoAll[sample], ".png", sep=""), width=3000, height=3000, res=50)
  #print(cloud(y~x+z, dataPlot3, panel.3d.cloud=panel.3dbars, col.facet='red', 
  #    xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
  #    par.settings = list(axis.line = list(col = "transparent")), aspect=c(0.3,0.3)))
  #dev.off()
  dataPlot3<-dataPlot3[order(dataPlot3$y, decreasing = TRUE),]
  write.csv(dataPlot3, file=paste("result/", sample,  patientInfoAll[sample], "_barPlot.csv", sep=""))
}


k <- read.table(text = 'x y z

TRGV4	1	TRBJ2-5',header=TRUE)

k$x <- factor(k$x)
k$z <- factor(k$z)

p1 <- cloud(y~x+z, k, panel.3d.cloud=panel.3dbars, col.facet='#99C945', 
            xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
            par.settings = list(axis.line = list(col = "transparent")), aspect=c(0.3,0.3))
p1

p1 <- cloud(y ~ x + z, k, panel.3d.cloud = panel.3dbars, col.facet = '#e64b35', 
            xbase = 0.8, ybase = 0.8, border = NA, 
            scales = list(arrows = FALSE, col = 1), 
            par.settings = list(axis.line = list(col = "transparent")), 
            aspect = c(0.3, 0.3))

p1 <- cloud(y ~ x + z, k, panel.3d.cloud = panel.3dbars, col.facet = '#e64b35', 
            xbase = 1, ybase = 1, border = NA, 
            scales = list(arrows = FALSE, col = 1), 
            par.settings = list(axis.line = list(col = "transparent"),
                                font.axis = list(cex = 0.3),   # 调整轴标签字体大小
                                font.main = list(cex = 1.5)),  # 调整主标题字体大小
            aspect = c(0.3, 0.3))




































