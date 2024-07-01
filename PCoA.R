rm(list=ls())
# Load package
library(vegan)
library(ggplot2)
library(ggthemes)
# Load data
setwd("D:/M/R/BYD")
otu_PcoA <- read.csv("spe_w.csv",sep=",",header = T,row.names = 1,check.names = F)
otu_PcoA<-t(otu_PcoA)
group1 <- read.csv("groupw.csv",header=TRUE,sep=",")
names(group1)[1] <- c('ID')
#pcoa
# vegdist函数，计算距离；method参数，选择距离类型
distance1 <- vegdist(otu_PcoA, method = 'bray')
# 对加权距离进行PCoA分析
pcoa1 <- cmdscale(distance1, k = (nrow(otu1) - 1), eig = TRUE)
## plot data
# 提取样本点坐标
plot_data1 <- data.frame({pcoa1$point})[1:2]

# 提取列名，便于后面操作。
plot_data1$ID <- rownames(plot_data1)
names(plot_data1)[1:2] <- c('PCoA1', 'PCoA2')

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig1 = pcoa1$eig

#为样本点坐标添加分组信息
plot_data1 <- merge(plot_data1, group1, by = 'ID', all.x = TRUE)
head(plot_data1)

# 计算加权bray-curtis距离
dune_dist1 <- vegdist(otu_PcoA, method="bray", binary=F)
dune_pcoa1 <- cmdscale(dune_dist1, k=(nrow(otu_PcoA) - 1), eig=T)

dune_pcoa_points1 <- as.data.frame(dune_pcoa1$points)
sum_eig1 <- sum(dune_pcoa1$eig)
eig_percent1 <- round(dune_pcoa1$eig/sum_eig1*100,1)
eig_percent1 
colnames(dune_pcoa_points1) <- paste0("PCoA", 1:3)

dune_pcoa_result1 <- cbind(dune_pcoa_points1, group1)

head(dune_pcoa_result1)
library(ggplot2)

ggplot(dune_pcoa_result1, aes(x=PCoA1, y=PCoA2, fill=Group,color=Group,shape = Group)) +
  geom_point(size=4) +
  stat_ellipse(level=0.95,size = 0.6)+
  scale_color_manual(values=c("#2878B5", "#8ECFC9"))+
  labs(x=paste("PCoA 1 (", eig_percent1[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent1[2], "%)", sep=""))+
  geom_hline(yintercept=0,linetype = 3,size = 0.8,color="gray") +
  geom_vline(xintercept=0,linetype = 3,size = 0.8,color="gray")+  theme_bw(base_size=16)+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
#-----------------------------------permonova----------
# 基于bray-curtis距离进行计算
dune.div1 <- adonis2(otu_PcoA ~ Group, data = group1, permutations = 999, method="bray")

dune.div1
library(ggalt)
dune_adonis1 <- paste0("R2= ",round(dune.div1$R2,2), "  p= ", dune.div1$`Pr(>F)`)

p <-ggplot(dune_pcoa_result1, aes(x=PCoA1, y=PCoA2, fill=Group,color=Group,shape = Group)) +
  geom_point(size=4) +
  stat_ellipse(level=0.95,size = 0.6)+
  scale_color_manual(values=c("#2878B5", "#8ECFC9"))+
  labs(x=paste("PCoA 1 (", eig_percent1[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent1[2], "%)", sep=""),
       title=dune_adonis1)+
  geom_hline(yintercept=0,linetype = 3,size = 0.8,color="gray") +
  geom_vline(xintercept=0,linetype = 3,size = 0.8,color="gray")+  theme_bw(base_size=16)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p

#------------------导出文本文件------------------------

write.table (dune_pcoa_result, file ="PCoA-result.csv",
             sep ="", row.names =TRUE, col.names =TRUE, quote =TRUE)
write.table (dune.div, file ="PCoA-pvalue.csv", 
             sep ="", row.names =TRUE, col.names =TRUE, quote =TRUE)
library("vegan")

distance.bray<-vegdist(otu1,method = 'bray')
distance.bray
anosim.result<-anosim(distance.bray,group1$Group,permutations = 999)

summary(anosim.result)

p2<-ggplot(dune_pcoa_result1, aes(x=PCoA1, y=PCoA2, fill=Group,color=Group,shape = Group))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+#隐藏网格线
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+#图中虚线
  guides(color=guide_legend(title=NULL))+#去除图例标题
  labs(x=paste("PCoA 1 (", eig_percent1[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent1[2], "%)", sep=""),
       title=dune_adonis1)+#将x、y轴标题改为贡献度
  stat_ellipse(level=0.95,size = 0.6)+
  scale_color_manual(values=c("#2878B5", "#8ECFC9")) +#点的颜色设置
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10))#修改y轴刻度标签文本


p2
