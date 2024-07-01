setwd("D:/M/R/BYD")
library(vegan)
library(ggplot2)
#OTU表格
df <- read.csv("spe_w.csv",sep=",",header = T,row.names = 1,check.names = F)
head(df)
#对OTU表格进行hellinger转化

#环境因子数据
env <- read.csv("wenv.csv",sep=",",header = T,row.names = 1,check.names = F)
head(env) 
#对环境因子数据进行log10对数转换
log_data<-log10(env+1)
print(decorana(t(df))) 
RDA <- rda(t(df),env,scale = T)
#提取数据
df_rda <- data.frame(RDA$CCA$u[,1:2],rownames(env))
colnames(df_rda)=c("RDA1","RDA2","samples")
# 提取物种得分
df_rda_score <- data.frame(RDA$CCA$v[,1:2])
#计算轴标签数据（=轴特征值/sum(所有轴的特征值)）
RDA1 =round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100,2)
RDA2 =round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100,2) 
#读入分组文件
group <- read.csv("wgroup.csv", sep=',', header=T)
#修改列名
colnames(group) <- c("samples","group")
#将绘图数据和分组合并
df_rda <- merge(df_rda,group,by="samples")
color=c("#1597A5","#FFC24B")#颜色变量
p1<-ggplot(data=df_rda,aes(x=RDA1,y=RDA2,
                           color=group))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=3,shape=16)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+#图中虚线
  geom_text(aes(label=samples, y=RDA2+0.03,x=RDA1+0.03,  vjust=0),size=3)+#添加数据点的标签
  # guides(color=guide_legend(title=NULL))+#去除图例标题
  labs(x=paste0("RDA1 (",RDA1,"%)"),
       y=paste0("RDA2 (",RDA2,"%)"))+#将x、y轴标题改为贡献度
  stat_ellipse(data=df_rda,
               level=0.95,
               linetype = 2,size=0.8,
               show.legend = T)+
  scale_color_manual(values = color) +#点的颜色设置
  scale_fill_manual(values = c("#1597A5","#FFC24B"))+
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
p1
#提取环境因子得分
df_rda_env <- RDA$CCA$biplot[,1:2]
df_rda_env <- as.data.frame(df_rda_env)
head(df_rda_env)
# 添加环境因子数据
p1+
  geom_segment(data=df_rda_env,aes(x=0,y=0,xend=df_rda_env[,1],yend=df_rda_env[,2]),
               color="black",size=0.5,
               arrow=arrow(angle = 35,length=unit(0.3,"cm")))+
  geom_text(data=df_rda_env,aes(x=df_rda_env[,1],y=df_rda_env[,2],
                                label=rownames(df_rda_env)),size=3.5,
            color="red",
            hjust="inward",
            vjust=0.5*(1-sign(df_rda_env[,2])))+
  theme(legend.position = "top") 
#蒙特卡洛置换检验评估环境因子与物种群落组成相关的显著性
write.csv(rdascore$sites,file="rda.sample.csv")
write.csv(rda$CCA$biplot,file="rda.env.csv")
write.csv(rdascore$species,file="rda.species.csv")
as.data.frame(rda$CCA$biplot[,1:2])
RDAS1 1]RDAS2 2]plotdata colnames(plotdata) "sample","RDAS1","RDAS2","group")rda1 1]/sum(rda$CCA$eig)*100,2)rda2 2]/sum(rda$CCA$eig)*100,2)
#RDA plot
plot_RDA   geom_label_repel(aes(RDAS1, RDAS2, label = sample), fill = "white", color = "black", box.padding = unit(0.6,"lines"),                    segment.colour = "grey50",label.padding = unit(0.35,"lines")) +geom_point(aes(fill = group, shape = group),size = 8)
+ scale_shape_manual(values = pch)+   scale_fill_manual(values = col)+   labs(title = "RDA Plot")
+    xlab(paste("CCA1 ( ",rda1,"%"," )", sep = "")) +    ylab(paste("CCA2 ( ",rda2,"%"," )", sep = "")) +   geom_segment(data = RDAE, aes(x = 0, y = 0, xend = RDAE[,1], yend = RDAE[,2]),                colour = "black", size = 1.2,                arrow = arrow(angle = 30, length = unit(0.4, "cm"))) +   geom_label_repel(data = RDAE, fill = "grey90",segment.colour = "black",                    aes(x = RDAE[,1], y = RDAE[,2], label = rownames(RDAE))) +   geom_vline(aes(xintercept = 0), linetype = "dotted") +   geom_hline(aes(yintercept = 0), linetype = "dotted") +   theme(panel.background = element_rect(fill = "white", colour = "black"),          panel.grid = element_blank(),         axis.title = element_text(color = "black", size = 18),         axis.ticks.length = unit(0.4,"lines"),         axis.ticks = element_line(color = "black"),         axis.line = element_line(colour = "black"),         axis.title.x = element_text(colour = "black", size = 18),         axis.title.y = element_text(colour="black", size = 18),         axis.text = element_text(colour = "black", size = 18),         legend.title = element_blank(),         legend.text = element_text(size = 18), legend.key = element_blank(),         plot.title = element_text(size = 22, colour = "black",                                    face = "bold", hjust = 0.5))png(filename="RDA.png",res=600,height=5400,width=7200)plot_RDAdev.off()#Permutation test envfit 999)r as.matrix(envfit$vectors$r)p as.matrix(envfit$vectors$pvals)env.p colnames(env.p) "r2","p-value")write.csv(as.data.frame(env.p),file="rdaenvfit.csv")
