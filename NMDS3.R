library(vegan)
library(ggplot2)
library(ggthemes)
# Load data
setwd("D:/M/R/BYD")
otu1 <- read.csv("spe4.csv",sep=",",header = T,row.names = 1,check.names = F)
otu<-t(otu1)
otu.distance <- vegdist(otu, method = 'jaccard')
#NMDS????????????vegan???е?metaMDS????
df_nmds <- metaMDS(otu.distance, k = 3)
#?????鿴??????עstress??points??species????ָ??
summary(df_nmds)
#Ӧ��????ֵ??<=0.2???���
df_nmds_stress <- df_nmds$stress
df_nmds_stress
#?????۲?ֵ??????????????????֮???Ĺ?ϵ????û?е??ֲ????߶ν?Զλ?ñ?ʾ?????ݿ???ʹ??NMDS????
stressplot(df_nmds)
#??ȡ??ͼ????
df_points <- as.data.frame(df_nmds$points)
#????samp1es??��
df_points$samples <- row.names(df_points)
#?޸?????
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)
p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+#ָ?????ݡ?X?ᡢY??
  geom_point(size=3)+#???Ƶ?ͼ???趨??С
  theme_bw()#????
p
#?????????ļ?
group <- read.csv("group4.csv",sep=",",header = T,check.names = F)
#?޸?????
colnames(group) <- c("samples","group")
#????ͼ???ݺͷ????ϲ?
df <- merge(df_points,group,by="samples")
head(df)
#ʹ??ggplot2????ͼ
color=c("#2878B5","#F18F01","#CD3700","#FFACACFF","#238E23")#??ɫ??��
ggplot(df, aes(x=NMDS1,y=NMDS2, fill=group, color=group)) +
  geom_point(size=4) +
  stat_ellipse(level=0.95,size = 0.8)+
  scale_color_manual(values=c("#2878B5", "#F18F01","#CD3700"))+
  labs(y="NMDS2",x="NMDS1")+ 
  theme(axis.text.y = element_text(size = 1))+
  geom_hline(yintercept=0,linetype = 3,size = 0.8,color="gray") +
  geom_vline(xintercept=0,linetype = 3,size = 0.8,color="gray")+  theme_bw(base_size=16)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))

p3<-ggplot(data=df,aes(x=NMDS1,y=NMDS2, fill=group, color=group))+#ָ?????ݡ?X?ᡢY?ᣬ??ɫ
  theme_bw()+#????????
  geom_point(size=3,shape=16)+#???Ƶ?ͼ???趨??С
  theme(panel.grid = element_blank())+#??????????
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+#ͼ??????
 #ȥ??ͼ??????
  labs(y="NMDS2",x="NMDS1")+#??x??y????????Ϊ???׶?
  scale_color_manual(values=c("#F18F01","#2878B5","#CD3700" )) +#??????ɫ????
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),#?޸?x???̶ȱ?ǩ?ı?
        axis.text.x=element_text(size=10))+#?޸?y???̶ȱ?ǩ?ı?
ggtitle(paste('Stress=',round(df_nmds_stress, 3)))

p3
library(vegan)
#anosim
anosim.result<-anosim(otu.distance,group$group,permutations = 999)

summary(anosim.result)

#adonis
adonis_result_dis = adonis2(otu.distance~group, group, permutations = 999)
adonis_result_dis

### 配对Adonis确定两两分组之间对物种组成差异的影响
library(pairwiseAdonis)
dune.pairwise.adonis <- pairwise.adonis(x=otu, factors=group$group, 
                                        sim.function = "vegdist",
                                        sim.method = "jaccard",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

library(ggpubr)
library(patchwork)
tab2 <- ggtexttable(dune.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                    theme = ttheme("blank"))%>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)%>%
  tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)  
tab2
tab1 <- data.frame(dune.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")])
  write.table(tab1, file = "permonova.txt",sep = "\t",row.names = F,col.names =F,quote = F)

#Shannon
  shannon <- read.csv("shannon.csv",sep=",",header = T,row.names = 1,check.names = F) 
  x <- shannon[,1]
  dune.pairwise.adonis <- pairwise.adonis(x=x, factors=shannon$group, 
                                          sim.function = "vegdist",
                                          sim.method = "jaccard",
                                          p.adjust.m = "BH",
                                          reduce = NULL,
                                          perm = 999)
  
  library(ggpubr)
  library(patchwork)
  tab2 <- ggtexttable(dune.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                      theme = ttheme("blank"))%>%
    tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)%>%
    tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)  
  tab2
  tab1 <- data.frame(dune.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")])
  write.table(tab1, file = "permonova.txt",sep = "\t",row.names = F,col.names =F,quote = F)
  
  

  #—————OLD—————————PERMANNOVA————————————

adonis(otu ~ Group,  # 矩阵 ~ 表型变量
       otu = group,  # 表型表格
       distance = "bray",   # 距离算法
       permutations = 999)  # 排列次数
# 循环所有算法
algo = c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis")
p=c()
for(title in algo)
{
  res = adonis2(otu~ Group,otu=group, 
               distance = title, permutations = 999)
  # 提取P value并保存
  p = c(p, res$aov.tab$Pr[1])
  print(title)
}
perm = data.frame(algo, p)
#第一个参数是需要导出的数据名称
#第二个参数是导出后新文件的名称
#第三个参数是指文件的分隔符
#导出数据和导入数据的参数类似，只是所使用的函数不同
write.table(dune.pairwise.adonis, "D:/M/R/BYD/pairwise_permonova.csv", sep=",")
#当然也可以直接使用write.csv()函数来实现
write.csv(mydata,"c:/mydata.csv")
#如果不想保存行名，可以设置row.names=F
#如果想在原文件后直接添加数据，可以使用append=T
write.csv(mydata,"c:/mydata.csv",row.names=F, append=T)