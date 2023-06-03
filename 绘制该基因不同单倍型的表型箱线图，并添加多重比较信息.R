################################################################################
# 分析目标：绘制该基因不同单倍型的表型箱线图，并添加多重比较信息???
# 分析思路???
#            1、提取这些基因中样本数大???5%的单倍型
#            2、获取这些单倍型对应的群体信息和表型信息
#            3、多重比???
#            4、绘制箱线图
#            5、添加多重比较结???
#
################################################################################


library(pacman)
p_load(tidyverse,data.table,ggpubr,agricolae,DescTools,dplyr,generics)
rm(list = ls())
dev.off()
gene <- c('Os01g0838100')
pop <- fread('pop.CSV')
tgw <- fread('PHE.csv')
#water <- fread('water-phe.CSV')
g=gene[1]
# for (g in gene) {
file1 <- as.character(paste(g,'.csv',sep = ''))
df <- fread(file1,sep = ',')
df_tgw <- tgw %>% inner_join(df,by = 'ID')
df_tgw1 <- df_tgw %>% inner_join(pop,by = 'ID')
temp <- df_tgw1[,-5]
temp$K5 <- rep('Whole',nrow(temp))
df_tgw2 <- rbind(temp,df_tgw1)
#df_water <- water %>% inner_join(df,by = 'ID')
#df_water1 <- df_water %>% inner_join(pop,by = 'ID')
#temp1 <- df_water1[,-5]
#temp1$K5 <- rep('Whole',nrow(temp1))
#df_water2 <- rbind(temp1,df_water1)
df_tgw2$Hap=as.factor(df_tgw2$Hap)
df_tgw2$K5=as.factor(df_tgw2$K5)
df_tgw2$K5 <- factor(df_tgw2$K5,levels = c('Whole','XI','GJ','ADM','Aus','Bas'),ordered = T)
#df_water2$K5 <- factor(df_water2$K5,levels = c('Whole','XI','GJ','ADM','Aus','Bas'),ordered = T)


# file2 <- as.character(paste(g,'sand.csv',sep = '_'))
# file3 <- as.character(paste(g,'water.csv',sep = '_'))
# fwrite(df_sand2,file2,quote = F,row.names = F)
# fwrite(df_water2,file3,quote = F,row.names = F)
 # }



# plot ==========================

p <- ggboxplot(df_tgw2, x = "K5", y = "TGW",
               bxp.errorbar = T,
               color = "Hap", palette = "nejm")+xlab('3K')
p
p1 <- ggboxplot(df_tgw2, x = "K5", y = "PL",
             bxp.errorbar = T,
             color = "Hap", palette = "nejm")+xlab('3K')
p1
#p2 <- ggboxplot(df_water2, x = "K5", y = "CLs",
 #               bxp.errorbar = T,
#               color = "Hap", palette = "nejm")+xlab('water')
#p3 <- ggboxplot(df_water2, x = "K5", y = "MLs",
 #               bxp.errorbar = T,
  #             color = "Hap", palette = "nejm")+xlab('water')
#p3
#grid.arrange(p,p1,
 #             p2,p3,
  #            ncol = 2)
ggarrange(p,p1,ncol = 1,nrow = 2,
          common.legend = T,
          legend = "top")
# no outlier ==============================
p <- ggboxplot(df_tgw2, x = "K5", y = "TGW",
               outlier.shape = NA,
               bxp.errorbar = T,
               color = "Hap", palette = "nejm")+xlab('3K')
p
p1 <- ggboxplot(df_tgw2, x = "K5", y = "PL",
               outlier.shape = NA,
               bxp.errorbar = T,
             color = "Hap", palette = "nejm")+xlab('3K')
#p2 <- ggboxplot(df_water2, x = "K5", y = "CLs",
 #               outlier.shape = NA,
  #              bxp.errorbar = T,
   #            color = "Hap", palette = "nejm")+xlab('water')
#p3 <- ggboxplot(df_water2, x = "K5", y = "MLs",
 #               outlier.shape = NA,
  #              bxp.errorbar = T,
   #            color = "Hap", palette = "nejm")+xlab('water')
ggarrange(
           p1,#p3, #Os07g
          p,#p2, # Os08g
          ncol = 1,nrow = 2,
          common.legend = T,
          legend = "top")


# 多重比较 ==========



model <- aov(TGW~Hap,data=df_tgw1)
## Os07 =======================
# Os07-sand
a=df_tgw2 %>% filter(K5 == 'Whole') %>% filter(Hap != 'Hap4')
b=df_tgw2 %>% filter(K5 == 'XI') %>% filter(Hap != 'Hap3' & Hap != 'Hap4')
c=df_tgw2 %>% filter(K5 == 'GJ') %>% filter(Hap != 'Hap1' & Hap != 'Hap4')

mod1 <- aov(TGW~Hap,data=a)
res <- LSD.test(mod1, 'Hap', p.adj = 'bonferroni')
print(res$groups)
mod2 <- aov(TGW~Hap,data=b)
res1 <- LSD.test(mod2, 'Hap', p.adj = 'bonferroni')
print(res1$groups)
mod3 <- aov(TGW~Hap,data=c)
res2 <- LSD.test(mod3, 'Hap', p.adj = 'bonferroni')
print(res2$groups)

# Os07-water
#d=df_water2 %>% filter(K5 == 'Whole') 
#e=df_water2 %>% filter(K5 == 'XI') %>% filter(Hap != 'Hap3')
#g=df_water2 %>% filter(K5 == 'GJ') %>% filter( Hap != 'Hap4')
#h=df_water2 %>% filter(K5 == 'ADM') %>% filter(Hap != 'Hap3' & Hap != 'Hap4')
#mod4 <- aov(MLs~Hap,data=d)
#res3 <- LSD.test(mod4, 'Hap', p.adj = 'bonferroni')
#print(res3$groups)
#mod5 <- aov(MLs~Hap,data=e)
#res4 <- LSD.test(mod5, 'Hap', p.adj = 'bonferroni')
#print(res4$groups)
#mod6 <- aov(MLs~Hap,data=g)
#res5 <- LSD.test(mod6, 'Hap', p.adj = 'bonferroni')
#print(res5$groups)
#mod7 <- aov(MLs~Hap,data=h)
#res6 <- LSD.test(mod7, 'Hap', p.adj = 'bonferroni')
#print(res6$groups)
## Os08 =======================
# sand
#a=df_sand2 %>% filter(K5 == 'Whole') 
#b=df_sand2 %>% filter(K5 == 'XI') %>% filter(Hap != 'Hap2')
#c=df_sand2 %>% filter(K5 == 'GJ') %>% filter(Hap != 'Hap4')

#mod1 <- aov(CLs~Hap,data=a)
#res <- LSD.test(mod1, 'Hap', p.adj = 'bonferroni')
#print(res$groups)
#mod2 <- aov(CLs~Hap,data=b)
#res1 <- LSD.test(mod2, 'Hap', p.adj = 'bonferroni')
#print(res1$groups)
#mod3 <- aov(CLs~Hap,data=c)
#res2 <- LSD.test(mod3, 'Hap', p.adj = 'bonferroni')
#print(res2$groups)

# water
#d=df_water2 %>% filter(K5 == 'Whole') 
#e=df_water2 %>% filter(K5 == 'XI') %>% filter(Hap != 'Hap2')
#f=df_water2 %>% filter(K5 == 'GJ') %>% filter( Hap != 'Hap4')
#mod1 <- aov(CLs~Hap,data=d)
#res <- LSD.test(mod1, 'Hap', p.adj = 'bonferroni')
#print(res$groups)
#mod2 <- aov(CLs~Hap,data=e)
#res1 <- LSD.test(mod2, 'Hap', p.adj = 'bonferroni')
#print(res1$groups)
#mod3 <- aov(CLs~Hap,data=f)
#res2 <- LSD.test(mod3, 'Hap', p.adj = 'bonferroni')
#print(res2$groups)

# check ====================================
a=df_tgw2 %>% filter(Hap=="Hap4") %>% filter(K5 =="Xian")
b=a[,-2]
c=na.omit(b)
ggboxplot(c, x = "K5", y = "TGW",color = "Hap", palette = "nejm")

# 8g
a=df_tgw2 %>%  filter(K5 =="Adm")
b=a[,-2]
c=na.omit(b)
ggboxplot(c, x = "K5", y = "TGW",color = "Hap", palette = "nejm")

