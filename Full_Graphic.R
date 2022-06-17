#library에 없는게 있으면 보통 install.packages("")로 설치하면 되고 안되는 경우엔 검색하면 바로 나옴
# library(ape)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(dplyr)
library(ggtree)
# library(tidyverse)
# library(FSA)
library(RColorBrewer)
library(vegan)
library(rstatix)
# library(metagMisc)
# library(MicrobeR)
# library(mia)
# library(lefser)
# library(SummarizedExperiment)
# library(microbiomeMarker)
library(pheatmap)
# library(Graphlan)
# library("corrplot")

theme_set(
  theme_bw()+
  theme(axis.line = element_line(size=1),
          axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
))

getwd()
setwd("AD_microbiome_data/All_AD")
#이건 자기가 설정한 디렉토리 잘 찾아서 진행하기


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### STACKED BAR GRAPH #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#########+ function for graphing ######
# finding order... 보고 적당히 수기로 해야할 듯

plot_rela_group <- function (out, level = l) {
  level <- level + 1
  Levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

  r <- data.frame(out)
  r <- r[order(r$total),]
  order <- r[,level]
  
  # phylum
  con <- cbind(r[,level], r[,c("CON_mean")], c("CON")) 
  ad <- cbind(r[,level], r[,c("AD_mean")], c("AD"))
  ap <- cbind(r[,level], r[,c("AP_mean")], c("AP"))

  if (!sum( as.numeric(con[,2]))==1|!sum( as.numeric(ad[,2]))==1|!sum( as.numeric(ap[,2]))==1){
    con <- rbind(con, c("others", 1-sum( as.numeric(con[,2])), "CON"))
    ad <- rbind(ad, c("others", 1-sum( as.numeric(ad[,2])), "AD"))
    ap <- rbind(ap, c("others", 1-sum( as.numeric(ap[,2])), "AP"))
  }
  
  colnames(con) <- c(Levels[level], "Relative_abundance", "group")
  colnames(ad) <- c(Levels[level], "Relative_abundance", "group")
  colnames(ap) <- c(Levels[level], "Relative_abundance", "group")
  
  mr <- data.frame(rbind(con, ad, ap))
  
  mr[,1] <- factor(mr[,1], levels = c("others", order))
  mr[,2] <- as.numeric(mr[,2])
  mr[,3] <- factor(mr[,3], levels = c("CON", "AD", "AP"))

  
  ggplot(mr, aes(fill=mr[,1], y=Relative_abundance, x=group)) +
    geom_bar(position="stack", stat = "identity")+
    ggtitle("Relative abundance stack bar plot") +
    labs(x = "", y="Relative abundance") +
    ylim (0, 1) +##여기 숫자로 원하는 크기로 조정ㅎ가능
    
    theme(axis.title = element_text(color="black", face="bold", size=10)) +
    theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
    scale_color_brewer(palette = 'Set3')+
    scale_fill_brewer(palette = "Set3")+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.line = element_line(size=1),
          axis.ticks = element_line(size=1),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.subtitle = element_text(hjust = 0.5)
    )
}


##### @others가 있는 stacked bar 만들기####

nopalette_plot_rela_group <- function (out, level = l) {
  
  Levels = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species")
  
  r <- data.frame(out)
  r <- r[order(r$total),]
  order <- r[,level]
  
  # phylum
  con <- cbind(r[,level], r[,c("CON_mean")], c("CON")) 
  ad <- cbind(r[,level], r[,c("AD_mean")], c("AD"))
  ap <- cbind(r[,level], r[,c("AP_mean")], c("AP"))
  
  if (!sum( as.numeric(con[,2]))==1|!sum( as.numeric(ad[,2]))==1|!sum( as.numeric(ap[,2]))==1){
    con <- rbind(con, c("others", 1-sum( as.numeric(con[,2])), "CON"))
    ad <- rbind(ad, c("others", 1-sum( as.numeric(ad[,2])), "AD"))
    ap <- rbind(ap, c("others", 1-sum( as.numeric(ap[,2])), "AP"))
  }
  
  colnames(con) <- c(Levels[level], "Relative_abundance", "group")
  colnames(ad) <- c(Levels[level], "Relative_abundance", "group")
  colnames(ap) <- c(Levels[level], "Relative_abundance", "group")
  
  mr <- data.frame(rbind(con, ad, ap))
  
  mr[,1] <- factor(mr[,1], levels = c("others", order)) 
  mr[,2] <- as.numeric(mr[,2])
  mr[,3] <- factor(mr[,3], levels = c("CON", "AD", "AP"))
  
  ggplot(mr, aes(fill=mr[,1], y=Relative_abundance, x=group)) +
    geom_bar(position="stack", stat = "identity")+
    ggtitle("Relative abundance stack bar plot") +
    labs(x = "", y="Relative abundance") +
    ylim (0, 1) +##여기 숫자로 원하는 크기로 조정ㅎ가능
    
    theme(axis.title = element_text(color="black", face="bold", size=10)) +
    theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
#    scale_color_brewer(palette = 'Set3')+
#    scale_fill_brewer(palette = "Set3")+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.line = element_line(size=1),
          axis.ticks = element_line(size=1),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.subtitle = element_text(hjust = 0.5)
    )
}



out_6g[2,]$Genus=!"g__uncultured"
out_cut_6g <- subset(out_6g, !Genus==" g__uncultured")
out_cut_7s <- subset(out_7s, !Species==" s__uncultured_bacterium" & #49
                             !Species==" s__unidentified" &#34
                             !Species==" s__uncultured_organism" #29
                       )
out_weird <- subset(out_cut_7s, !Species==" s__uncultured_Clostridiales"&
                    !Species==" s__Clostridiales_bacterium"&
                    !Species==" s__Clostridium_sp."
                    )

cut_rela_group <- function (out, cut = c) {
  
  r <- data.frame(out)
  r <- r[rev(order(r$total)),]
  
  sum <- c("others", "others", "others", "others", "others", "others", "others")
  sum <- data.frame(sum)
  
    for (i in 8:ncol(r)) {
    sum <- rbind(sum, sum(as.numeric(others[,i])))
  }
  t_sum <- t(sum)
  colnames(t_sum) <- colnames(r)
  or <- rbind(r[1:cut,], t_sum)
  
  return(or)
}





plot_rela_group(out_2p, level = 2)
ggsave("relative_abd/Group_Phylum.png", width=6, height=5, units="in", device = "png")

plot_rela_group(out_3c , level = 3)
ggsave("relative_abd/Group_Class.png", width=6, height=5, units="in", device = "png")


nopalette_plot_rela_group(out_4o , level = 4)
ggsave("relative_abd/Sample_Order.png", width=10, height=5, units="in", device = "png")

nopalette_plot_rela_group(out_5f , level = 5)
ggsave("relative_abd/Sample_Family.png", width=10, height=5, units="in", device = "png")

nopalette_plot_rela_group(out_cut_6g , level = 6)
ggsave("relative_abd/Sample_Genus.png", width=14, height=5, units="in", device = "png")

nopalette_plot_rela_group(out_weird , level = 7)
ggsave("relative_abd/Sample_Family.png", width=14, height=5, units="in", device = "png")

cut_5f <- cut_rela_group(out_5f, cut = 10)
plot_rela_group(cut_5f, level = 5)

sum(as.numeric(data.frame(out_3c)$total))


# ########## ++ Order level ###########
# 
# plot_bar(relaphyseq, fill = "Order") +
#   geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack") +
#   labs(x = "", y="Relative abundance") +
#   ggtitle("Relative abundance stack bar plot by Treatment") +
#   theme(axis.title = element_text(color="black", face="bold", size=10)) +
#   theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
#   #scale_color_brewer(palette = 'Set3')+
#   #scale_fill_brewer(palette = "Set3")+
#   scale_y_continuous(expand = c(0,0))+
#   theme_bw()+
#   theme(axis.line = element_line(size=1),
#         axis.ticks = element_line(size=1),
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5, face="bold"),
#         plot.subtitle = element_text(hjust = 0.5)
#   )



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### ALPHA DIVERSITY #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##richness1 <- estimate_richness(rphyseq)
##richness2 <- estimate_richness(rphyseq2)
richness <- estimate_richness(rphyseq)

richness$group <- metadata$group
richness$group <- factor (richness$group, levels = c("CON", "AD", "AP"))

#간단하게 전체적으로 살펴보는 방법
##plot_richness(rphyseq, sortby = META$group)
plot_richness(physeq, sortby = META$group)

##ppplot에서 군별 paired 분석에 대한 pvalue 표기하려면 pair를 정해줘야함.
pair <- list (c("AD", "AP"), c("AD", "CON"), c("AP", "CON"))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + Chao1 ##############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
dunn_Chao1 <- dunn_test(data = richness, Chao1 ~ group)

#draw barplot

ggplot(data=richness, aes(x=group, y=Chao1)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Chao1', x= ' ', y= ''
      # , tag = "A"
       ) +
  geom_point(aes(fill=group, col=group))+
#  geom_jitter()+
  ylim (75, 220) + ##여기 숫자로 원하는 크기로 조정ㅎ가능
  stat_compare_means(method = "kruskal.test", label.y = 210) +  # Add global p-value
  stat_pvalue_manual(dunn_Chao1, #이게....
                     y.position = c(180, 200, 190)) +
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position= "none",
  )


ggsave("alpha_div/Chao1.png", width=3, height=3, units="in", device = "png")

# #one way anova
# one.way <- aov(Chao1 ~ META$group, data = richness)
# summary(one.way)
# 
# kruskal.test(Chao1 ~ META$group, data = richness)
# 
# #pairwise.wilcox.test(chao1_data$chao1, chao1_data$Group,
# #                    p.adjust.method = "BH")
# 
# dunnTest(Chao1 ~ META$group, data = richness,
#          method="bonferroni")




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + Shannon ##############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dunn_Shannon <- dunn_test(data = richness, Shannon ~ group)

#draw barplot

ggplot(data=richness, aes(x=group, y=Shannon)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Shannon', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (1.5, 4.5) + ##여기 숫자로 원하는 크기로 조정ㅎ가능
  stat_compare_means(method = "kruskal.test", label.y = 4.4) +  # Add global p-value
  stat_pvalue_manual(dunn_Shannon, 
                     y.position = c(3.8, 4.2, 4.0)) +
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position= "none",
  )

ggsave("alpha_div/Shannon.png", width=3, height=3, units="in", device = "png")




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + InvSimpson ##############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dunn_InvSimpson <- dunn_test(data = richness, InvSimpson ~ group)

#draw barplot

ggplot(data=richness, aes(x=group, y=InvSimpson)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Inverse Simpson', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (1.5, 25) + ##여기 숫자로 원하는 크기로 조정ㅎ가능
  stat_compare_means(method = "kruskal.test", label.y = 24) +  # Add global p-value
  stat_pvalue_manual(dunn_InvSimpson, 
                     y.position = c(18, 22, 20)) +
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position= "none",
  )

ggsave("alpha_div/InvSimpson.png", width=3, height=3, units="in", device = "png")




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + Fisher ##############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dunn_Fisher <- dunn_test(data = richness, Fisher ~ group)

#draw barplot

ggplot(data=richness, aes(x=group, y=Fisher)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Fisher', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (10, 30) + ##여기 숫자로 원하는 크기로 조정ㅎ가능
  stat_compare_means(method = "kruskal.test", label.y = 29) +  # Add global p-value
  stat_pvalue_manual(dunn_Fisher, 
                     y.position = c(26, 28, 27)) +
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position= "none",
  )

ggsave("alpha_div/Fisher.png", width=3, height=3, units="in", device = "png")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### BETA DIVERSITY #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#https://github.com/joey711/phyloseq/issues/1046

#set.seed(134) ###이거 뭐지

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + PCoA  - bray method #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###****calculating bray curtis distance matrix****##
PCoA_bray <- ordinate(relaphyseq, 
                         method ="PCoA", #defalt 가 DCA라는데
                         distance = "bray")

#making a data frame from the sample_data
sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test




PCoA_bray_d <- phyloseq::distance(relaphyseq, method = "bray")
adonis2_bray <- adonis2(PCoA_bray_d ~ META$group, data = sampledf2) #0.007
adonis2(PCoA_bray_d ~ META$group, data = sampledf2)
anosim_bray <- anosim(PCoA_bray_d, grouping = META$group)
##그리기...
#theme_set(theme_bw())

plot_ordination(relaphyseq, PCoA_bray, color = "group")+
  # stat_ellipse(colour = "transparent", level=0.3, #level 얼마로 해야하는지 확인.
  #              alpha=0.3, #색진하기 정도
  #              geom = "polygon", aes(fill = group))+  
  stat_conf_ellipse(colour = "transparent", #level=0.95, #level 얼마로 해야하는지 확인.
                    alpha=0.3, #색진하기 정도
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Bray",
       caption = paste("PERMANOVA = ", adonis2_bray$`Pr(>F)`[1],
                           ", ANOSIM = ", anosim_bray$signif))

ggsave("beta_div/bray.png", width=4, height=3, units="in", device = "png")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + PCoA  - unweighted #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#calculating unifrac curtis distance matrix
PCoA_unifrac <- ordinate(relaphyseq, 
                      method ="PCoA", #defalt 가 DCA라는데
                      distance = "unifrac")

#making a data frame from the sample_data

sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test

PCoA_unifrac_d <- phyloseq::distance(relaphyseq, method = "unifrac")
adonis2_unifrac <- adonis2(PCoA_unifrac_d ~ META$group, data = sampledf2) #0.007
anosim_unifrac <- anosim(PCoA_unifrac_d, grouping = META$group)
##그리기...
#theme_set(theme_bw())

plot_ordination(relaphyseq, PCoA_unifrac, color = "group")+
  # stat_ellipse(colour = "transparent", level=0.3, #level 얼마로 해야하는지 확인.
  #              alpha=0.3, #색진하기 정도
  #              geom = "polygon", aes(fill = group))+  
  stat_conf_ellipse(colour = "transparent", #level=0.3, #level 얼마로 해야하는지 확인.
                    alpha=0.3, #색진하기 정도
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Unweighted Unifrac",
       caption = paste("PERMANOVA = ", adonis2_unifrac$`Pr(>F)`[1],
                       ", ANOSIM = ", anosim_unifrac$signif))

ggsave("beta_div/unifrac.png", width=4, height=3, units="in", device = "png")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + PCoA  - weighted #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#calculating wunifrac curtis distance matrix
PCoA_wunifrac <- ordinate(relaphyseq, 
                         method ="PCoA", #defalt 가 DCA라는데
                         distance = "wunifrac")

#making a data frame from the sample_data

sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test

PCoA_wunifrac_d <- phyloseq::distance(relaphyseq, method = "wunifrac")
adonis2_wunifrac <- adonis2(PCoA_wunifrac_d ~ META$group, data = sampledf2) #0.007
anosim_wunifrac <- anosim(PCoA_wunifrac_d, grouping = META$group)
##그리기...
#theme_set(theme_bw())

plot_ordination(relaphyseq, PCoA_wunifrac, color = "group")+
  # stat_ellipse(colour = "transparent", level=0.3, #level 얼마로 해야하는지 확인.
  #              alpha=0.3, #색진하기 정도
  #              geom = "polygon", aes(fill = group))+  
  stat_conf_ellipse(colour = "transparent", #level=0.3, #level 얼마로 해야하는지 확인.
                    alpha=0.3, #색진하기 정도
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Weighted Unifrac",
       caption = paste("PERMANOVA = ", adonis2_wunifrac$`Pr(>F)`[1],
                       ", ANOSIM = ", anosim_wunifrac$signif))

ggsave("beta_div/wunifrac.png", width=4, height=3, units="in", device = "png")






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### FINDING SIGNIFICANT SPECIES #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + LEfSe #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# relasum <- makeTreeSummarizedExperimentFromPhyloseq(relaphyseq)
# unique(relasum$SampleType)
# 
# SummarizedExperiment(relaphyseq)
# lefser(relasum, groupCol = "group")
#*****lefse에서 2개 군만 선별해서 돌리ㄹ는 법을 찾아야함-> 해결 못했음****###
#*****lefse와 maaslin에서도 동일 data set으로 돌아가려면 rarefaction을 R로 넘어오기 전에 해야할듯??****###
#Summarized Experiment로 아예 그룹을 빼고 physeq를 만들어봐야할 긋?



lefse <- run_lefse(relaphyseq,"group")

lefse_data <-data.frame(marker_table(lefse))
plot_ef_bar(lefse)+
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5))


ggsave("other/LEfSe.png", width=10, height=10, units="in", device = "png")




####+tax로 합치고 돌리면####
rela_tax <- tax_glom(relaphyseq, taxrank = "Species")

lefseOTU <- run_lefse(rela_tax,"group")

lefse_data_OTU <-data.frame(marker_table(lefseOTU))
plot_ef_bar(lefseOTU)+
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5))


ggsave("other/LEfSeOTU.png", width=10, height=6, units="in", device = "png")

####+tip으로 합치고 돌리면####
rela_tip <- tip_glom(relaphyseq)

lefseOTUtip <- run_lefse(rela_tip,"group")

lefse_data_OTUtip <-data.frame(marker_table(lefseOTUtip))
plot_ef_bar(lefseOTUtip)+
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5))


ggsave("other/LEfSeOTUtip.png", width=10, height=7, units="in", device = "png")








# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### @ Cladogram #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


##*****cladogram 오류 햐결 안됨...****##
# library(ggtree)
# ####https://yulab-smu.top/treedata-book/chapter4.html
# 
# fan.angle ()
# ggtree(relaphyseq, layout = 'circular', 
#        branch.length = 3, ladderize) #우리가 원하는게 circular 인건 맞는 듯
# 
# 
# data(lefse)
# lefse_small <- phyloseq::subset_taxa(
#   lefse,
#   Phylum %in% c("Firmicutes")
# )
# plot_cladogram(lefse, color = c(AD ="red", AP = "blue", CON = "green"), only_marker = TRUE,  clade_label_level = 4)
# 
# 
# plot_cladogram(mm_lefse, color = c(Healthy = "darkgreen", Tumor = "red")) +
#   theme(plot.margin = margin(0, 0, 0, 0))

rela_otu <- data.frame(relaphyseq@otu_table@.Data)
rela_otu$OTUID <- rownames(rela_otu)
m_rela_table <- merge(rela_otu, tax, by="OTUID")
m_rela_table <- m_rela_table[-1]
m_rela_table <- m_rela_table[-18]

write.table(m_rela_table,	file	=	"phyloseq/m_rela_table.tsv",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)
# 여기서 부터는 수기로 편집.
#           Taxonomy 열 -> ; 로 부분되어 있는 것 |로 구분 
# Group 행
# Sample 행





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + Heatmap + taxonomy #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

plot_heatmap(
  lefse,
  transform = c("log10", "log10p", "identity"),
  cluster_marker = FALSE,
  cluster_sample = FALSE,
  markers = NULL,
  label_level = 1,
  max_label_len = 60,
  sample_label = FALSE,
  scale_by_row = FALSE,
  annotation_col = NULL,
  'group'
)

plot_ef_dot(lefse)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### @ MAASLIN #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#
#library("devtools") #github에서 다운받으려면 필요한 것 
#LEfSe identifies those data features that are distinct between a pair of metadatums (e.g. differences between two sampling sites, two clinical outcomes, two biochemical markers, two modalities, etc.).  MaAsLin extends the functionality of LEfSe to identify associations between data features and multiple metadata factors, which can be discrete and/or continuous and can include time series data

# library("yingtools2")
# library("tidyverse")
# library("Maaslin2")
# 
# #https://github.com/biobakery/maaslin#markdown-header-input-files
# #*****위의 형식을 보고 input file 만들기****
# 
# maaslin <- Maaslin2(input_data = 'MaAsLin/feature_table_L7.txt',# sample과 feature만 있는 것, 이미 rarefaction 되었음
#                     input_metadata = 'MaAsLin/metadata.txt', #group 정보 등
#                     analysis_method = 'LM', #기본 값이 LM인데 다른 방법들의 차이를 확인할수가 없음
#                     output = 'MaAsLin/maaslin', #output folder 이름
#                     reference = c('Group;AD;AP;CON')) #다변수인 경우 나열할 것.
# 
# 
# 
# 
# 
# 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### CORRELATION GRAPHS  #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##################### + heatmap  #############################

#if there's NA in the data, can't graph. so remove NAs
coord_NAx <- data.frame(coord_bug_bug[complete.cases(coord_bug_bug),])


#colnames(coord_NAx) <- c("Phylum", "Genus", "R", "p", "fdr", "p_", "fdr_")
coord_NAx$R <- as.numeric(coord_NAx$R) #needs to be numeric to be colored in geom_tile

#sub_coord_NAx <- coord_NAx[c(1:10,80:89,159:168),]


#그림 오래 걸림...
ggplot(sub_coord_NAx, aes(x=Phylum, y=Genus, fill=R))+
  geom_tile()+
  scale_fill_gradient2(low="skyblue", mid = "white", high = "pink")+
  geom_text( aes(label= ifelse(fdr_=="ns",p_,fdr_)), 
             color = ifelse(sub_coord_NAx$fdr_=="ns","grey","black"),
             alpha = ifelse(sub_coord_NAx$p_=="ns", 0, 1)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  
  coord_fixed() #keep them square, not rectangle
  

ggsave("output_correlation/sub.png", width=20, height=20, units="in", device = "png")

rowSums(is.na(coord_bug_behav[,3:5]))==0




#if there's NA in the data, can't graph. so remove NAs
coord_NAx <- data.frame(coord_bug_self[complete.cases(coord_bug_self),])

coord_NAx$R <- as.numeric(coord_NAx$R) #needs to be numeric to be colored in geom_tile

# sub_coord_NAx <- coord_NAx[c(1:10,80:89,159:168),]


#그림 오래 걸림...

ggplot(coord_NAx, aes(x=name1, y=name2, fill=R))+
  geom_tile()+
  scale_fill_gradient2(low="skyblue", mid = "white", high = "pink")+
  geom_text( aes(label= ifelse(fdr_=="ns",p_,fdr_)), 
             color = ifelse(coord_NAx$fdr_=="ns","grey","black"),
             alpha = ifelse(coord_NAx$p_=="ns", 0, 1)) +
  #geom_text( aes(label= ifelse(name1==name2, name, "")))+
  scale_y_discrete()+
  scale_x_discrete(limits=rev)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_fixed() #keep them square, not rectangle
#글씨는 따로 넣어야 할 듯...

ggsave("output_correlation/bug_self_10.png", width=5, height=5, units="in", device = "png")

rowSums(is.na(coord_bug_behav[,3:5]))==0




#MicrobeR::Microbiome.Heatmap	
#microbiomeMarker::plot_heatmap	








#***correlation matrix****

library("corrplot")

correlation <- cor(df_otu_20)
# cor(df_otu_20)
pheat<- pheatmap(correlation, 
         color = colorRampPalette(c("skyblue", "white", "pink"))(50),
         border_color = "grey",
                   scale = "none", 
                   drop_levels = TRUE, #아무 영향 없음
                   angle_col = "90",
                   show_colnames = F
)

#ㅇㅒ는 ggsave로 저장이 안됨. 그래서 저장 함수를 설계 해줘야함

save_pheatmap_pdf <- function(x, filename, width=8, height=5.5) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
   }

save_pheatmap_pdf(pheat, "tax/heatmap.pdf")


# #이건 다른 ㅇ히트맵.. rela, correlation 아님
# phyloseq::plot_heatmap(relaphyseq20)
# #***얘는 다른 heatmap에서 활용할 수 있을 듯 ***
# 
# #*** 얘는 삼각형이 가능한 heatmap***
# corrplot(correlation,
#          method = 'shade',
#          type = 'lower'
#          )
# 
# ggsave("Other/correlation.png", width=10, height=10, units="in", device = "png")
# 
#   

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############# + Spearman correlation #################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#1 

for (i in 1: 20){
spearman_1 <-cbind(richness$Shannon, df_otu_20[i])
microbname <- colnames(df_otu_20[i])
colnames(spearman_1) <- c('Shannon', 'microb')

paste(colnames(df_otu_20[i]))                         
                          
fit <- lm(Shannon ~ microb, data = spearman_1)
summary(fit)
#dev.new(width=5, height=5, unit="px")

ggplot(spearman_1,aes(microb, Shannon))+
  geom_point(col = 'pink')+
  #  stat_summary(fun.data=mean_cl_normal) +
  geom_smooth(method='lm', col = 'pink', fill = 'pink')+
  xlab('Relative Abundance')+
  labs(title = paste(colnames(df_otu_20[i])), 
       caption = paste(" adj R^2 = ",signif(summary(fit)$adj.r.squared, 5),
                         #                   "\n", "Intercept =",signif(fit$coef[[1]],5 ),
                         #                   "Slope =",signif(fit$coef[[2]], 5), "\n",
                         ", P =",signif(summary(fit)$coef[2,4], 5)))
## adj R는 음수가 나오기도 한다....
#https://www.researchgate.net/post/Interpretation_of_negative_Adjusted_R_squared_R22
#이런 경우에는 그냥 0으로 간주한다.

ggsave(width = 3.5, height = 3.5, dpi = 300, units = "in",
       filename = paste(i, colnames(df_otu_20[i]), ".png"),
       path = "tax/",
       device = "png"
       )
}








