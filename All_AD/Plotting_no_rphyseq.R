getwd()


library(tidyverse)
library(FSA)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(ggplot2)
library(phyloseq)
library(ape)
library(rstatix)



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



otu	<-	read.table(file	=	"otu_table.txt",	header	=	TRUE)
tax	<-	read.table(file	=	"taxonomy.tsv",	sep	=	'\t',	header	=	TRUE)
merged_file	<-	merge(otu,	tax,	by.x=c("OTUID"),	by.y=c("OTUID"))
write.table(merged_file,	file	=	"combined_otu_tax.tsv",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)


otu_table	=	read.csv("otu_matrix.csv",	sep=",",	row.names=1)
otu_table	=	as.matrix(otu_table)


taxonomy	=	read.csv("taxonomy.csv",	sep=",",	row.names=1)
taxonomy	=	as.matrix(taxonomy)

metadata	=	read.table("metadata.tsv",	row.names=1)
colnames(metadata)<-metadata[1,]
metadata <- metadata[-1,]

phy_tree	=	read_tree("tree.nwk")

OTU	=	otu_table(otu_table,	taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META	=	sample_data(metadata)

physeq <- phyloseq(OTU, TAX, META, phy_tree)


relaphyseq <- transform_sample_counts(physeq, function(otu) {otu/sum(otu)})







###relaseq stacked barplot

plot_bar(relaphyseq, fill = "Class") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack") +
  labs(x = "", y="Relative abundance") +
  ggtitle("Relative abundance stack bar plot by Treatment") +
  theme(axis.title = element_text(color="black", face="bold", size=10)) +
  theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )



plot_bar(relaphyseq, fill = "Order") +
  geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack") +
  labs(x = "", y="Relative abundance") +
  ggtitle("Relative abundance stack bar plot by Treatment") +
  theme(axis.title = element_text(color="black", face="bold", size=10)) +
  theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
  scale_color_brewer(palette = 'Set3')+
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )









###alpha diversity

richness <- estimate_richness(physeq)
richness$group <- metadata$group


plot_richness(physeq, sortby = META$group)
pair <- list (c("AD", "AP"), c("AD", "CON"), c("AP", "CON"))


###Chao1

dunn_Chao1 <- dunn_test(data = richness, Chao1 ~ group)


ggplot(data=richness, aes(x=group, y=Chao1)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Chao1', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (75, 220) +
  stat_compare_means(method = "anova", label.y = 210) +
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



###Shannon
dunn_Shannon <- dunn_test(data = richness, Shannon ~ group)


ggplot(data=richness, aes(x=group, y=Shannon)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Shannon', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (1.5, 4.5) +
  stat_compare_means(method = "anova", label.y = 4.4) +
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



###InvSimpson
dunn_InvSimpson <- dunn_test(data = richness, InvSimpson ~ group)

ggplot(data=richness, aes(x=group, y=InvSimpson)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Inverse Simpson', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (1.5, 25) +
  stat_compare_means(method = "anova", label.y = 24) +
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




###Fisher
dunn_Fisher <- dunn_test(data = richness, Fisher ~ group)

ggplot(data=richness, aes(x=group, y=Fisher)) +
  geom_boxplot(alpha = 0.5, aes(fill=group, col=group)) +
  labs(title= 'Fisher', x= ' ', y= ''
       # , tag = "A"
  ) +
  geom_point(aes(fill=group, col=group))+
  #  geom_jitter()+
  ylim (10, 30) +
  stat_compare_means(method = "anova", label.y = 29) +
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





###PCoA bray method
PCoA_bray <- ordinate(relaphyseq, 
                      method ="PCoA",
                      distance = "bray")

sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test
PCoA_bray_d <- phyloseq::distance(relaphyseq, method = "bray")
adonis2_bray <- adonis2(PCoA_bray_d ~ META$group, data = sampledf2) #0.007
adonis2(PCoA_bray_d ~ META$group, data = sampledf2)
anosim_bray <- anosim(PCoA_bray_d, grouping = META$group)


plot_ordination(relaphyseq, PCoA_bray, color = "group")+
  # stat_ellipse(colour = "transparent", level=0.3,
  #              alpha=0.3,
  #              geom = "polygon", aes(fill = group))+  
  stat_conf_ellipse(colour = "transparent", #level=0.3,
                    alpha=0.5,
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Bray",
       caption = paste("PERMANOVA = ", adonis2_bray$`Pr(>F)`[1],
                       ", ANOSIM = ", anosim_bray$signif))





###PCoA unweighted
PCoA_unifrac <- ordinate(relaphyseq, 
                         method ="PCoA",
                         distance = "unifrac")

sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test
PCoA_unifrac_d <- phyloseq::distance(relaphyseq, method = "unifrac")
adonis2_unifrac <- adonis2(PCoA_unifrac_d ~ META$group, data = sampledf2) #0.007
anosim_unifrac <- anosim(PCoA_unifrac_d, grouping = META$group)


plot_ordination(relaphyseq, PCoA_unifrac, color = "group")+
  # stat_ellipse(colour = "transparent", level=0.3, #level 얼마로 해야하는지 확인.
  #              alpha=0.3, #색진하기 정도
  #              geom = "polygon", aes(fill = group))+  
  stat_conf_ellipse(colour = "transparent", #level=0.3, #level 얼마로 해야하는지 확인.
                    alpha=0.5, #색진하기 정도
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Unweighted Unifrac",
       caption = paste("PERMANOVA = ", adonis2_unifrac$`Pr(>F)`[1],
                       ", ANOSIM = ", anosim_unifrac$signif))




###PCoA weighted
PCoA_wunifrac <- ordinate(relaphyseq, 
                          method ="PCoA", #defalt 가 DCA라는데
                          distance = "wunifrac")

sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test
PCoA_wunifrac_d <- phyloseq::distance(relaphyseq, method = "wunifrac")
adonis2_wunifrac <- adonis2(PCoA_wunifrac_d ~ META$group, data = sampledf2) #0.007
anosim_wunifrac <- anosim(PCoA_wunifrac_d, grouping = META$group)


plot_ordination(relaphyseq, PCoA_wunifrac, color = "group")+
  # stat_ellipse(colour = "transparent", level=0.3, #level 얼마로 해야하는지 확인.
  #              alpha=0.3, #색진하기 정도
  #              geom = "polygon", aes(fill = group))+  
  stat_conf_ellipse(colour = "transparent", #level=0.3, #level 얼마로 해야하는지 확인.
                    alpha=0.5, #색진하기 정도
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Weighted Unifrac",
       caption = paste("PERMANOVA = ", adonis2_wunifrac$`Pr(>F)`[1],
                       ", ANOSIM = ", anosim_wunifrac$signif))




library(mia)
library(lefser)
library(tidyverse)
library(SummarizedExperiment)
library(microbiomeMarker)


###LEfSe
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




###Heatmap + Taxonomy

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




###PhyloHeatmap 아래에 복사해둔 코드를 저장해야할 듯...
topN <- 20
sort(taxa_sums(relaphyseq), TRUE)
most_abundant_taxa <- sort(taxa_sums(relaphyseq), TRUE)[1:topN]
print(most_abundant_taxa)
#physeq20 <- prune_taxa(names(most_abundant_taxa), physeq)
relaphyseq20 <- prune_taxa(names(most_abundant_taxa), relaphyseq) #얘는 안됨....
#뭔가 이 prune 명령어가 subgroup 가능하게 하는 것 같음

table_20 <- t(as.matrix(relaphyseq20@otu_table@.Data))
tax_table_20 <- data.frame(relaphyseq20@tax_table)
tax_table_20

df_table_20 <- data.frame(table_20)

colnames(table_20) <-c('g__Lactobacillus', 's__Lactobacillus_intestinalis', 'g__Helicobacter', 's__uncultured_Bacteroidales', 
                       's__uncultured_Bacteroidales', 's__uncultured_Bacteroidales', 's__uncultured_bacterium', 's__uncultured_bacterium', 
                       'g__Muribaculaceae', 's__unidentified', 'g__Bacteroides', 'g__Bacteroides', 
                       's__uncultured_bacterium', 's__uncultured_bacterium', 'g__Lactobacillus', 'g__Lactobacillus', 
                       'c__Bacilli', 'c__Bacilli', 'c__Bacilli', 's__uncultured_Clostridiales')




###correlation matrix
library(corrplot)

correlation <- cor(table_20)
cor(df_table_20)
pheatmap(correlation, 
         color = colorRampPalette(c("skyblue", "white", "pink"))(50),
         border_color = "grey",
         scale = "none", 
         drop_levels = TRUE, #아무 영향 없음
         angle_col = "90",
         show_colnames = F
)

phyloseq::plot_heatmap(relaphyseq20)

corrplot(correlation,
         method = 'shade',
         type = 'lower'
)




###spearman
for (i in 1: 20){
  spearman_1 <-cbind(richness$Shannon, df_table_20[i])
  microbname <- colnames(df_table_20[i])
  colnames(spearman_1) <- c('Shannon', 'microb')
  
  paste(colnames(df_table_20[i]))                         
  
  fit <- lm(Shannon ~ microb, data = spearman_1)
  summary(fit)
  #dev.new(width=5, height=5, unit="px")
  
  ggplot(spearman_1,aes(microb, Shannon))+
    geom_point(col = 'pink')+
    #  stat_summary(fun.data=mean_cl_normal) +
    geom_smooth(method='lm', col = 'pink', fill = 'pink')+
    xlab('Relative Abundance')+
    labs(title = paste(colnames(df_table_20[i])), 
         caption = paste(" adj R^2 = ",signif(summary(fit)$adj.r.squared, 5),
                         #                   "\n", "Intercept =",signif(fit$coef[[1]],5 ),
                         #                   "Slope =",signif(fit$coef[[2]], 5), "\n",
                         ", P =",signif(summary(fit)$coef[2,4], 5)))
}
