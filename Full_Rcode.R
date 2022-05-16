#library에 없는게 있으면 보통 install.packages("")로 설치하면 되고 안되는 경우엔 검색하면 바로 나옴
library(tidyverse)
library("FSA")
library("RColorBrewer")
library("ggpubr")
library(vegan)
library("ggplot2")
library("phyloseq")
library("ape")

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


########################################################################
###########################IMPORT DATA##################################
########################################################################
#https://nephele.niaid.nih.gov/user_guide_phyloseq_tutorial/
#https://github.com/joey711/phyloseq/issues/821
#아래 깃허브 주소에 pdf 올려주신 분이 있는데 그거 따라하면 됨.

#******* 수기로 #OTU ID -> OTUID로 변경해야함, 또한 taxonomy라고 써져있는 부분 지워야함.*******
#	read	in	OTU	table	
otu	<-	read.table(file	=	"phyloseq/DADA2_table.txt",	header	=	TRUE)
head(otu)

#******* 역시 수기로 OTUID, not Feature ID*******
#	read	in	taxonomy	table #taxonomy 에 Kingdom... 이런 식으로 일일히 column name 작성해야함
tax	<-	read.table(file	=	"phyloseq/taxonomy.tsv",	sep	=	'\t',	header	=	TRUE)
head(tax)

#	merge	files	
merged_file	<-	merge(otu,	tax,	by.x	=	c("OTUID"),	by.y=c("OTUID"))
head(merged_file)

#	note:	number	of	rows	should	equal	your	shortest	Sile	length,	drops	taxonomy	for	OTUs	that	don’t	exist	in	your	OTU	table
#	output	merged	.txt	Pile
write.table(merged_file,	file	=	"phyloseq/combined_otu_tax.tsv",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)

#	It	seems	tedious	but	you	need	to	open	the	merged	.txt	file	in	excel	and	split into	two	files:	one	for	taxonomy	(containing	only	the	columns	OTUID	and taxonomic	info)	and	the	other	for	the	OTU	matrix	(containing	only	OTUID	and abundances	in	each	sample).	Note:	for	the	taxonomy	file,	you	need	to	use	data —>	text-to-columns	in	Excel	and	separate	on	semicolon	to	get	columns	for kingdom,	phylum,	class,	etc…	once	you	make	these	two	separate	files	in	excel, save	each	as	a	.csv
#******* OTUID, taxonomic	info -> taxonomy.csv 로 저장, taxonomy ;->,로 ㅎ변경하여 column으로 만들기*******
#******* OTUID, abundance -> OTU	matrix.csv로 저장*******


#	Step	5,	Finally,	upload	all	of	your	files	into	phyloseq	in	R!



#	read	in	otu	table
otu_table	=	read.csv("phyloseq/otu_matrix.csv",	sep=",",	row.names=1)
otu_table	=	as.matrix(otu_table)

#	read	in	taxonomy
#	seperated	by	kingdom	phylum	class	order	family	genus	species
# taxonomy 에 Kingdom... 이런 식으로 일일히 column name 작성해야함
taxonomy	=	read.csv("phyloseq/taxonomy.csv",	sep=",",	row.names=1)
taxonomy	=	as.matrix(taxonomy)

#	read	in	metadata	
#	variables	=	???
metadata	=	read.table("phyloseq/metadata.tsv",	row.names=1)
colnames(metadata)<-metadata[1,]
metadata <- metadata[-1,]

#	read	in	tree
phy_tree	=	read_tree("phyloseq/tree.nwk")

#	import	as	phyloseq	objects
OTU	=	otu_table(otu_table,	taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META	=	sample_data(metadata)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

#	make	sure	files	have	the	same	sample	names	
sample_names(OTU)
sample_names(META)

#	merge	into	one	phyloseq	object
physeq <- phyloseq(OTU, TAX, META, phy_tree)
physeq

#	Now,	continue	to	analysis	in	phyloseq!




########Ready for Microbiome########





########################################################################
##################### Rarefy feature table #############################
########################################################################

#rarefy란 같은 수의 read로 맞추는 것 
#기존의 qiime에서 가장 작은 수의 sample이 rarefaction(?) graph에서 평평하게 유지ㅎ되는지 확인해야함
#*******우리 샘플의 경우 13410개가 최소였음*******
rphyseq <- rarefy_even_depth(physeq) #기본값이 최솟값으로 맞추는 것
rphyseq2 <- rarefy_even_depth(physeq, sample.size = 13410) #기본값이 최솟값으로 맞추는 것





########################################################################
##################### relative abundance #############################
########################################################################

#relative abundance로 변환
relaphyseq <- transform_sample_counts(rphyseq, function(otu) {otu/sum(otu)}) #relative abundance 구하는 것으로 예상



########################################################################
#****stacked bar graph of relative abundance***#
########################################################################

#***group 합ㄱ치는게 안됨 ***#
plot_bar(relaphyseq, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y="Relative abundance") +
  ggtitle("Relative abundance stack bar plot by Treatment") +
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  theme(axis.title = element_text(color="black", face="bold", size=10)) +
  theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )


plot_bar(relaphyseq, fill = "Class") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack") +
  labs(x = "", y="Relative abundance") +
  ggtitle("Relative abundance stack bar plot by Treatment") +
  theme(axis.title = element_text(color="black", face="bold", size=10)) +
  theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
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
  #scale_color_brewer(palette = 'Pastel1')+
  #scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )

plot_bar(relaphyseq, fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") +
  labs(x = "", y="Relative abundance") +
  ggtitle("Relative abundance stack bar plot by Treatment") +
  theme(axis.title = element_text(color="black", face="bold", size=10)) +
  theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
  #scale_color_brewer(palette = 'Pastel1')+
  #scale_fill_brewer(palette = "Pastel1")+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )

#****stacked bar graph of relative abundance PER GROUP ***#
plot_bar(relaphyseq, fill = "Class") +
  geom_bar(aes(color = Class, fill = Class), stat = "identity", position = "stack", 
           group = META$group) +
  labs(x = "", y="Relative abundance") +
  ggtitle("Relative abundance stack bar plot by Treatment") +
  theme(axis.title = element_text(color="black", face="bold", size=10)) +
  theme(plot.title = element_text(color="black", face = "bold", size =12, hjust = 0.5))+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )
























########################################################################
##################### alpha diversity #############################
########################################################################

richness <- estimate_richness(rphyseq)
richness2 <- estimate_richness(rphyseq2)
richness0 <- estimate_richness(physeq)
richness$group <- metadata$group

#간단하게 전체적으로 살펴보는 방법
plot_richness(rphyseq, sortby = META$group)

##ppplot에서 군별 paired 분석에 대한 pvalue 표기하려면 pair를 정해줘야함.
pair <- list (c("AD", "AP"), c("AD", "CON"), c("AP", "CON"))



##########################################
############### Chao1 ##############
##########################################

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
  stat_compare_means(method = "anova", label.y = 210) +  # Add global p-value
  stat_pvalue_manual(dunn, 
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




##########################################
############### Shannon ##############
##########################################

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
  stat_compare_means(method = "anova", label.y = 4.4) +  # Add global p-value
  stat_pvalue_manual(dunn, 
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




##########################################
############### InvSimpson ##############
##########################################

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
  stat_compare_means(method = "anova", label.y = 24) +  # Add global p-value
  stat_pvalue_manual(dunn, 
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




##########################################
############### Fisher ##############
##########################################

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
  stat_compare_means(method = "anova", label.y = 29) +  # Add global p-value
  stat_pvalue_manual(dunn, 
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






########################################################################
##################### Beta diversity #############################
########################################################################


#https://github.com/joey711/phyloseq/issues/1046

#set.seed(134) ###이거 뭐지


##################### PCoA  - bray method #############################

#calculating bray curtis distance matrix
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
  stat_conf_ellipse(colour = "transparent", #level=0.3, #level 얼마로 해야하는지 확인.
                    alpha=0.3, #색진하기 정도
                    geom = "polygon", aes(fill = group))+ 
  scale_color_brewer(palette = 'Pastel1')+
  scale_fill_brewer(palette = "Pastel1")+
  labs(title = "Bray",
       caption = paste("PERMANOVA = ", adonis2_bray$`Pr(>F)`[1],
                           ", ANOSIM = ", anosim_bray$signif))

ggsave("beta_div/bray.png", width=4, height=3, units="in", device = "png")


##################### PCoA  - unweighted #############################


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




##################### PCoA  - weighted #############################

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














########################################################################
##################### LEfSe #############################
########################################################################
library("mia")
library("lefser")
library("tidyverse")
library("SummarizedExperiment")
library("microbiomeMarker")

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


ggsave("Other/LEfSe.png", width=10, height=10, units="in", device = "png")


########################################################################
#####################  Cladogram #############################
########################################################################

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

########################################################################
#####################  Heatmap + taxonomy #############################
########################################################################

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

########################################################################
#####################  MAASLIN #############################
########################################################################

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

########################################################################
#####################  correlation  #############################
########################################################################

#****가장 많은 종 20개 찾기******
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

colnames(table_20) <-c('g__Faecalibaculum', 'g__Lactobacillus_1', 's__Lactobacillus_intestinalis', 'c__Bacilli_1', 
                       'c__Bacilli_2', 'c__Bacilli_3', 'g__Lactobacillus_2', 'g__Lactobacillus_3', 
                       'g__Bacteroides_1', 'g__Bacteroides_2', 'g__Muribaculaceae_s_Unid', 'g__Muribaculaceae', 
                       's__uncultured_bacterium_1', 's__uncultured_Bacteroidales_1', 's__uncultured_Bacteroidales_2', 's__uncultured_Bacteroidales_3', 
                       's__uncultured_bacterium_2', 'g__Helicobacter', 'f__Lachnospiraceae', 's__uncultured_Clostridiales')

#tax <- physeq20@tax_table@.Data

#***correlation matrix****
library("corrplot")

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
# 
# 
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#     stopifnot(!missing(x))
#     stopifnot(!missing(filename))
#     pdf(filename, width=width, height=height)
#     grid::grid.newpage()
#     grid::grid.draw(x$gtable)
#     dev.off()
#   }
#   save_pheatmap_pdf(xx, "test.pdf")
#   
#   
#   
# ggsave("Other/correlations.png", width=10, height=7, units="in", device = "png")


#이건 다른 ㅇ히트맵.. rela, correlation 아님
phyloseq::plot_heatmap(relaphyseq20)
#***얘는 다른 heatmap에서 활용할 수 있을 듯 ***

#*** 얘는 삼각형이 가능한 heatmap***
corrplot(correlation,
         method = 'shade',
         type = 'lower'
         )

ggsave("Other/correlation.png", width=10, height=10, units="in", device = "png")

  

########################################################################
#*****Spearman correlation****#
########################################################################

#1 

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
## adj R는 음수가 나오기도 한다....
#https://www.researchgate.net/post/Interpretation_of_negative_Adjusted_R_squared_R22
#이런 경우에는 그냥 0으로 간주한다.

ggsave(width = 3.5, height = 3.5, dpi = 300, units = "in",
       filename = paste(colnames(df_table_20[i]),".png"),
       path = "correlation/",
       device = "png"
       )
}







