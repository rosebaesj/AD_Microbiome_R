library(phyloseq)
library(tidyverse)
library("FSA")

#############Testing Line###################

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
write.table(merged_file,	file	=	"combined_otu_tax.txt",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)

#	It	seems	tedious	but	you	need	to	open	the	merged	.txt	file	in	excel	and	split into	two	files:	one	for	taxonomy	(containing	only	the	columns	OTUID	and taxonomic	info)	and	the	other	for	the	OTU	matrix	(containing	only	OTUID	and abundances	in	each	sample).	Note:	for	the	taxonomy	file,	you	need	to	use	data —>	text-to-columns	in	Excel	and	separate	on	semicolon	to	get	columns	for kingdom,	phylum,	class,	etc…	once	you	make	these	two	separate	files	in	excel, save	each	as	a	.csv
#******* OTUID, taxonomic	info -> taxonomy.csv 로 저장, taxonomy ;->,로 ㅎ변경하여 column으로 만들기*******
#******* OTUID, abundance -> OTU	matrix.csv로 저장*******


#	Step	5,	Finally,	upload	all	of	your	files	into	phyloseq	in	R!

library("ggplot2")
library("phyloseq")
library("ape")

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
physeq	=	phyloseq(OTU,	TAX,	META,	phy_tree)
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
##################### alpha diversity #############################
########################################################################

richness <- estimate_richness(rphyseq)
richness2 <- estimate_richness(rphyseq2)
richness0 <- estimate_richness(physeq)

#간단하게 전체적으로 살펴보는 방법
plot_richness(rphyseq, sortby = META$group)

##########################################
############### Chao1 ##############
##########################################


#draw barplot
ggplot(data=richness, aes(x=META$group, y=Chao1)) +
  geom_boxplot(fill=c("blue","red","green")) +
  labs(title= 'chao1', x= ' ', y= '', tag = "A") +
  geom_point()

#one way anova
one.way <- aov(Chao1 ~ META$group, data = richness)
summary(one.way)

kruskal.test(Chao1 ~ META$group, data = richness)

#pairwise.wilcox.test(chao1_data$chao1, chao1_data$Group,
#                    p.adjust.method = "BH")

dunnTest(Chao1 ~ META$group, data = richness,
         method="bonferroni")

##############  Shannon index ##############

#draw barplot
ggplot(data=richness, aes(x=META$group, y=Shannon)) +
  geom_boxplot(fill=c("blue","red","green")) +
  labs(title= 'Shannon', x= ' ', y= '', tag = "A") +
  geom_point()

#one way anova
one.way <- aov(Shannon ~ META$group, data = richness)
summary(one.way)

kruskal.test(Shannon ~ META$group, data = richness)

#pairwise.wilcox.test(chao1_data$chao1, chao1_data$Group,
#                     p.adjust.method = "BH")

#dunnTest(chao1 ~ Group, data = chao1_data,
#         method="bonferroni")

############### inverse simpson ##############

#draw barplot
ggplot(data=richness, aes(x=META$group, y=InvSimpson)) +
  geom_boxplot(fill=c("blue","red","green")) +
  labs(title= 'InvSimpson', x= ' ', y= '', tag = "A") +
  geom_point()

#one way anova
one.way <- aov(InvSimpson ~ META$group, data = richness)
summary(one.way)

kruskal.test(InvSimpson ~ META$group, data = richness)

#pairwise.wilcox.test(chao1_data$chao1, chao1_data$Group,
#                     p.adjust.method = "BH")

#dunnTest(chao1 ~ Group, data = chao1_data,
#         method="bonferroni")

########################################################################
##################### relative abundance #############################
########################################################################
#*******모름 아직 못함******
# works: from phyloseq object to relative abundance otu table
table(tax_table(rphyseq)[, "Phylum"])

rel_abund <- transform_sample_counts(rphyseq, function(x){x / sum(x)})

phylum_rel <- tax_glom(rel_abund, "Phylum")
taxa_names(ps_phylum_rel) <- tax_table(phylum_rel)[, "Phylum"]
rel_table <- as(otu_table(ps_phylum_rel), "matrix")



##아래에 있음



########################################################################
##################### PCoA  - bray method #############################
########################################################################


library(vegan)

#https://github.com/joey711/phyloseq/issues/1046

#relative abundance로 변환
relaphyseq <- transform_sample_counts(rphyseq, function(otu) {otu/sum(otu)}) #relative abundance 구하는 것으로 예상

#set.seed(134) ###이거 뭐지
##*******일단 bray로 선택. 근데 뭘로 할건지 한번 정하긴 해야함.*****#

#calculating bray curtis distance matrix
PCoA_bray <- phyloseq::distance(relaphyseq, method = "bray")

#making a data frame from the sample_data
sampledf <- data.frame(sample_data(rphyseq))
sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test
adonis2(PCoA_bray ~ META$group, data = sampledf) #0.005 
#sampledf에서 불러오는게 안돼서 그냥 META에서 부름
adonis2(PCoA_bray ~ META$group, data = sampledf2) #0.007
anosim(PCoA_bray, grouping = META$group)
##그리기...

plot <- plot_ordination(rphyseq, PCoA, color = "group")

library("ggplot2")

plot + 
  stat_ellipse(level=0.3, alpha=0.2, fill = TRUE) + #level 얼마로 해야하는지 확인.
  theme_bw()

########################################################################
##################### PCoA  - unweighted #############################
########################################################################

#calculating bray curtis distance matrix
PCoA_unifrac_d <- phyloseq::distance(relaphyseq, method = "unifrac")
PCoA_unifrac <- ordinate(relaphyseq, method ="PCoA", distance = "unifrac")

#making a data frame from the sample_data
sampledf <- data.frame(sample_data(rphyseq))
sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test
adonis2(PCoA_unifrac_d ~ META$group, data = sampledf) #0.005 
#sampledf에서 불러오는게 안돼서 그냥 META에서 부름
adonis2(PCoA_unifrac_d ~ META$group, data = sampledf2) #0.007

anosim(PCoA_unifrac_d, grouping = META$group)
##그리기...

plot <- plot_ordination(rphyseq, PCoA_unifrac, color = "group")
plot

library("ggplot2")

plot + 
  stat_ellipse(level=0.3, alpha=0.2, fill = TRUE) + #level 얼마로 해야하는지 확인.
  theme_bw()

########################################################################
##################### PCoA  - weighted #############################
########################################################################

#calculating bray curtis distance matrix
PCoA_wunifrac_d <- phyloseq::distance(relaphyseq, method = "wunifrac")
PCoA_wunifrac <- ordinate(relaphyseq, method ="PCoA", distance = "wunifrac")

#making a data frame from the sample_data
sampledf <- data.frame(sample_data(rphyseq))
sampledf2 <- data.frame(sample_data(relaphyseq))


#running adonis test
adonis2(PCoA_wunifrac_d ~ META$group, data = sampledf) #0.005 
#sampledf에서 불러오는게 안돼서 그냥 META에서 부름
adonis2(PCoA_wunifrac_d ~ META$group, data = sampledf2) #0.007

anosim(PCoA_wunifrac_d, grouping = META$group)


##그리기...

plot <- plot_ordination(rphyseq, PCoA_wunifrac, color = "group")
plot

library("ggplot2")

plot + 
  stat_ellipse(level=0.3, alpha=0.2, fill = TRUE) + #level 얼마로 해야하는지 확인.
  theme_bw()

########################################################################
##################### LEfSe #############################
########################################################################
library("mia")
library("lefser")
library("tidyverse")
library("SummarizedExperiment")

relasum <- makeTreeSummarizedExperimentFromPhyloseq(relaphyseq)
unique(relasum$SampleType)

SummarizedExperiment(relaphyseq)
lefser(relasum, groupCol = "group")
#*****lefse에서 2개 군만 선별해서 돌리ㄹ는 법을 찾아야함-> 해결 못했음****###
#*****lefse와 maaslin에서도 동일 data set으로 돌아가려면 rarefaction을 R로 넘어오기 전에 해야할듯??****###
#Summarized Experiment로 아예 그룹을 빼고 physeq를 만들어봐야할 긋?

library("microbiomeMarker")
lefse <- run_lefse(relaphyseq,"group")

lefse_data <-data.frame(marker_table(lefse))
plot_ef_bar(lefse)


##*****cladogram 오류 햐결 안됨...****##
library(ggtree)
####https://yulab-smu.top/treedata-book/chapter4.html

fan.angle ()
ggtree(relaphyseq, layout = 'circular', 
       branch.length = 3, ladderize) #우리가 원하는게 circular 인건 맞는 듯


data(lefse)
lefse_small <- phyloseq::subset_taxa(
  lefse,
  Phylum %in% c("Firmicutes")
)
plot_cladogram(lefse, color = c(AD ="red", AP = "blue", CON = "green"), only_marker = TRUE,  clade_label_level = 4)


plot_cladogram(mm_lefse, color = c(Healthy = "darkgreen", Tumor = "red")) +
  theme(plot.margin = margin(0, 0, 0, 0))


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

#****논문에서 쓴것은 pheatmap package***#
library("pheatmap")
relaphyseq@otu_table@.Data
pheatmap(df_num_scale, main = "pheatmap default")

b <- run_test_two_groups(relaphyseq, 'group')

########################################################################
#####################  MAASLIN #############################
########################################################################

#
#library("devtools") #github에서 다운받으려면 필요한 것 
#LEfSe identifies those data features that are distinct between a pair of metadatums (e.g. differences between two sampling sites, two clinical outcomes, two biochemical markers, two modalities, etc.).  MaAsLin extends the functionality of LEfSe to identify associations between data features and multiple metadata factors, which can be discrete and/or continuous and can include time series data
library("yingtools2")
library("tidyverse")
library("Maaslin2")

#https://github.com/biobakery/maaslin#markdown-header-input-files
#*****위의 형식을 보고 input file 만들기****

maaslin <- Maaslin2(input_data = 'MaAsLin/feature_table_L7.txt',# sample과 feature만 있는 것, 이미 rarefaction 되었음
                    input_metadata = 'MaAsLin/metadata.txt', #group 정보 등
                    analysis_method = 'LM', #기본 값이 LM인데 다른 방법들의 차이를 확인할수가 없음
                    output = 'MaAsLin/maaslin', #output folder 이름
                    reference = c('Group;AD;AP;CON')) #다변수인 경우 나열할 것.





########################################################################
#####################  correlation  #############################
########################################################################

#****가장 많은 종 20개 찾기******
topN <- 20
most_abundant_taxa <- sort(taxa_sums(relaphyseq), TRUE)[1:topN]
print(most_abundant_taxa)
#physeq20 <- prune_taxa(names(most_abundant_taxa), physeq)
relaphyseq20 <- prune_taxa(names(most_abundant_taxa), relaphyseq) #얘는 안됨....
#뭔가 이 prune 명령어가 subgroup 가능하게 하는 것 같음

table_20 <- t(as.matrix(relaphyseq20@otu_table@.Data))
colnames(table_20) <-c('g__Faecalibaculum', 'g__Lactobacillus_1', 's__Lactobacillus_intestinalis', 'c__Bacilli_1', 
                       'c__Bacilli_2', 'c__Bacilli_3', 'g__Lactobacillus_2', 'g__Lactobacillus_3', 
                       'g__Bacteroides_1', 'g__Bacteroides_2', 'Unid_g__Muribaculaceae', 'g__Muribaculaceae', 
                       's__uncultured_bacterium_1', 's__uncultured_Bacteroidales_1', 's__uncultured_Bacteroidales_2', 's__uncultured_Bacteroidales_3', 
                       's__uncultured_bacterium_2', 's__uncultured_bacterium_3', 'g__Helicobacter', 's__uncultured_Clostridiales')

#tax <- physeq20@tax_table@.Data

#***correlation matrix****

correlation <- cor(table_20)
pheatmap(correlation, scale = "none", 
         drop_levels = TRUE, #아무 영향 없음
         angle_col = "90",
)

########################################################################
#*****Spearman correlation****#
########################################################################

spearman_1 <-data.frame(richness$Shannon)
spearman_1 <-cbind(spearman_1, table_20[,3])
t_spearman_1 <- data.frame(t(spearman_1))
colnames(spearman_1) <- c('Shannon', 's__Lactobacillus_intestinalis')

fit <- lm(Shannon ~ s__Lactobacillus_intestinalis, data = spearman_1)

dev.new(width=5, height=5, unit="px")


p <- ggplot(spearman_1,aes(s__Lactobacillus_intestinalis,Shannon))+
  geom_point(col = 'pink')+
  #  stat_summary(fun.data=mean_cl_normal) +
  geom_smooth(method='lm', col = 'pink', fill = 'pink')+
  xlab('Relative Abundance')+
  labs(title = "s__Lactobacillus_intestinalis", 
       subtitle  = paste(" R2 = ",signif(summary(fit)$adj.r.squared, 5),
                         #                   "\n", "Intercept =",signif(fit$coef[[1]],5 ),
                         #                   "Slope =",signif(fit$coef[[2]], 5), "\n",
                         ", P =",signif(summary(fit)$coef[2,4], 5)))+
  theme_bw()+
  theme(axis.line = element_line(size=1),
        axis.ticks = element_line(size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5)
  )
#수기로 correlation_graphs 폴더를 생성
ggsave(plot = p, width = 3, height = 3, dpi = 300, 
       filename = "s__Lactobacillus_intestinalis.pdf", 
       path = 'phyloseq/correlation_graphs')

########################################################################
#*****Spearman correlation****#
########################################################################







