#library에 없는게 있으면 보통 install.packages("")로 설치하면 되고 안되는 경우엔 검색하면 바로 나옴
library(ape)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(plyr)
library(dplyr)
library(tidyverse)
library(FSA)
library(RColorBrewer)
library(vegan)
library(rstatix)
library(metagMisc)
library(MicrobeR)
library(mia)
library(lefser)
library(SummarizedExperiment)
library(microbiomeMarker)
library(pheatmap)
library(berryFunctions)
library(matrixTests)
library(rstatix)
library(dunn.test)


getwd()
setwd("All_AD")
#이건 자기가 설정한 디렉토리 잘 찾아서 진행하기


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###########################IMPORT DATA##################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#https://nephele.niaid.nih.gov/user_guide_phyloseq_tutorial/
#https://github.com/joey711/phyloseq/issues/821
#아래 깃허브 주소에 pdf 올려주신 분이 있는데 그거 따라하면 됨.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#####+ import from qiime ############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#******* 수기로 #OTU ID -> OTUID로 변경해야함, 또한 taxonomy라고 써져있는 부분 지워야함.*******
#	read	in	OTU	table	
##otu	<- read.table(file	=	"phyloseq/DADA2_table.txt",	header	=	TRUE)
otu<- read.table(file="qiime/otu_table.txt",header = TRUE)

#******* 역시 수기로 OTUID, not Feature ID*******
#	read	in	taxonomy	table #taxonomy 에 Kingdom... 이런 식으로 일일히 column name 작성해야함
##tax	<-	read.table(file	=	"phyloseq/taxonomy.tsv",	sep	=	'\t',	header	=	TRUE)
tax<-read.table(file="qiime/taxonomy.tsv",sep='\t',header=TRUE)

#	merge	files	
merged_file<-merge(otu,tax,by.x=c("OTUID"),by.y=c("OTUID"))

#	note:	number	of	rows	should	equal	your	shortest	Sile	length,	drops	taxonomy	for	OTUs	that	don’t	exist	in	your	OTU	table
#	output	merged	.txt	Pile
##write.table(merged_file,	file	=	"phyloseq/combined_otu_tax.tsv",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)
write.table(merged_file,file="qiime/combined_otu_tax.tsv",sep='\t',col.names=TRUE,row.names=FALSE)

#	It	seems	tedious	but	you	need	to	open	the	merged	.txt	file	in	excel	and	split into	two	files:	one	for	taxonomy	(containing	only	the	columns	OTUID	and taxonomic	info)	and	the	other	for	the	OTU	matrix	(containing	only	OTUID	and abundances	in	each	sample).	Note:	for	the	taxonomy	file,	you	need	to	use	data —>	text-to-columns	in	Excel	and	separate	on	semicolon	to	get	columns	for kingdom,	phylum,	class,	etc…	once	you	make	these	two	separate	files	in	excel, save	each	as	a	.csv
#******* OTUID, taxonomic	info -> taxonomy.csv 로 저장, taxonomy ;->,로 ㅎ변경하여 column으로 만들기*******
#******* OTUID, abundance -> OTU	matrix.csv로 저장*******


#	Step	5,	Finally,	upload	all	of	your	files	into	phyloseq	in	R!


#	read	in	otu	table
##otu_table	=	read.csv("phyloseq/otu_matrix.csv",	sep=",",	row.names=1)
otu_table=read.csv("qiime/otu_matrix.csv",sep=",",row.names=1)
otu_table=as.matrix(otu_table)
# ##수기로 만들어야 함
# otu_g_table=read.csv("qiime/otu_g_matrix.csv",sep=",",row.names=1)
# otu_g_table = as.matrix(otu_g_table)

#	read	in	taxonomy
#	seperated	by	kingdom	phylum	class	order	family	genus	species
# taxonomy 에 Kingdom... 이런 식으로 일일히 column name 작성해야함
##taxonomy	=	read.csv("phyloseq/taxonomy.csv",	sep=",",	row.names=1)
taxonomy=read.csv("qiime/taxonomy.csv",sep=",",row.names=1)
taxonomy=as.matrix(taxonomy)

#	read	in	metadata	
#	variables	=	???
##metadata	=	read.table("phyloseq/metadata.tsv",	row.names=1)
metadata=read.table("qiime/metadata.tsv",row.names=1)
colnames(metadata)<-metadata[1,]
metadata <- metadata[-1,]

# metadata_g = read.table("qiime/metadata_g.tsv",row.names=1)
# colnames(metadata_g)<-metadata_g[1,]
# metadata_g <- metadata_g[-1,]


#	read	in	tree
##phy_tree	=	read_tree("phyloseq/tree.nwk")
phy_tree	=	read_tree("qiime/tree.nwk")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###### + make phyloseq objects #######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

OTU=otu_table(otu_table,taxa_are_rows=TRUE)
#OTUg=otu_table(otu_g_table, taxa_are_rows=TRUE)
TAX=tax_table(taxonomy)
META=sample_data(metadata)
#METAg=sample_data(metadata_g)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

#	make	sure	files	have	the	same	sample	names	
sample_names(OTU)
sample_names(META)

#	merge	into	one	phyloseq	object
#physeq <- phyloseq(OTU, TAX, META, phy_tree) ##원칙상 순서는 이게 아닌데 걍 되는듯
rphyseq <- phyloseq(OTU, phy_tree, TAX, META)
#rphyseq_g <- phyloseq(OTUg, phy_tree, TAX, METAg)
#	Now,	continue	to	analysis	in	phyloseq!

#Ready for Microbiome analysis



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + rarefy feature table #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#*******Rarefy를 이미 해서 qiime에서 가져오는 것으로 변경하였음******
#rarefy란 같은 수의 read로 맞추는 것 
#기존의 qiime에서 가장 작은 수의 sample이 rarefaction(?) graph에서 평평하게 유지ㅎ되는지 확인해야함
#우리 샘플의 경우 13410개가 최소였음*******<<과거형. 다시 checkup 해야함
##rphyseq <- rarefy_even_depth(physeq) #기본값이 최솟값으로 맞추는 것
##rphyseq2 <- rarefy_even_depth(physeq, sample.size = 13410) #기본값이 최솟값으로 맞추는 것



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + relative abundance #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#relative abundance로 변환
##relaphyseq <- transform_sample_counts(rphyseq, function(otu) {otu/sum(otu)}) #relative abundance 구하는 것으로 예상
relaphyseq <- transform_sample_counts(rphyseq, function(otu) {otu/sum(otu)})
# relaphyseq_g <- transform_sample_counts(rphyseq_g, function(otu) {otu/sum(otu)})




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + taxonomic agglomeration #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# 얘는 중복되는 AVS를 OTU로 변환하는 것, 완전 별도....
rela_tax_s7 <- tax_glom(relaphyseq, taxrank = "Species")
rela_tax_g6 <- tax_glom(relaphyseq, taxrank = "Genus")
rela_tax_f5 <- tax_glom(relaphyseq, taxrank = "Family")
rela_tax_o4 <- tax_glom(relaphyseq, taxrank = "Order")
rela_tax_c3 <- tax_glom(relaphyseq, taxrank = "Class")
rela_tax_p2 <- tax_glom(relaphyseq, taxrank = "Phylum")
rela_tax_k1 <- tax_glom(relaphyseq, taxrank = "Kingdom")



#### 다른 시도를 해보았당 처참히 실패했당!
# Phylum_table <- tax_glom(relaphyseq, taxrank = "Phylum")
# o<- data.frame(Phylum_table@otu_table)
# o<- rownames_to_column(o)
# p<- data.frame(Phylum_table@tax_table[,2])
# p<- rownames_to_column(p)
# phylum <- left_join(o, p, by = c("rowname"="rowname"))
# phylum <- phylum[,-1]
# phylum[nrow(phylum)+1,] <- colnames(phylum)
# phylum <- column_to_rownames(phylum, "Phylum")
# 
# ggplot(phylum, x="Phylum")+
#   geom_bar(position="stack", stat="identity")
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# data <- data.frame(specie,condition,value)


# 얘는 가까운 taxa를 합치는 것. OTU와 같은 개념 
# 이건 많은 수가 taxonomic assignment가 안된 경우에 활용해볼 수 있다.
# taxa 수가 많은 경우 많이 느리다. 그리고 tree가 반드시 있어야 한다.
rela_tip <- tip_glom(relaphyseq)
# 애는 훨씬 줄었음 57개가 됐음



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### RELATIVE ABUNDANCE (for STACKED BAR GRAPH) #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##### + out put, 개체 별로 꺼내야 한다 ######

# Extract abundance matrix from the phyloseq object

out_rela <- function(phyloseq, filepath){
  o = as(otu_table(phyloseq), "matrix")
  o = as.data.frame(o)
  t = as(tax_table(phyloseq), "matrix")
  t = as.data.frame(t)
  rownames(o)
  rownames(t)
  
  j = merge(t, o, by = 0)
  write.table(j, file = filepath, sep='\t',col.names=TRUE,row.names=FALSE)
  return(j)
}

out_asv <- out_rela(relaphyseq, 'tables/asv.tsv')
out_s7 <- out_rela(rela_tax_s7, 'tables/Species.tsv')
out_g6 <- out_rela(rela_tax_g6, "tables/Genus.tsv")
out_f5 <- out_rela(rela_tax_f5, "tables/Family.tsv")
out_o4 <- out_rela(rela_tax_o4, "tables/Order.tsv")
out_c3 <- out_rela(rela_tax_c3, "tables/Class.tsv")
out_p2 <- out_rela(rela_tax_p2, "tables/Phylum.tsv")
#out_k1 <- out_rela(rela_tax_k1, "tables/Kingdom.tsv")



out_asv <- rbind(c(NA, NA, NA, NA, NA, NA, NA, NA, 
                  sample_data(relaphyseq)$group), out_asv)
#rownames(out_all) <- out_all[,1] 겹치는게 있어서 안됨 ㅠㅠ
rownames(out_asv)[1] <- "group"


tax_asv <- out_asv[,1:8]
otu_asv <- out_asv[,9:ncol(out_asv)]
t_asv <- data.frame(t(otu_asv))


#####+ taxonomy ######

out_tax <- rbind(out_p2, out_c3, out_o4, out_f5, out_g6, out_s7) 
#Kingdom이 위에 있어야 구별이 잘됨
#얘네는 그룹 정보가 빠져있음. 

out_tax <- rbind(c(NA, NA, NA, NA, NA, NA, NA, NA, 
                   sample_data(relaphyseq)$group), out_tax)
rownames(out_tax)[1] <- "group"

tax_tax <- out_tax[,1:8]
otu_tax <- out_tax[,9:ncol(out_tax)]
t_tax <- data.frame(t(otu_tax))



##### + statistics asv ######

kwp <- NULL
dunnt <- NULL

for (i in 1:ncol(t_asv)){
  
  #Kruskal test
  kw <- kruskal.test(t_asv[,i] ~ group, t_asv)
  kw$p.value
  kwp <- rbind(kwp, kw$p.value)
  
  #Dunn Test, between groups
  gg <- data.frame(cbind(t_asv[,i], t_asv$group))
  colnames(gg) <- c('a', 'group')
  dunn <- dunn_test (data=gg, a ~ group, p.adjust.method = "bonferroni")
  dunn$p.adj
  dunnt <- rbind(dunnt, t(dunn$p.adj))
  
}
#오래걸림


colnames(kwp) <- c('Kruskal_test')
View(dunn) #여기에서 나오는 그룹 순서대로 적으면 됨
colnames(dunnt) <- c('Dunn_AD-AP', 'Dunn_AD-CON', 'Dunn_AP-CON')
stat_out_asv <- cbind(out_asv, kwp, dunnt)

##### + statistics tax ######

kwp <- NULL
dunnt <- NULL

for (i in 1:ncol(t_tax)){
  
  #Kruskal test
  kw <- kruskal.test(t_tax[,i] ~ group, t_tax)
  kw$p.value
  kwp <- rbind(kwp, kw$p.value)
  
  #Dunn Test, between groups
  gg <- data.frame(cbind(t_tax[,i], t_tax$group))
  colnames(gg) <- c('a', 'group')
  dunn <- dunn_test (data=gg, a ~ group, p.adjust.method = "bonferroni")
  dunn$p.adj
  dunnt <- rbind(dunnt, t(dunn$p.adj))
  
}
#오래걸림 
#Kingdom이 포함되면 너무 많으면 오류 나는 것 같음


colnames(kwp) <- c('Kruskal_test')
View(dunn) #여기에서 나오는 그룹 순서대로 적으면 됨
colnames(dunnt) <- c('Dunn_AD-AP', 'Dunn_AD-CON', 'Dunn_AP-CON')
stat_out_tax <- cbind(out_tax, kwp, dunnt)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### ALPHA DIVERSITY #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# transposed compared to 
richness <- data.frame(estimate_richness(rphyseq))
richness <- rbind(metadata$group, richness)
#richness$group <- factor (richness$group, levels = c("CON", "AD", "AP"))


#간단하게 전체적으로 살펴보는 방법
##plot_richness(rphyseq, sortby = META$group)
plot_richness(physeq, sortby = META$group)

##ppplot에서 군별 paired 분석에 대한 pvalue 표기하려면 pair를 정해줘야함.
pair <- list (c("AD", "AP"), c("AD", "CON"), c("AP", "CON"))

t_richness <- as.data.frame(t(richness))

kwp <- NULL
dunnt <- NULL

for (i in 1:ncol(richness)){
  
  #Kruskal test
  kw <- kruskal.test(richness[,i] ~ group, richness)
  kw$p.value
  kwp <- rbind(kwp, kw$p.value)
  
  #Dunn Test, between groups
  gg <- data.frame(cbind(richness[,i], richness$group))
  colnames(gg) <- c('a', 'group')
  dunn <- dunn_test (data=gg, a ~ group, p.adjust.method = "bonferroni")
  dunnt <- rbind(dunnt, t(dunn$p.adj))
  
}

colnames(kwp) <- c('Kruskal_test')
View(dunn) #여기에서 나오는 그룹 순서대로 적으면 됨
colnames(dunnt) <- c('Dunn_AD-AP', 'Dunn_AD-CON', 'Dunn_AP-CON')
t_richness <- cbind(t_richness, kwp, dunnt)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### BETA DIVERSITY #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



















#https://github.com/joey711/phyloseq/issues/1046

#set.seed(134) ###이거 뭐지

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### + PCoA  - bray method #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
##################### + correlation  #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#****가장 많은 종 20개 찾기******
topN <- 20


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### ++ ASV => relaphyseq  #############################
most_abd <- sort(taxa_sums(relaphyseq), TRUE)[1:topN]
print(most_abd)
#physeq20 <- prune_taxa(names(most_abd), physeq)
relaphyseq20 <- prune_taxa(names(most_abd), relaphyseq) 


##################### ++ tax_glom => rela_tax  #############################
most_abd <- sort(taxa_sums(rela_tax), TRUE)[1:topN]
print(most_abd)
#physeq20 <- prune_taxa(names(most_abd), physeq)
relaphyseq20 <- prune_taxa(names(most_abd), rela_tax) 

most_abd <- sort(taxa_sums(rela_tax_g), TRUE)[1:topN]
print(most_abd)
#physeq20 <- prune_taxa(names(most_abd), physeq)
relaphyseq20 <- prune_taxa(names(most_abd), rela_tax_g) 

##################### ++ tip_glom => rela_tip  #############################
most_abd <- sort(taxa_sums(rela_tip), TRUE)[1:topN]
print(most_abd)
#physeq20 <- prune_taxa(names(most_abd), physeq)
relaphyseq20 <- prune_taxa(names(most_abd), rela_tip) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

otu_20 <- t(as.matrix(relaphyseq20@otu_table@.Data))
tax_20 <- data.frame(relaphyseq20@tax_table)
View(tax_20)

df_otu_20 <- data.frame(otu_20)
colnames(df_otu_20) ### 에서 나온 값들 찾아서 이름 넣기

## colnames(otu_20) <-c('g__Faecalibaculum', 'g__Lactobacillus_1', 's__Lactobacillus_intestinalis', 'c__Bacilli_1', 
##                       'c__Bacilli_2', 'c__Bacilli_3', 'g__Lactobacillus_2', 'g__Lactobacillus_3', 
##                       'g__Bacteroides_1', 'g__Bacteroides_2', 'g__Muribaculaceae_s_Unid', 'g__Muribaculaceae', 
##                       's__uncultured_bacterium_1', 's__uncultured_Bacteroidales_1', 's__uncultured_Bacteroidales_2', 's__uncultured_Bacteroidales_3', 
##                       's__uncultured_bacterium_2', 'g__Helicobacter', 'f__Lachnospiraceae', 's__uncultured_Clostridiales')

# colnames(df_otu_20) <-c('1g__Lactobacillus', '2s__Lactobacillus_intestinalis', '3g__Helicobacter', '4s__uncultured_Bacteroidales', 
#                        '5s__uncultured_Bacteroidales', '6s__uncultured_Bacteroidales', '7s__uncultured_bacterium', '8s__uncultured_bacterium', 
#                        '9g__Muribaculaceae', '10s__unidentified', '11g__Bacteroides', '12g__Bacteroides', 
#                        '13s__uncultured_bacterium', '14s__uncultured_bacterium', '15g__Lactobacillus', '16g__Lactobacillus', 
#                        '17c__Bacilli', '18c__Bacilli', '19c__Bacilli', '20s__uncultured_Clostridiales')
# ##이름 겹치면 그래프 덮어쓰기 되니까 주의

colnames(df_otu_20) <- tax_20$Genus

#tax <- physeq20@tax_table@.Data

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








