library(phyloseq)
library(tidyverse)
library("FSA")
library(ggplot2)
library(reshape2)
library(phyloseq)
library(ape)





########################################################################
###########################IMPORT DATA##################################
########################################################################


#******* 수기로 #OTU ID -> OTUID로 변경해야함, 또한 taxonomy라고 써져있는 부분 지워야함.*******
#	read	in	OTU	table	
otu	<-	read.table(file	=	"DADA2_table.txt",	header	=	TRUE)


#******* 역시 수기로 OTUID, not Feature ID*******
#	read	in	taxonomy	table #taxonomy 에 Kingdom... 이런 식으로 일일히 column name 작성해야함
tax	<-	read.table(file	=	"taxonomy.txt",	sep	=	'\t',	header	=	TRUE)


#	merge	files	
merged_file	<-	merge(otu,	tax,	by.x	=	c("OTUID"),	by.y=c("OTUID"))


#	note:	number	of	rows	should	equal	your	shortest	Sile	length,	drops	taxonomy	for	OTUs	that	don’t	exist	in	your	OTU	table
#	output	merged	.txt	Pile
write.table(merged_file,	file	=	"combined_otu_tax.txt",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)

#	It	seems	tedious	but	you	need	to	open	the	merged	.txt	file	in	excel	and	split into	two	files:	one	for	taxonomy	(containing	only	the	columns	OTUID	and taxonomic	info)	and	the	other	for	the	OTU	matrix	(containing	only	OTUID	and abundances	in	each	sample).	Note:	for	the	taxonomy	file,	you	need	to	use	data —>	text-to-columns	in	Excel	and	separate	on	semicolon	to	get	columns	for kingdom,	phylum,	class,	etc…	once	you	make	these	two	separate	files	in	excel, save	each	as	a	.csv
#******* OTUID, taxonomic	info -> taxonomy.csv 로 저장, taxonomy ;->,로 ㅎ변경하여 column으로 만들기*******
#******* OTUID, abundance -> OTU	matrix.csv로 저장*******


#	Step	5,	Finally,	upload	all	of	your	files	into	phyloseq	in	R!

#	read	in	otu	table
otu_table	=	read.csv("otu_matrix.csv",	sep=",",	row.names=1)
otu_table	=	as.matrix(otu_table)


#	read	in	taxonomy
#	seperated	by	kingdom	phylum	class	order	family	genus	species
# taxonomy 에 Kingdom... 이런 식으로 일일히 column name 작성해야함
taxonomy	=	read.csv("taxonomy.csv",	sep=",",	row.names=1)
taxonomy	=	as.matrix(taxonomy)


#	read	in	metadata	
#	variables	=	???
metadata	=	read.table("metadata.tsv",	row.names=1)
colnames(metadata)<-metadata[1,]
metadata <- metadata[-1,]


#	read	in	tree
phy_tree	=	read_tree("tree.nwk")

#	import	as	phyloseq	objects
OTU	=	otu_table(otu_table,	taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META	=	sample_data(metadata)
#	(tree	was	already	imported	as	a	phyloseq	object)


#	merge	into	one	phyloseq	object
physeq	=	phyloseq(OTU,	TAX,	META,	phy_tree)


#	Now,	continue	to	analysis	in	phyloseq!




########Welcome to Github########



########################################################################
##################### Rarefy feature table #############################
########################################################################

#rarefy란 같은 수의 read로 맞추는 것 
#기존의 qiime에서 가장 작은 수의 sample이 rarefaction(?) graph에서 평평하게 유지ㅎ되는지 확인해야함
#*******우리 샘플의 경우 13410개가 최소였음*******
rphyseq <- rarefy_even_depth(physeq) #기본값이 최솟값으로 맞추는 것
rphyseq2 <- rarefy_even_depth(physeq, sample.size = 13410) #기본값이 최솟값으로 맞추는 것



relaphyseq <- transform_sample_counts(rphyseq, function(otu) {otu/sum(otu)})








### All samples manual color palette
plot_bar(relaphyseq, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#8E14F3", "#4B66FF", "#3CF3FF", "#1FEC22", "#A1FE13", "#FBC417", "#EC6F1E", "#F03309", "#FE6666", "#FA217A", "#FC7CF0")) +
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




###All samples color palette
plot_bar(gphyseq, fill = "Phylum") +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +
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









##Group별 stacked bar plot


###기본 셋팅은 앞의 것과 같음. 다만 아직은 직접 손 보는 형태. 분명... 뭔가 있을텐데.
gotu <- read.table(file = "0302_table.txt", header = TRUE)
gtax	<-	read.table(file	=	"gtaxonomy.txt",	sep	=	'\t',	header	=	TRUE)
gmerged_file	<-	merge(gotu,	gtax,	by.x	=	c("OTUID"),	by.y=c("OTUID"))
write.table(merged_file,	file	=	"combined_gotu_tax.txt",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)

gotu_table = read.csv("gotu_matrix_group.csv", sep=",", row.names=1)
gotu_table	=	as.matrix(gotu_table)

gtaxonomy	=	read.csv("gtaxonomy.csv",	sep=",",	row.names=1)
gtaxonomy	=	as.matrix(gtaxonomy)

gmetadata	=	read.table("gmetadata.tsv",	row.names=1)
colnames(gmetadata)<-gmetadata[1,]
gmetadata <- gmetadata[-1,]

phy_tree	=	read_tree("tree.nwk")

GOTU	=	otu_table(gotu_table,	taxa_are_rows	=	TRUE)
GTAX	=	tax_table(gtaxonomy)
GMETA	=	sample_data(gmetadata)

gphyseq = phyloseq(GOTU, GTAX, GMETA, phy_tree)

grelaphyseq <- transform_sample_counts(gphyseq, function(otu) {otu/sum(otu)})



plot_bar(grelaphyseq, fill = "Class") +
  geom_bar(aes(fill = Class), stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +
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
















########Metadata를 이용해서 rarefaction된 정보들 sum하기.

otu	<-	read.table(file	=	"0302_table.txt",	header	=	TRUE)
tax	<-	read.table(file	=	"taxonomy.txt",	sep	=	'\t',	header	=	TRUE)
merged_file	<-	merge(otu,	tax,	by.x	=	c("OTUID"),	by.y=c("OTUID"))
write.table(merged_file,	file	=	"combined_otu_tax.txt",	sep	=	'\t',	col.names	= TRUE,	row.names	=	FALSE)


otu_table	=	read.csv("otu_matrix.csv",	sep=",",	row.names=1)
otu_table	=	as.matrix(otu_table)


taxonomy	=	read.csv("taxonomy.csv",	sep=",",	row.names=1)
taxonomy	=	as.matrix(taxonomy)


metadata	=	read.table("metadata.tsv",	row.names=1)
colnames(metadata)<-metadata[1,]
metadata <- metadata[-1,]


phy_tree	=	read_tree("tree.nwk")

Metadata<-
  etadata %>% 
  left_join(gotu)


OTU	=	otu_table(otu_table,	taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy)
META	=	sample_data(metadata)


relaphyseq <- transform_sample_counts(rphyseq, function(otu) {otu/sum(otu)})







