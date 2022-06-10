#BEHAVIOR

# library("phyloseq")
# library("ggplot2") 
# library("microbiomeMarker")
# library(ape)
# library(ggpubr)
# library(plyr)
# library(dplyr)
# library(tidyverse)
# library(FSA)
# library(RColorBrewer)
# library(vegan)
# library(rstatix)
# library(metagMisc)
# library(MicrobeR)
# library(mia)
#library(SummarizedExperiment)
# library(pheatmap)
# library(berryFunctions)
# library(matrixTests)
# library(rstatix)
# library(philr)
library(stringr)

setwd("All_AD")
getwd()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######################### IMPORT DATA ###########################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#+ data set은 다음과 같은 사항을 잘 지켜서 만들 것
#+ #+ ctrl+f blank 띄어쓰기 를 찾아서 모두 change into "_"
#+ single row and column matrix로 만든다
#+ no merged cell
#+ =CONCAT(A2, A3) merges cell texts.
#+ no empty cells....  비어도 될수도 있긴 한데 없는 것을 가정해서 진행함...
#+ 간혹 오류가 있을 수 있으므로 원하는 값만 있는 matrix를 복사, new sheet에 paste 해서
#+ save as .txt

behavior <- read.table(file="input_behavior/AD-acu-behavior_edit.txt", 
                      sep='\t', header = TRUE)
rownames(behavior) <- behavior[,1]

#불필요한 데이터 지우기
behavior <- behavior [,c(-1, -2)]

# if features are column names and sample names are rows, transpose them
behavior <- data.frame(t(behavior))

metadata=read.table("input_behavior/metadata.tsv", header=TRUE)
group <- metadata$group

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### KDMS STATISTICS #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#*** Full_table.R line 191 kdms function needed 미리 불러야함 ***#

stat_behav <- kdms (behavior, group)
#시간 걸림.
stat_behavior <- cbind (rownames(behavior), stat_behav, behavior)

write.table(stat_behavior, file ='output_behavior/tables/stat_behavior.tsv' , 
            sep='\t', col.names=TRUE,row.names=FALSE)







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### CORRELATION ANALYSIS #############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#+ correlation analysis requires features in column names 
#+ and sample names in row names
#+ which means the opposite of micorobiome feature tables

#+ you can add up any type of transposed data as matrix and simply run cor()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + import data for correlation analysis #####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#+ stat output 들에서 stat 값만 제외하고 복사해서 input_correlation으로 저장!
#+ >이 과정은 추후에 각 process 에 자동화해서 넣는 것이 좋을 듯
#+ feature가 rowname에 , sample이 colname에 있는 것들을 


####### + microbiome features ##########


#read it from results if you have erased "out_asv" or any other similar objects from environment
#out_asv <- read.table("tables/rela_asv.tsv", header = TRUE)

cor_asv <- out_asv[,c(20:ncol(out_asv))]
rownames(cor_asv) <- out_asv[,1]
cor_2p <- out_2p[,c(20:ncol(out_2p))]
rownames(cor_2p) <- out_2p[,3]
cor_3c <- out_3c[,c(20:ncol(out_3c))]
rownames(cor_3c) <- out_3c[,4]
cor_4o <- out_4o[,c(20:ncol(out_4o))]
rownames(cor_4o) <- out_4o[,5]
cor_5f <- out_5f[,c(20:ncol(out_5f))]
rownames(cor_5f) <- out_5f[,6]

############ +++ Species Genus 하나의 행으로 합치기 for species ####
gf_6g <- str_c(out_6g$Genus, "_in_", out_6g$Family)
cor_6g <- out_6g[,c(20:ncol(out_6g))]
rownames(cor_6g) <- gf_6g

# ############ @+++ Species Genus 하나의 행으로 합치기 for species ####
# sg_7s <- str_c(out_7s$Species, "_in_", out_7s$Genus)
# cor_7s <- out_7s[,c(20:ncol(out_7s))]
# rownames(cor_7s) <- sg_7s
### @ 이래도 같아서 해결이 아됨

####### + behavior features ##########

cor_behavior <- behavior

###### + diversity features ##########
# read it from results if you have erased "alpha_div", "stat_beta_div" objects from environment
# alpha_div <- read.table("tables/stat_alpha_div.tsv", header = TRUE)
cor_alpha <- alpha_div[,c(12:ncol(alpha_div))]
rownames(cor_alpha) <- alpha_div[,1]





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + make feature matrix #####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cor <- rbind(cor_6g, cor_behavior) #whatever you want




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### run correlation analysis #####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

n_cor <- matrix(as.numeric(unlist(cor)), nrow = nrow(cor))
colnames(n_cor)<- colnames(cor)
rownames(n_cor)<- rownames(cor)
t_cor <- t(n_cor)
pearson <- cor(t_cor,  method = "pearson")

bug_behav_cor <- pearson[1:nrow(cor_6g), (nrow(cor_6g)+1):nrow(pearson)]

write.table(bug_behav_cor, file = "tables/bug_behav_cor.tsv", sep='\t',
            col.names=TRUE,row.names=FALSE)



##
