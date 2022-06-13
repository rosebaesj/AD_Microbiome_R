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
library(vegan)
# library(rstatix)
# library(metagMisc)
# library(MicrobeR)
# library(mia)
#library(SummarizedExperiment)
# library(pheatmap)
# library(berryFunctions)
# library(matrixTests)
library(rstatix)
# library(philr)
library(stringr)
library(psych)


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

#function to make taxonomic to needed format
cor_features <- function(data, sample_rows, feature_names){
  c <- data[,c(sample_rows)]
  c <- matrix(as.numeric(unlist(c)), nrow = nrow(c))
  rownames(c) <- feature_names
  colnames(c) <- colnames(data[,c(sample_rows)])
  c <- t(c)
  return(c)
}

#read it from results if you have erased "out_asv" or any other similar objects from environment
# out_asv <- read.table("tables/rela_asv.tsv", header = TRUE)
# out_2p <- read.table("tables/rela_Phylum.tsv", header = TRUE)
# out_3c <- read.table("tables/rela_Class.tsv", header = TRUE)
# out_4o <- read.table("tables/rela_Order.tsv", header = TRUE)
# out_5f <- read.table("tables/rela_Family.tsv", header = TRUE)
# out_6g <- read.table("tables/rela_Genus.tsv", header = TRUE)
# out_7s <- read.table("tables/rela_Species.tsv", header = TRUE)


cor_asv <- cor_features(out_asv, 20:ncol(out_asv), out_asv[,1])

cor_2p <- cor_features(out_2p, 20:ncol(out_2p), out_2p[,3])
cor_3c <- cor_features(out_3c, 20:ncol(out_3c), out_3c[,4])
cor_4o <- cor_features(out_4o, 20:ncol(out_4o), out_4o[,5])
cor_5f <- cor_features(out_5f, 20:ncol(out_5f), out_5f[,6])

############ +++ Species Genus 하나의 행으로 합치기 for species ####
gf_6g <- str_c(out_6g$Genus, "_in_", out_6g$Family)
cor_6g <- cor_features(out_6g, 20:ncol(out_6g), gf_6g)


# ############ @+++ Species Genus 하나의 행으로 합치기 for species ####
# sg_7s <- str_c(out_7s$Species, "_in_", out_7s$Genus)
# cor_7s <- out_7s[,c(20:ncol(out_7s))]
# rownames(cor_7s) <- sg_7s
### @ 이래도 같아서 해결이 아됨

####### + behavior features ##########

cor_behavior <- cor_features(behavior, 1:ncol(behavior), rownames(behavior))


###### + diversity features ##########
# read it from results if you have erased "alpha_div", "stat_beta_div" objects from environment
# alpha_div <- read.table("tables/stat_alpha_div.tsv", header = TRUE)
cor_alpha <- cor_features(alpha_div, 12:ncol(alpha_div), alpha_div[,1])





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### + make feature matrix #####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# cor <- cbind(cor_6g, cor_behavior) #whatever you want




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### run correlation analysis #####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


############### + Pearson with false discovery rated adjusted #####################
bug_behav_fdr <- corr.test(cor_6g, cor_behavior, method = "pearson", adjust = "fdr")
bug_alpha_fdr <- corr.test(cor_6g, cor_alpha, method = "pearson", adjust = "fdr")
bug_bug_fdr <- corr.test(cor_6g, cor_2p, method = "pearson", adjust = "fdr")
#this all so analysis non-adjusted p values
#한 feature set만 넣으면 
r <- bug_behav_fdr$r
fdr <- bug_behav_fdr$p.adj
none <- bug_behav_fdr$p
#this is correlation matrix


####### + make coordination matrix ###########

coord_cor <- function(cor_fdr){
  coord = NULL
  for (i in 1:ncol(cor_fdr$r)) {
    for (j in 1: nrow(cor_fdr$r)){
      #get stars for p
      p <- cor_fdr$p[j,i]
      if (is.na(p)) { ps <- c("NA")
      } else if (p<=0.001) { ps <- c("***")
      } else if (0.001<=p&&p<0.01) { ps <- c("**")
      } else if (0.01<=p&&p<0.05) { ps <- c("*")
      } else ps <- c("ns")
      
      #get stars for fdr adjusted p
      f<- cor_fdr$p.adj[j,i]
      if (is.na(f)) { fdrs <- c("NA")
      } else if (f<=0.001) { fdrs <- c("***")
      } else if (0.001<=f&&f<0.01) { fdrs <- c("**")
      } else if (0.01<=f&&f<0.05) { fdrs <- c("*")
      } else fdrs <- c("ns")
      
      coord = rbind(coord, c(colnames(cor_fdr$r)[i], rownames(cor_fdr$r)[j], 
                             cor_fdr$r[j,i], p, f, ps, fdrs))
    }
  }
  return(coord)
}

coord_bug_behav <- coord_cor(bug_behav_fdr)
coord_bug_alpha <- coord_cor(bug_alpha_fdr)
coord_bug_bug <- coord_cor(bug_bug_fdr)
  
####### + make triangular coordination matrix ###########
  


##
