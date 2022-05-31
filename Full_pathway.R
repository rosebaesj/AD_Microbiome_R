#PICRUST
library("phyloseq")
library("ggplot2") 
library("microbiomeMarker")


########################################################################
############################################################
#* 기존 PICRUST1에서는 바로 pathway collapse가 python 코드로 실행이 가능했었다
#* 하지만 PICRUST2는 지원하지 않음...
#* 그래서 엄청나게 자잘한 여러개의 pathway로 구성되는 문제점이 있음
#* 그래서 아얘 R에서 pathway collapse를 코딩해야함... ㅠㅠ
#* 그래도 git 상에 PICRUST 개발자들에의해 제공되는 (보증되지는 않는) 코드가 있어서 아래 시행함
########################################################################



setwd("All_AD")



########################################################################
#*** categorize_by_function.py 을 대체할 R function****
########################################################################

### Reproducing the categorize by function (level 3) functionality in plain-text tables.
### Doing this because adding a column of KEGG Pathways to a table and then converting
### that table to BIOM is difficult.

categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
  
  colnames(out_pathway) <- c("pathway", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}




##### L2

categorize_by_function_l2 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
  
  colnames(out_pathway) <- c("pathway", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}









########################################################################
#* 실제 실행 CATEGORIZATION ****
#*#*** INPUT ****
########################################################################

### Example commands:
### Read in BRITE hierarchy per KO.
#*KEGG categorization database 정보.
#*https://www.dropbox.com/s/a5o4li0irsqupt3/picrust1_KO_BRITE_map.tsv?dl=1.
#*위의 주소에서 다운받을 수 있음

kegg_brite_map <- read.table("picrust/picrust1_KO_BRITE_map.tsv",
                             header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)

#*data 가져오기. K 코드로 되어있는거면 됨. description 필요 없음
#*QIIME에서 PICRUST2 실행, picrust2_pipeline.py \실행 하면 나옴.
### When reading in tab-delimited file of KO predictions (PICRUSt2 output):
test_ko <- read.table("picrust/pred_metagenome_unstrat.tsv", header=TRUE, sep="\t", row.names=1)


### Alternatively, when reading in legacy TSV BIOM file (PICRUSt1 output): 
### test_ko <- read.table("/path/to/test_ko.tsv",
###                       header=TRUE, sep="\t", row.names=1, skip=1, comment.char="")
### if(length(which(colnames(test_ko) == "KEGG_Pathways")) > 0)) {
###     test_ko <- test_ko[, -which(colnames(test_ko) == "KEGG_Pathways")]
### }




########################################################################
#*** 실제 실행 CATEGORIZATION ****
#*#*** 얘는 L3로 output 하는 법 ****
########################################################################


### Run function to categorize all KOs by level 3 in BRITE hierarchy.
test_ko_L3 <- categorize_by_function_l3(test_ko, kegg_brite_map)

test_ko_L3_sorted <- test_ko_L3[rownames(orig_ko_L3), ]
#
#
### Commands that could be used to compare the KO levels from this function with the actual output of categorize_by_function.py:
# orig_ko_L3 <- read.table("/path/to/test_ko_L3.tsv",
#                          header=TRUE, sep="\t", row.names=1, skip=1, comment.char="", quote="")
# 
# orig_ko_L3 <- orig_ko_L3[, -which(colnames(orig_ko_L3) == "KEGG_Pathways")]
# 
# orig_ko_L3 <- orig_ko_L3[-which(rowSums(orig_ko_L3) == 0),]
#
#
### The below command will be True when the output is exactly the same.
# identical(test_ko_L3_sorted, orig_ko_L3)







########################################################################
#*** LEFSE ****
########################################################################


##수기로 첫번째 column(code만 있는 것) 삭제하여 txt로 저장.
##1.  PATHWAY with METACYC

test_ko_L3 <- rownames_to_column(test_ko_L3, "pathway")

write.table(test_ko_L3,"picrust2/test_ko_L3.tsv", row.names = FALSE, sep = "\t")

pathway <- import_picrust2("picrust2/test_ko_L3.tsv",
                           "phyloseq/metadata.tsv", #metadata
                           trait = "KO") #pathway data base type

pathway_lefse <- run_lefse(pathway,"group")

pathway_lefse_data <-data.frame(marker_table(pathway_lefse))
plot_ef_bar(pathway_lefse)+
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


ggsave("picrust2/pathway_LEfSe.png", width=8, height=6, units="in", device = "png")


########################################################################
#*** Cladogram? ****
########################################################################

#* pathway analysis도 cladogram (3level) 가능한가봄. 하지만. 방법은 모름
#* 

