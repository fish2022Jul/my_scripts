library("tidyverse")
library("org.Anigrocauda.eg.db")
library("dplyr")
library("DESeq2")
library("clusterProfiler")
library("KEGGREST")
library("topGO")
library("stringr")
library("pathview")


#setwd("e:/database/black_tail_genome/2023.hypoxia/Result/03.expression/gene/")

read_num<-as_tibble(read.csv("count.csv"))
brain_count<-read_num %>% column_to_rownames(.,var="ID") %>% select(ends_with("B"))
heart_count<-read_num %>% column_to_rownames(.,var="ID") %>% select(ends_with("H"))
colnames(brain_count) <- c("0h_1B","0h_2B","0h_3B","1.5h_1B","1.5h_2B","1.5h_3B","Re1.5h_1B","Re1.5h_2B","Re1.5h_3B","3h_1B","3h_2B","3h_3B","Re3h_1B","Re3h_2B","Re3h_3B","4.5h_1B","4.5h_2B","4.5h_3B","Re4.5h_1B","Re4.5h_2B","Re4.5h_3B")
colnames(heart_count) <- c("0h_1H","0h_2H","0h_3H","1.5h_1H","1.5h_2H","1.5h_3H","Re1.5h_1H","Re1.5h_2H","Re1.5h_3H","3h_1H","3h_2H","3h_3H","Re3h_1H","Re3h_2H","Re3h_3H","4.5h_1H","4.5h_2H","4.5h_3H","Re4.5h_1H","Re4.5h_2H","Re4.5h_3H")


source("./enrich_topGO_KEGGREST.R")

#DE analysis for heart tissue

#DE analysis for brain tissue
#rna_groups(brain_count,"1.5h","4.5h")
#print("maker 6!")
rna_groups(brain_count,"Re1.5h","Re3h")
print("maker 7!")

rna_groups(brain_count,"Re1.5h","Re4.5h")
print("maker 8!")

#rna_groups(brain_count,"3h","Re3h")
#rna_groups(brain_count,"Re3h","Re4.5h")
#rna_groups(brain_count,"1.5h","Re1.5h")
#rna_groups(brain_count,"4.5h","Re4.5h")
print("All done!")
