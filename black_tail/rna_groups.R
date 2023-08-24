
rna_groups <- function(count,treat,control){
count2<-count %>% dplyr::select(starts_with(c(control,treat))) %>% filter(rowSums(. > 5) > 2)

condition<-factor(rep(c(control,treat),each=3))

col2<-data.frame(colnames(count2),condition)
dds <- DESeqDataSetFromMatrix(countData = round(count2),colData=col2,design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition",treat,control))
res <- res[order(res$padj),]
 diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
 #diff_gene <- row.names(diff_gene)
 idx <- which(dds$condition %in% c(treat,control))
diff_data <- merge(as.data.frame(diff_gene), as.data.frame(counts(dds, normalized=TRUE)[,idx]),by="row.names",sort=FALSE)
#resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)[,idx]),by="row.names",sort=FALSE)
#file_name="DEGs_"+treat+"_vs_"+control
file_name=paste("DEGs",treat,"vs",control,sep="_")
dir_pos=paste(treat,control,sep="_")
if(!dir.exists(dir_pos))
    {dir.create(dir_pos)
    }
setwd(dir_pos)
write.csv(count2,file = paste0(file_name,".bg.csv"))
write.csv(diff_data,file = paste0(file_name,".csv"),row.names = F)



rld <- rlog(dds)
tiff(file=paste0(file_name,"_PCA",".tiff"),res=300,unit="mm",width=160,height =100,compression = "lzw")
myp1<-plotPCA(rld, intgroup="condition")
print(myp1)
dev.off()



# Enrichment analysis
genes <- as.character(rownames(diff_gene))
bg<-read.csv(file=paste0(file_name,".bg.csv"))[,"X"]
term2gene <-read.csv("../../00.data/all_genes.go.csv")
term2gene<-term2gene[,c(2,1)]
term2name<-read.csv("../../00.data/term2name.csv")

ego<-enricher(gene=genes,  pvalueCutoff = 0.05,pAdjustMethod = "BH",universe=bg,minGSSize = 2,maxGSSize = 500,qvalueCutoff = 0.2,TERM2GENE=term2gene,TERM2NAME = term2name)
#ego <- enrichGO(gene= genes, OrgDb=org.Anigrocauda.eg.db,universe=as.character(colnames(count2)), keyType= 'GID',ont = "ALL", pAdjustMethod = "none", minGSSize =2,
#                pvalueCutoff  = 0.05,qvalueCutoff  = 0.2)

if(is.character(ego$ID) & length(ego$ID)!=0){
    write.csv(ego, file=paste0(file_name,"_GO",".csv"),  quote = FALSE, row.names = FALSE)
    tiff(file=paste0(file_name,"GO_bar",".tiff"),res=300,unit="mm",width=160,height =100,compression = "lzw")

    myp2<-barplot(ego, showCategory=20, title="EnrichmentGO")
    print(myp2)
    dev.off()
}



#KEGG   FOR NON-MODEL   https://zhuanlan.zhihu.com/p/561522453
#高级可视化     https://www.jianshu.com/p/5a7eaa8fa4b8


kegg_all<-read.csv("../kegg_all.csv")
kegg_rich <- enricher(gene = genes,  #待富集的基因列表
    TERM2GENE = kegg_all[c('PathwayID', 'SeqID')],  
    TERM2NAME = kegg_all[c('PathwayID', 'Description')], 
    universe=as.character(colnames(count2)),minGSSize =2,
    pAdjustMethod = 'none', pvalueCutoff = 0.05,  qvalueCutoff = 0.2)
kegg_rich
tiff(file=paste0(file_name,".KEGG",".tiff"),res=300,unit="mm",width=160,height =100,compression = "lzw")
if(is.character(kegg_rich$ID) & length(kegg_rich$ID) !=0){
    myp3<-barplot(kegg_rich, showCategory=20,title="Enrichment_KEGG")
    print(myp3)
    write.csv(kegg_rich, file=paste0(file_name,"_KEGG",".csv"),  quote = FALSE, row.names = FALSE)
}
#kegg_res <- enrichKEGG(gene=genes, organism='Anigrocauda', keyType="GID")
dev.off()


######################

setwd("../")

}

#
