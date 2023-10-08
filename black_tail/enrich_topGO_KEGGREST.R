
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
str2<-deparse(substitute(count))
tissue=str_split(str2,"_")[[1]][1]
dir_pos=paste(tissue,treat,control,sep="_")
print(dir_pos)
if(!dir.exists(dir_pos))
    {dir.create(dir_pos)
    }
setwd(dir_pos)
write.csv(count2,file = paste0(file_name,".bg.csv"))
write.csv(diff_data,file = paste0(file_name,".csv"),row.names = F)
save.image(file=paste0(file_name,".Rdata"))


rld <- rlog(dds)
tiff(file=paste0(file_name,"_PCA",".tiff"),res=300,unit="mm",width=160,height =100,compression = "lzw")
myp1<-plotPCA(rld, intgroup="condition")
print(myp1)
dev.off()



# topGO Enrichment analysis
deg<-read.csv(file = paste0(file_name,".csv"))
pv2<-deg$pvalue
names(pv2)<-deg$Row.names
go2gene<-readMappings("../../00.data/go2genes.txt")

selection <- function(x) TRUE
GOdata <- new("topGOdata", ontology="BP", allGenes=pv2,annot=annFUN.GO2genes, GO2genes=go2gene,geneSel=selection, nodeSize=3)
#finally got it ready!!!

resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = 15, numChar = 120)
write.csv(tab,file = paste0(file_name,".topGO.csv"),row.names = F)

# Assuming 'resultTable' is the output from GenTable
significantGOs <- tab$GO.ID
genesForAllTerms <- genesInTerm(GOdata, significantGOs)


df <- data.frame(
		   GO = rep(names(genesForAllTerms), sapply(genesForAllTerms, length)),
		     Gene = unlist(genesForAllTerms, use.names = FALSE)
		   )

colnames(df)<-c("GO","geneID")


b2z_data <- read.table("../../00.data/black2zebra_gene_name.txt", header =FALSE, sep = " ", stringsAsFactors = FALSE,fill=TRUE)
colnames(b2z_data)<- c("geneID", "col2", "geneName")
merged_df <- df %>%
	  left_join(b2z_data[, c("geneID", "geneName")], by = "geneID",relationship = "many-to-many")

  merged_df$geneName[is.na(merged_df$geneName)] <- ""

  write.csv(merged_df, "merged_output.csv", row.names = FALSE)

# Create the rainbow bar plot with ggplot2
 num_terms <- nrow(tab)
tab$Term<- stringr::str_wrap(tab$Term, width = 30)
plot <- ggplot(tab, aes(x = reorder(Term,-as.numeric(raw.p.value)), y = -log10(as.numeric(raw.p.value)), fill = Term)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rainbow(num_terms)) +  # Use rainbow palette for colors
  coord_flip() +
  labs(y = "-log10(p-value)", x = "GO Term", title = "TopGO Enrichment Analysis") +
  theme_minimal() +
  guides(fill = "none")

tiff(filename = "GO.bar.tiff",res = 300,width = 160,height = 100,compression = "lzw",units = "mm")
par(mgp=c(1.5,0.5,0),cex.axis=0.7,mar=c(2.5,2.5,2,0.5)+0.1)
print(plot)
dev.off()


#KEGG   FOR NON-MODEL   https://zhuanlan.zhihu.com/p/561522453
#高级可视化     https://www.jianshu.com/p/5a7eaa8fa4b8



kegg_all<-read.csv("../kegg_all.csv")
TERM2GENE = kegg_all[c('PathwayID', 'SeqID')]
result_list <- split(TERM2GENE$SeqID, TERM2GENE$PathwayID)

TERM2NAME = kegg_all[c('PathwayID', 'C')]
pathways.list<- split(TERM2NAME$C, TERM2NAME$PathwayID)
p2.list<-lapply(pathways.list, unique)

pVals.by.pathway <- t(sapply(names(result_list),
	function(pathway) {
	pathway.genes <- result_list[[pathway]]
	list.genes.in.pathway <- intersect(names(pv2), pathway.genes)
	list.genes.not.in.pathway <- setdiff(names(pv2), list.genes.in.pathway)
	scores.in.pathway <- pv2[list.genes.in.pathway]
	scores.not.in.pathway <- pv2[list.genes.not.in.pathway]
    if (length(scores.in.pathway) > 0){
  	 p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
	} else{													        p.value <- NA
		}
	return(list(p.value = p.value, 
                    Annotated = length(list.genes.in.pathway),     
                    Genes = list.genes.in.pathway))
			}
				))
oo<-as.data.frame(pVals.by.pathway)
oo<-subset(oo,p.value<0.1)
oo$pathway.name <- p2.list[rownames(oo)]
oo$p.value=as.numeric(oo$p.value)
oo<-oo[order(oo$p.value),]
list_cols <- sapply(oo, class) == "list"
list_cols
oo[list_cols] <- lapply(oo[list_cols], function(col) {
		       sapply(col, function(cell) {
		             paste(unlist(cell), collapse = ",")
							       })
		   })
oo<-oo[,c("p.value","pathway.name","Annotated","Genes")]
write.csv(oo,file=paste0(file_name,".KEGG.csv"),row.names = F)



###draw KEGG barplot
if(nrow(oo)>0){
plotdat <- oo
plotdat$Annotated <- as.numeric(plotdat$Annotated)
plotdat$p.value <- as.numeric(plotdat$p.value)
pp <- ggplot(plotdat, aes(x = -log10(p.value), y = reorder(pathway.name, -log10(p.value)), fill = Annotated)) +
  geom_bar(stat = "identity") +
  labs(x = "-log10(P-Value)", y = NULL, fill = "# Genes") +
  scale_fill_gradient(low = "red", high = "blue")
tiff(filename = "KEGG.bar.tiff",res = 300,width = 160,height = 100,compression = "lzw",units = "mm")
par(mgp=c(1.5,0.5,0),cex.axis=0.7,mar=c(2.5,2.5,2,0.5)+0.1)
print(pp)
dev.off()
}

######################
print(paste0(file_name,"is done!"))
setwd("../")

}

#
