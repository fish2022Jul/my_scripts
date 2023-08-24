
deg<-read.csv("DEGs_0h_vs_Re4.5h.csv")
pv2<-deg$pvalue
names(pv2)<-deg$Row.names


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
		} else{
			p.value <- NA
		}
		return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
	}
))


outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- p2.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat$pathway.name <- as.character(outdat$pathway.name)
outdat$pathway.name <- stringr::str_wrap(outdat$pathway.name, width = 30)
outdat <- outdat[order(outdat$p.value),]
write.csv(outdat,file = paste0(file_name,".KEGG.csv"),row.names = F)


###draw KEGG barplot
library(ggplot2)
library(stringr)
plotdat <- outdat[1:10,]

pp<-ggplot(plotdat, aes(x = -log10(p.value), y = reorder(pathway.name, -log10(p.value)), fill = Annotated)) + geom_bar(stat = "identity") + labs(x = "-log10(P-Value)", y = NULL, fill = "# Genes") + scale_fill_gradient(low = "red", high = "blue")
tiff(filename = "KEGG.bar.tiff",res = 300,width = 160,height = 100,compression = "lzw",units = "mm")
par(mgp=c(1.5,0.5,0),cex.axis=0.7,mar=c(2.5,2.5,2,0.5)+0.1)
print(pp)
dev.off()


