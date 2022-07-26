library("edgeR")
library(RColorBrewer)
count<-read.table("filtered.brain.count",head=F,row.name=1)
group <-  c(1,1,1,2,2,2)
colnames(count)<-c("Brain.S2","Brain.S3","Brain.S5","Brain.S6","Brain.S7","Brain.S8")
y<- DGEList(counts=count,group=group)

#看样品间相关性PCA
#plotMDS(y,col=c("red","red","red","blue","blue","blue"),main="MDS for mC in Brain Gene")
#dev.off()

y <-estimateDisp(y)
et<- exactTest(y)

result = topTags(et, n = nrow(et$table))$table
result2 <-subset(result,PValue<0.01)

test<-merge(count,result2,by="row.names",sort=F)
test2<-test[order(test$FDR),]
#获得差异甲基化位点
write.csv(test2,file = "brain.DMS.csv")

#看DMSs在染色体上的分布
test3<-subset(test2,FDR<0.05)
library("stringr")
test3$Chr<-str_extract(test3$Row.names,"Chr\\d+")
table(test3$Chr)
plot(table(test3$Chr),col =brewer.pal(8, "Dark2"),xlab="Chr",ylab="Num of DMCs")
