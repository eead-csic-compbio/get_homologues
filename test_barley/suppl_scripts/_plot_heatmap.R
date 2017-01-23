library(gplots)
library(cluster)

data = read.table(file="accessory_stats_min5.tab",header=F,sep="\t")
rnames <- data[,7]
mat_data <- data.matrix(data[,1:3])
rownames(mat_data) <- rnames
#colnames(mat_data) <- c("73","SC","Hs")
palette <- colorRampPalette(c("lightyellow","red"))(n = 20)
pdf("accessory_stats_min5.pdf")
heatmap.2(mat_data,cellnote=mat_data,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",dendrogram="none",margins =c(1,25),colsep=c(0:3),rowsep=c(0:40),sepwidth=c(0.01,0.01), sepcolor="grey", col=palette, lhei = c(1,10),labCol='')
dev.off()
