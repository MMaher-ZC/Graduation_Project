res <- read.csv("F:\\ALL_MIRNA.csv", header=TRUE)
#res <- read.csv("F:\\DEGS_DESEQ2_Two_Tailed22.csv", header=TRUE)

head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(min(res$log2FoldChange, na.rm = T)-1,
                                           max(res$log2FoldChange,na.rm=T)+1),
               ylim=c(-log10(max(res$padj,na.rm=T))-0.5,-log10(min(res$padj,na.rm=T))+0.5),
                      col="grey"))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)

with(subset(res,padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=20, col="darkblue"))
with(subset(res, padj<.05 & log2FoldChange>.5), points(log2FoldChange, -log10(padj), pch=20, col="darkred"))
legend("right",legend=c("Upregulated",'Downregulated',"Not Significant"),col=c("darkred",'darkblue',"grey"),fill =c('darkred',"darkblue",'grey') )


     