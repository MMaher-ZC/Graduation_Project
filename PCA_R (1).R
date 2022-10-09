#Data_Read
Data_read=read.csv(file = 'GSE157103_genes.ec.tsv', sep = '\t',row.names=1)
Namess=names(Data_read)
Data_read=as.matrix(Data_read)
for (i in 1:length(Namess)){
  if (grepl("NC",Namess[i], fixed=TRUE)){
    Namess[i]="NC"
  }
  else{
    Namess[i]="C"
  }
}

library(htmltools)
library(ggplot2)
library( "DESeq2" )
library(stringr)

library(janitor)

Data_read_DESEQ2=read.csv(file = 'GSE157103_genes.ec.tsv', sep = '\t',row.names=1)
Data_read_DESEQ2=round(Data_read_DESEQ2)
Condition=names(Data_read_DESEQ2)
Names=names(Data_read_DESEQ2)
Data_read_DESEQ2=as.matrix(Data_read_DESEQ2)
mode(Data_read_DESEQ2) <- "integer"


RRNA=read.csv(file = 'group-1379.csv')$Approved.symbol

K=0
GGenes=row.names(Data_read_DESEQ2)
inter=intersect(GGenes,RRNA)



Data_read_DESEQ2=Data_read_DESEQ2[!(row.names(Data_read_DESEQ2)%in%inter),]


#MetaData


for (i in 1:length(Condition)){
  if (grepl("NC",Condition[i], fixed=TRUE)){
    Condition[i]="NC"
  }
  else{
    Condition[i]="C"
  }
}


Matrix_Data=read.delim(file = 'GSE157103_series_matrix.txt', sep = '\t',skip=30)
Gender=Matrix_Data[10,2:127]
SAMPLES=colnames(Gender)
Gender=str_remove_all(Gender, "Sex: ")
Gender[115]="female"
Gender=as.data.frame(Gender)
colnames(Gender)=c("Gender")
rownames(Gender)=SAMPLES


Age=Matrix_Data[9,2:127]
x <- gregexpr("[0-9]+", Age)
Age <- unlist(regmatches(Age, x))
Age=as.numeric(unlist(Age))
Age

age=Age[order(Age)]
age


for (i in 1:length(Age)){
  if (Age[i]<60){
    Age[i]="Less Than 60"
  }
  else{
    Age[i]='60 or More '
  }
}

Age=as.data.frame(Age)
colnames(Age)=c("Age")
rownames(Age)=SAMPLES

Charlson=Matrix_Data[13,2:127]
x <- gregexpr("[0-9]+", Charlson)
Charlson <- unlist(regmatches(Charlson, x))
Charlson=as.numeric(unlist(Charlson))
for (i in 1:length(Charlson)){
  if (Charlson[i]<6){
    Charlson[i]="Mild-Moderate"
  }
  else{
    Charlson[i]='Severe'
  }
}

Charlson=as.data.frame(Charlson)
colnames(Charlson)=c("Charlson")
rownames(Charlson)=SAMPLES




MetaDataAll=cbind(Age,Gender,Condition,Charlson)
MetaDataAll$Age=factor(MetaDataAll$Age)
MetaDataAll$Gender=factor(MetaDataAll$Gender)
MetaDataAll$Condition=factor(MetaDataAll$Condition)
MetaDataAll$Charlson=factor(MetaDataAll$Charlson)

#DeSEQ2 and VST

dds1=DESeqDataSetFromMatrix(countData=Data_read_DESEQ2,colData=MetaDataAll,design = ~Condition)
dds1$Condition=relevel(dds1$Condition,"NC")
dds1 <- DESeq(dds1)

vsdata <- vst(dds1, blind=FALSE)




pcaData <- plotPCA(vsdata, intgroup = c( "Age",'Condition'), returnData = TRUE) # vsd and plotPCA are part of DESeq2 package, nothing with my example below. 
percentVar <- round(100 * attr(pcaData, "percentVar")) 


ggplot(pcaData, aes(x = PC1, y = PC2, color = Age , shape = Condition)) + 
  geom_point(size =3.5, aes(fill=Age)) + 
  geom_point(size =3.5) + 
  scale_shape_manual(values=c(24,22)) + scale_color_manual(values=c("red",'darkblue'),aesthetics = c("colour","fill"))+
  xlab(paste0('PC1: ', percentVar[1], '% variance')) + 
  ylab(paste0('PC2: ', percentVar[2], '% variance')) + guides(color=guide_legend(order=1),fill=guide_legend(order=1)) 




cor()
