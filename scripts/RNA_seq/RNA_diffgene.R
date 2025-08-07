library(DESeq2)
library(tidyverse)
library(corrplot)
library(dplyr)
library(getopt)
library(tibble)
library(ggrepel)
library(pheatmap)
library(textshape)
#print("Usage: Rscript reads_quality.R [-D Dir]")
options <- c("read", "r", 2, "character", "read file",
             "control","C", 2, "character","control name",
             "experiment","E", 2, "character", "experiment",
             "control_name", "c", 2, "character", "experiment",
             "experiment_name", "e", 2, "character", "experiment",
             "fold", "f", 2, "numeric", "log2 fold change",
             "method", "m", 2, "character", "padj or pvalue",
             "pvalue","p", 2, "numeric", "p value",
             "output","o", 2, "character","output dir")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
dir.create(opt$output, recursive = TRUE)

count <- read.csv(opt$read, sep="\t", header = T, check.names = F)
count <- count[!duplicated(count$Geneid),]
count <- count[!is.na(count$Geneid),]

#colnames(count)<- gsub("[.]", "-", colnames(count))
count <- count[rowSums(count[,9:ncol(count)])>0,]
count <- column_to_rownames(count, "Geneid") #textshape::
meta_count <- count[,1:7]
reads_count <- count[,8:ncol(count)]
ctr <- strsplit(opt$control, ",", fixed = TRUE)[[1]]
exp <- strsplit(opt$experiment, ",", fixed = TRUE)[[1]]
ctr_num <- length(ctr)
exp_num <- length(exp)

raw_count <- cbind(reads_count[,ctr], reads_count[,exp])
raw_count <- raw_count[rowSums(raw_count)>0,]
condition <- factor(c(rep(opt$control_name,ctr_num), rep(opt$experiment_name,exp_num)))
coldata <- data.frame(row.names=colnames(raw_count), condition)
dds <- DESeqDataSetFromMatrix(countData=raw_count, colData=coldata, design=~condition)
rld <- rlogTransformation(dds)
dds_size <- estimateSizeFactors(dds)

normalized_count = counts(dds_size, normalized = TRUE) # normalized count
normalized_count <- as.data.frame(normalized_count)

dds <- DESeq(dds)
resdata <- results(dds)
head(resdata)
summary(resdata)
resdata <- as.data.frame(resdata)
resdata <- rownames_to_column(resdata, "GeneID")

CTRL <-rowMeans((normalized_count[,ctr])) %>% as.data.frame()
colnames(CTRL) <- paste0("mean(", opt$control_name, ")")
EXP <-rowMeans((normalized_count[,exp])) %>% as.data.frame()
colnames(EXP) <- paste0("mean(", opt$experiment_name, ")")
MEAN <- cbind(CTRL, EXP)
MEAN <- round(MEAN, 2)
MEAN <- rownames_to_column(MEAN, "GeneID")

resdata <- merge(MEAN, resdata, by = "GeneID")
meta_count <- rownames_to_column(meta_count, "GeneID")
resdata <- merge(meta_count, resdata, by="GeneID") #
resdata <- resdata[,c(1,2,9:ncol(resdata))]
if(resdata[1,3] > resdata[1,4]  & resdata[1,6] > 0){
  resdata[,6] <- -resdata[,6]
}
if(resdata[1,4] > resdata[1,3]  & resdata[1,6] < 0){
  resdata[,6] <- -resdata[,6]
}
#------------------------------------------------------------
#write out data
write.table(resdata, file=paste0(opt$output,"/All_gene.xls"),sep = "\t", row.names=F,quote=F)
if(opt$method == "padj"){
  diff_gene <-subset(resdata, padj < opt$pvalue & abs(log2FoldChange) > opt$fold)
}else if(opt$method == "pvalue"){
  diff_gene <-subset(resdata, pvalue < opt$pvalue & abs(log2FoldChange) > opt$fold)
}else if(opt$method == "FDR"){
  diff_gene <-subset(resdata, FDR < opt$pvalue & abs(log2FoldChange) > opt$fold)
}
write.table(diff_gene, file=paste0(opt$output,"/diff_gene.xls"),sep = "\t",row.names=F,quote=F)

#------------------------------------------------------------
#get normalized count
normalized_count <- rownames_to_column(normalized_count, "GeneID")
normalized_count <- merge(meta_count, normalized_count, by="GeneID")
normalized_count <- normalized_count[,c(1:2,9:ncol(normalized_count))]
write.table(normalized_count, file=paste0(opt$output,"/All_normalized_count.xls"),sep = "\t", row.names=F,quote=F)

#------------------------------------------------------------
#PCA plot
pcaData <- plotPCA(rld, returnData=TRUE)
#pcaData[,3:4] <- as.factor(pcaData[,3:4])
percentVar <- round(100 * attr(pcaData, "percentVar"))
pc1 <- as.character(percentVar[1])
pc2 <- as.character(percentVar[2])

if(ctr_num > 3 & exp_num >3){
  p <- ggplot(pcaData, aes(PC1, PC2, color=group)) + #, shape=name
    geom_point(size=3) +
    geom_text_repel(label = pcaData$name, size = 5)+
    stat_ellipse(aes(color = group), level = 0.90, show.legend = FALSE)+
    xlab(paste0("PC1: ", pc1, "% variance")) +
    ylab(paste0("PC2: ", pc2, "% variance")) +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
}else{
  p <- ggplot(pcaData, aes(PC1, PC2, color=group)) + #, shape=name
    geom_point(size=3) +
    geom_text_repel(label = pcaData$name, size = 5)+
    #stat_ellipse()+
    #stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)+
    xlab(paste0("PC1: ",pc1,"% variance")) +
    ylab(paste0("PC2: ",pc2,"% variance")) +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"))
}

ggsave(p, file=paste0(opt$output,"/PCA.pdf"))
ggsave(p, file=paste0(opt$output,"/PCA.png"))

#------------------------------------------------------------
#sample correlation
vsd <- vst(dds_size, blind = T)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat, method = "pearson")

pdf(paste0(opt$output,"/correlation_pearson.pdf"))
pheatmap(vsd_cor,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = F)
dev.off()

png(paste0(opt$output,"/correlation_pearson.png"))
pheatmap(vsd_cor,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = F)
dev.off()
