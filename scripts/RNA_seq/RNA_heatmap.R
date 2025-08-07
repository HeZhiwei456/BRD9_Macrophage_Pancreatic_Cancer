library(pheatmap)
library(dplyr)
library(getopt)
library(ggplot2)
library(tibble)
library(textshape)
#print("Usage: Rscript reads_quality.R [-D Dir]")
options <- c("reads", "r", 2, "character", "reads file",
             "gene","g", 2, "character", "diffgene file",
             "top20","t", 2, "character", "top 20 genes",
             "control","C", 2, "character","control",
             "experiment","E", 2, "character", "experiment",
             "control_name","c", 2, "character", "control name",
             "experiment_name", "e", 2, "character", "experiment name",
             "output","o", 2, "character","output dir",
             "name","n", 2, "character", "output name")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
dir.create(opt$output, recursive = TRUE)

#-----------------------------------------------------
#get normalized count
#opt$count = "/home/zm/Analysis/RNA-seq/count/All_normalized_count.xls"
count <- read.table(opt$reads, sep="\t", header = T, check.names = F)
gene <- read.table(opt$gene, sep="\t", header =T)
up_gene <- gene %>% filter(log2FoldChange > 0) %>% top_n(n = 20, wt = -log10(pvalue))
down_gene <- gene %>% filter(log2FoldChange < 0)  %>% top_n(n = 20, wt = -log10(pvalue))
top20 <- rbind(up_gene, down_gene)
if (opt$top20 == "T"){
  gene <- top20
}

gene <- gene[,c(1:2)]
count <- merge(count, gene, by=colnames(gene))
#head(count)
ctr <- strsplit(opt$control, ",", fixed = TRUE)[[1]]
exp <- strsplit(opt$experiment, ",", fixed = TRUE)[[1]]
ctr_num <- length(ctr)
exp_num <- length(exp)
count <- count[!duplicated(count$SYMBOL),]
count <- na.omit(count)
count <- column_to_rownames(count, "SYMBOL")
normalized_count <- cbind(count[,ctr], count[,exp])
head(normalized_count)
#-----------------------------------------------------
#annotation information
annotation <- data.frame(
  type = c(rep(opt$control_name, ctr_num), rep(opt$experiment_name, exp_num)),
  sample = colnames(normalized_count))
rownames(annotation) <- colnames(normalized_count)

#-----------------------------------------------------
#plot with pheatmap
bk <- c(seq(-3, 3,by=0.01))
color = c(colorRampPalette(colors = c("#57C3F3", "#F7F398"))(length(bk)/2),colorRampPalette(colors = c("#F7F398","#E95C59"))(length(bk)/2))
rownames = F
if (opt$top20 == "T"){
  rownames = T
}

pdf(paste0(opt$output,"/", opt$name, ".pdf"))
pheatmap(
  normalized_count,
  scale = "row",
  annotation_col = annotation,
  show_rownames = rownames, 
  show_colnames = T,
  #color = color,
  #breaks = bk,
  cluster_rows = T,
  cluster_cols = FALSE,
  annotation_legend = T
)
dev.off()

png(paste0(opt$output,"/", opt$name, ".png"))
pheatmap(
  normalized_count,
  scale = "row",
  annotation_col = annotation,
  show_rownames = rownames, 
  show_colnames = T,
  #color = color,
  #breaks = bk,
  cluster_rows = T,
  cluster_cols = FALSE,
  annotation_legend = T
)
dev.off()
