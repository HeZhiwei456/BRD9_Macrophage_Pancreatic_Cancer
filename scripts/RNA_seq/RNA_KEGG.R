library(pathview)
library(gage)
library(gageData)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(MetBrewer)
library(ggnewscale)
library(clusterProfiler)
library(DOSE)
library(stringr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(RColorBrewer)
library(getopt)
library(GseaVis)
library(stringr)
options <- c("gene", "g", 2, "character", "diffgene file",
             "species", "s", 2, "character", "human, mouse or rat",
             "Number", "N", 2, "numeric", "show number",
             "output","o", 2, "character", "output dir",
             "name","n", 2, "character", "output name")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
dir.create(opt$output, recursive = TRUE)
gene <- read.csv(opt$gene, sep="\t", header = T)

if (opt$species == "human" | opt$species == "hs"){
  ORG <- org.Hs.eg.db
  organism = "hsa"
}else if(opt$species == "mouse" | opt$species == "mm"){
  ORG <- org.Mm.eg.db
  organism = "mmu"
}else if(opt$species == "rat"){
  ORG <- org.Rn.eg.db
  organism = 'rno'
}else{
  print("species should be human, mouse or rat")
  quit()
}

entrezID = select(ORG, keys=unique(gene$SYMBOL), columns = "ENTREZID", keytype="SYMBOL")
kegg <- enrichKEGG(gene = entrezID$ENTREZID, keyType = 'kegg',organism = organism,pvalueCutoff = 0.05,pAdjustMethod  = "BH",qvalueCutoff  = 0.05)
if(opt$species == "mouse"){
  kegg@result$Description <- str_split_fixed(kegg@result$Description, " - Mus",2)[,1]
}
length <- as.data.frame(kegg)

if(nrow(length) < 5){
	kegg <- enrichKEGG(gene = entrezID$ENTREZID, keyType = 'kegg',organism = organism,pvalueCutoff = 1,qvalueCutoff  = 1)
}
length <- as.data.frame(kegg)
if(nrow(length) < opt$Number){
  opt$Number = nrow(length)
}
pdf(paste0(opt$output,"/", opt$name, ".pdf"))
dotplot(kegg, showCategory = opt$Number)
dev.off()
symbol_kegg <- setReadable(kegg, OrgDb=ORG, keyType="ENTREZID")
write.table(symbol_kegg, file=paste0(opt$output,"/", opt$name,"_enrichment.xls"),sep="\t",quote=F,row.names = F)      

png(paste0(opt$output,"/", opt$name, ".png"))
dotplot(kegg,showCategory = opt$Number)
dev.off()

df <- as.data.frame(symbol_kegg)
df <- df[,c(4,5,7,8,9,10,11)]
df$GeneRatio <- as.numeric(str_split_fixed(df$GeneRatio, "/", 2)[,1])/as.numeric(str_split_fixed(df$GeneRatio, "/", 2)[,2])

geneID <- str_split(df$geneID, "/")
num <- 10
ID <- NULL
for(i in 1:length(geneID)){
  Genes <- geneID[[i]]
  Genes_number <- length(Genes)
  if(Genes_number < num){
    j = Genes_number
  }else{
    j = num
  }
  for(z in 1:j){
    if(z == 1){
      GENE <- Genes[1]
    }else{
      GENE <- paste0(GENE, "/", Genes[z])
    }
  }
  ID[i] <- GENE
}
df$geneID <- ID

df <- df[order(df$Count, decreasing = T),]
df <- top_n(df, 10, -pvalue)
df$Description <- factor(df$Description, levels = df$Description)
max_value <- max(df$GeneRatio)
p <- ggplot(data = df, aes(x = GeneRatio, y = rev(Description))) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8, fill = "#6BB9D2") +
  geom_point(aes(size = Count, color = p.adjust, fill = p.adjust)) +
  scale_fill_gradientn(colors = met.brewer("Cassatt1")) +
  scale_color_gradientn(colors = met.brewer("Cassatt1")) +
  scale_x_continuous(expand = c(0,0), limits = c(0, max_value*1.1)) +
  labs(x = "GeneRatio", y = "", title = "KEGG Pathway enrichment") +
  geom_text(size = 4, aes(x = 0, label = Description), hjust = 0) +
  geom_text(size = 3, aes(x = 0, label = geneID), hjust = 0, vjust = 4) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0,"cm"),
        plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))
ggsave(p, file=paste0(opt$output,"/", opt$name, "_v2.pdf"))
ggsave(p, file=paste0(opt$output,"/", opt$name, "_v2.png"))

kegg_df <- data.frame(symbol_kegg) %>% head(10)
p <- sankeyGoPlot(goData = kegg_df)
ggsave(p, file=paste0(opt$output,"/", opt$name, "_v3.pdf"))
ggsave(p, file=paste0(opt$output,"/", opt$name, "_v3.png"))

file.remove('Rplots.pdf')

#pdf(paste0(opt$output,"/KEGG_bar.pdf"))
#barplot(kegg,showCategory = 10)
#dev.off()

#png(paste0(opt$output,"/KEGG_bar.png"))
#barplot(kegg,showCategory = 10)
#dev.off()

#Up_entrezID = select(ORG, keys=up_gene$ENSEMBL, columns = "ENTREZID", keytype="ENSEMBL")
#Down_entrezID = select(ORG, keys=down_gene$ENSEMBL, columns = "ENTREZID", keytype="ENSEMBL")
#UD_gene <- list(Up = Up_entrezID$ENTREZID, Down = Down_entrezID$ENTREZID)
#compKEGG <- compareCluster(geneCluster=UD_gene,fun= "enrichKEGG",organism = organism,pvalueCutoff  = 0.05,pAdjustMethod = "BH")
#pdf(paste0(opt$output,"/Fig8.pdf"))
#dotplot(compKEGG, showCategory = 10, font.size = 7, title = "KEGG Pathway Enrichment Analysis")
#dev.off()
#compsymbol_kegg <- setReadable(compKEGG, OrgDb=ORG, keyType="ENTREZID")
#write.table(compsymbol_kegg, file=paste0(opt$output, "/compKEGG_enrichment.xls"), quote=F, row.names=F, sep= "\t")
