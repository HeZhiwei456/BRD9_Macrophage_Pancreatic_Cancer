library(pathview)
library(ggplot2)
library(gage)
library(gageData)
library(dplyr)
library(MetBrewer)
library(ggnewscale)
library(clusterProfiler)
library(DOSE)
library(stringr)
library(RColorBrewer)
library(getopt)
options <- c("gene", "g", 2, "character", "diffgene file",
             "species", "s", 2, "character", "human, mouse or rat",
             "output", "o", 2, "character", "output dir",
             "name","n", 2, "character", "output name")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
dir.create(opt$output, recursive = TRUE)
gene <- read.csv(opt$gene, sep="\t", header = T)
if (opt$species == "human" | opt$species == "hs"){
  library(org.Hs.eg.db)
  ORG <- org.Hs.eg.db
}else if(opt$species == "mouse" | opt$species == "mm"){
  library(org.Mm.eg.db)
  ORG <- org.Mm.eg.db
}else if(opt$species == "rat"){
  library(org.Rn.eg.db)
  ORG <- org.Rn.eg.db
}else{
  print("species should be human, mouse or rat")
  quit()
}

entrezID = select(ORG, keys=unique(gene$SYMBOL), columns = "ENTREZID", keytype="SYMBOL")
GO <- enrichGO(gene = entrezID$ENTREZID,OrgDb = ORG,pvalueCutoff =0.05,qvalueCutoff = 0.05,ont="all",readable =T)
num <- nrow(as.data.frame(GO))

if(num == 0){
  GO <- enrichGO(gene = entrezID$ENTREZID,OrgDb = ORG,pvalueCutoff =1,qvalueCutoff = 1,ont="all",readable =T)
}
pdf(paste0(opt$output,"/", opt$name, ".pdf"))
dotplot(GO,showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free', space = "free")
dev.off()
write.table(GO, file=paste0(opt$output,"/", opt$name,"_enrichment.xls"),sep="\t",quote=F,row.names = F)      

png(paste0(opt$output,"/", opt$name, ".png"))
dotplot(GO,showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free', space = "free")
dev.off()

#UD_gene <- list(Up = up_gene$SYMBOL, Down = down_gene$SYMBOL)
#CompGO <- compareCluster(UD_gene, fun="enrichGO", OrgDb = ORG, keyType = "SYMBOL", ont = "BP", pvalueCutoff=0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05) 
#pdf(paste0(opt$output,"/Fig6.pdf"))
#dotplot(CompGO, showCategory = 10, includeAll=TRUE)
#dev.off()
#write.table(CompGO, file=paste0(opt$output,"/CompGO_enrichment.xls"),sep="\t",quote=F,row.names = F)      
GO_item <- as.data.frame(GO)
geneID <- str_split(GO_item$geneID, "/")
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
GO_item$geneID <- ID

enrich <- GO_item %>% group_by(ONTOLOGY) %>% top_n(n = 5, wt = -pvalue) %>% filter(ONTOLOGY %in% c("BP", "CC", "MF")) %>% as.data.frame()
dt <- enrich
dt <- dt[order(dt$ONTOLOGY),]
dt$Description <- factor(dt$Description, levels = dt$Description)
num <- nrow(dt)
max_value <- max(-log10(dt$pvalue))
p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description))) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8, fill = "#6BB9D2") +
  geom_point(aes(size = Count, color = p.adjust, fill = p.adjust)) +
  scale_fill_gradientn(colors = met.brewer("Cassatt1")) +
  scale_color_gradientn(colors = met.brewer("Cassatt1")) +
  scale_x_continuous(expand = c(0,0), limits = c(0, max_value*1.1)) +
  labs(x = bquote(-Log[10](pvalue)), y = "", title = "GO Pathway enrichment") +
  geom_text(size = 3.8, aes(x = 0.05, label = Description), hjust = 0) +
  geom_text(size = 2.5, aes(x = 0.05, label = geneID), hjust = 0, vjust = 3.5) +
 
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
        axis.title = element_text(color = "black", size = 10),
        axis.text = element_text(color = "black",size = 10),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0,"cm"),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        legend.title = element_text(color = "black",size = 10),
        legend.text = element_text(color = "black", size = 10),
        plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))+
  facet_grid(ONTOLOGY~., scale='free',space = "free")
ggsave(p, file=paste0(opt$output,"/", opt$name, "_v2.pdf"), width = 6, height = num*0.5)
ggsave(p, file=paste0(opt$output,"/", opt$name, "_v2.png"), width = 6, height = num*0.5)
