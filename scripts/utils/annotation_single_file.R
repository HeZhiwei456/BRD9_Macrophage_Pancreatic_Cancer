library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(ggupset)
library(ggplot2)
library(stringr)
library(getopt)
library(dplyr)
library(DOSE)
library(GseaVis)
options <- c("input", "i", 2, "character", "input dir",
             "genome","g", 2,"character","human or mouse",
             "output","o", 2, "character", "output dir",
             "name","n", 2, "character", "output name")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

#opt$input="/home/biology/Analysis/ATAC-seq/20241102_庄铨南_复旦大学附属儿科医院/diffpeak/WT_S34/diffpeak.xls"

if(!file.exists(opt$output)){
  dir.create(opt$output)
}

name <- str_split_fixed(opt$input,"[.]", n=2)[1]
name <- str_split(name,"/")[[1]]
name <- name[length(name)]
peak <- readPeakFile(opt$input, header = T)
COLnames <- colnames(as.data.frame(peak))
print(COLnames)
if(COLnames[7] == "X."){
  peak <- readPeakFile(opt$input, header = F)
  print(head(peak))
}
#opt$species='hg19'
if (opt$genome=="hg38"){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(EnsDb.Hsapiens.v86)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ORG <- "org.Hs.eg.db"
  organism = "hsa"
}else if(opt$genome=="hg19"){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(EnsDb.Hsapiens.v75)
  library(org.Hs.eg.db)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  ORG <- "org.Hs.eg.db"
  organism = "hsa"
}else if(opt$genome=="mm10"){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(EnsDb.Mmusculus.v75)
  library(org.Mm.eg.db)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- "org.Mm.eg.db"
  organism = "mmu"
}

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb=ORG)

pdf(paste0(opt$output,"/",opt$name, "_annoPie.pdf"))
plotAnnoPie(peakAnno)
dev.off()

png(paste0(opt$output,"/",opt$name, "_annoPie.png"))
plotAnnoPie(peakAnno)
dev.off()

pdf(paste0(opt$output,"/",opt$name, "_annoUpset.pdf"))
upsetplot(peakAnno, vennpie=F)
dev.off()

png(paste0(opt$output,"/",opt$name, "_annoUpset.png"))
upsetplot(peakAnno, vennpie=F)
dev.off()

peakAnno <- as.data.frame(peakAnno)
if(opt$genome=="hg19"){
  trans <- bitr(peakAnno$geneId, fromType = "ENTREZID", toType = c("ENSEMBL","SYMBOL","GENENAME"), OrgDb =  ORG)
  peakAnno <- merge(peakAnno, trans, by.x = "geneId", by.y = "ENTREZID")
  peakAnno <- peakAnno[,c(2:20,1,21:25)]
  colnames(peakAnno) <- c("seqnames","start","end","widths","strand","width","strands","Conc","Conc_Treatment","Conc_Control","Fold","p.value","FDR","annotation","geneChr","geneStart","geneEnd","geneLength","geneStrand","geneId","transcriptId","distanceToTSS","ENSEMBL","SYMBOL","GENENAME")
}else{
  colnames(peakAnno) <- c("seqnames","start","end","widths","strand","width","strands","Conc","Conc_Treatment","Conc_Control","Fold","p.value","FDR","annotation","geneChr","geneStart","geneEnd","geneLength","geneStrand","geneId","transcriptId","distanceToTSS","ENSEMBL","SYMBOL","GENENAME")
}
peakAnno$start <- peakAnno$start - 1
peakAnno <- peakAnno[,c(1:3,11:ncol(peakAnno))]
write.table(peakAnno, file=paste0(opt$output, "/", opt$name, "_anno.xls"), sep="\t", row.names=F, quote = F)

promoter_gene <- peakAnno[grep("Promoter", peakAnno$annotation),]$geneId
compKEGG <- enrichKEGG(promoter_gene,
                       organism = organism,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH")

length <- as.data.frame(compKEGG)
if(nrow(length) > 3){
  p <- dotplot(compKEGG, showCategory = 10, font.size = 7, title = "KEGG Pathway Enrichment Analysis")
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_KEGG.pdf"))
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_KEGG.png"))
  symbol_kegg <- setReadable(compKEGG, OrgDb=ORG, keyType="ENTREZID")
  write.table(symbol_kegg, file=paste0(opt$output, "/", opt$name, "_KEGG.xls"), row.names=F, sep= "\t")
}else{
  compKEGG <- enrichKEGG(promoter_gene,
                         organism = organism,
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1,
                         minGSSize = 2,
                         pAdjustMethod = "BH")
  
  p <- dotplot(compKEGG, showCategory = 10, font.size = 7, title = "KEGG Pathway Enrichment Analysis")
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_KEGG.pdf"))
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_KEGG.png"))
  symbol_kegg <- setReadable(compKEGG, OrgDb=ORG, keyType="ENTREZID")
  write.table(symbol_kegg, file=paste0(opt$output, "/", opt$name, "_KEGG.xls"), row.names=F, sep= "\t")
}

kegg_df <- data.frame(symbol_kegg) %>% head(10)
p <- sankeyGoPlot(goData = kegg_df)
ggsave(p, file=paste0(opt$output, "/", opt$name, "_KEGG_v2.pdf"))
ggsave(p, file=paste0(opt$output, "/", opt$name, "_KEGG_v2.png"))

GO <- enrichGO(promoter_gene,
               OrgDb = ORG,
               keyType = "ENTREZID",
               ont = "ALL",
               pvalueCutoff=0.05,
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05)

length <- as.data.frame(GO)
if(nrow(length) > 0){
  p <- dotplot(GO, showCategory = 7,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_GO.pdf"))
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_GO.png"))
  symbol_go <- setReadable(GO, OrgDb=ORG, keyType="ENTREZID")
  write.table(symbol_go, file=paste0(opt$output, "/", opt$name, "_GO.xls"), row.names=F, sep= "\t")
}else{
  GO <- enrichGO(promoter_gene,
                 OrgDb = ORG,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pvalueCutoff=1,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 1)
  p <- dotplot(GO, showCategory = 7,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_GO.pdf"))
  ggsave(p, file=paste0(opt$output, "/", opt$name, "_GO.png"))
  symbol_go <- setReadable(GO, OrgDb=ORG, keyType="ENTREZID")
  write.table(symbol_go, file=paste0(opt$output, "/", opt$name, "_GO.xls"), row.names=F, sep= "\t")
  
}
  
