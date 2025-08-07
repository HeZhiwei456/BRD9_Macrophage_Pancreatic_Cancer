library(ggplot2)
library(ggseqlogo)
library(patchwork)
library(limma)
library(dplyr)
library(getopt)
options <- c("input", "i", 2, "character", "input dir",
             "number", "n", 2, "numeric", "motif number",
             "output","o", 2, "character", "output dir")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if (!file.exists(opt$output)){
  dir.create(opt$output)
}

file <- dir(opt$input)
file <- file[grep("info|logo|RV|similar", file, invert = T)]
number <- length(file)

if (number < opt$number){
  opt$number = number
}

plot = list()
Dfs <- NULL
for(i in 1:opt$number){
  mFile <- read.csv(paste0(opt$input, "/motif", i, ".motif"), sep = "\t", header = T, check.names = F)
  
  gene <- strsplit2(colnames(mFile)[2], "[_:(/]")[,2]
  log10pvalue <- colnames(mFile)[4] %>% as.numeric()
  sites <- strsplit2(colnames(mFile)[6], "[:.]")[,2]  %>% as.numeric()
  data <- c(gene, log10pvalue, sites)
  Dfs <- rbind(Dfs, data)
  
  mFile <- mFile[,1:4]
  mFile <- t(mFile) %>% as.matrix()
  rownames(mFile) <- c("A","C","G","T")
  
  plot[gene] <- list(mFile)
}

logo <- ggplot() + geom_logo(plot, method = "prob") + theme_logo() + 
  scale_y_continuous(breaks = seq(0, 1, 1), limits = c(0, 1), expand = c(0,0))+ 
  #theme(plot.title = element_text(size = 1, angle = 90))+
  facet_wrap(~seq_group,ncol = 1,scales = "fixed",shrink = F,strip.position = "right", dir = "h")

ggsave(logo, file=paste0(opt$outpu,"/motif_logo.pdf"), width=5, height = (opt$number)*0.5)
ggsave(logo, file=paste0(opt$outpu,"/motif_logo.png"), width=5, height = (opt$number)*0.5)

Dfs <- as.data.frame(Dfs)
colnames(Dfs) <- c("gene","log10pvalue","sites")
Dfs$log10pvalue <- -as.numeric(Dfs$log10pvalue)
Dfs$log10pvalue <- round(Dfs$log10pvalue, 2)
text_pvalue <- (max(Dfs$log10pvalue))* 0.1
Dfs$sites <- as.numeric(Dfs$sites)
text_sites <- (max(Dfs$sites))* 0.1
Dfs$gene <- factor(Dfs$gene, levels = rev(Dfs$gene))

pvalue <- ggplot(data = Dfs, aes(x = log10pvalue, y = gene)) +
  geom_col(aes(fill = gene)) +
  geom_text(data = Dfs, aes(x = log10pvalue + text_pvalue, label = log10pvalue)) +
  #scale_x_continuous(expand = c(0,0)) +
  xlab("-log10(pvalue)") +
  ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave(pvalue, file=paste0(opt$outpu,"/motif_pvalue.pdf"), width=5, height = (opt$number)*0.5)
ggsave(pvalue, file=paste0(opt$outpu,"/motif_pvalue.png"), width=5, height = (opt$number)*0.5)

sites <- ggplot(data = Dfs, aes(x = sites, y = gene)) +
  geom_col(aes(fill = gene)) +
  geom_text(data = Dfs, aes(x = sites + text_sites, label = sites)) +
  #scale_x_continuous(expand = c(0,0)) +
  xlab("sites") +
  ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave(sites, file=paste0(opt$outpu,"/motif_sites.pdf"), width=5, height = (opt$number)*0.5)
ggsave(sites, file=paste0(opt$outpu,"/motif_sites.png"), width=5, height = (opt$number)*0.5)


