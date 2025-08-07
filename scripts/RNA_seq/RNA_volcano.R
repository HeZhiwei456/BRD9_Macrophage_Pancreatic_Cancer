library(dplyr)
library(ggplot2)
library(ggrepel)
library(getopt)
library(stringr)
library(ggrastr)
options <- c("gene", "g", 1, "character", "gene file",
             "pvalue","p", 1, "numeric", "pvalue",
             "fold","f", 1, "numeric","log2fold",
             "xlim","x", 2, "numeric", "xlim",
             "ylim","y", 2, "numeric","ylim",
             "output","o", 1, "character", "output dir",
             "name","n", 1, "character", "output name")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
gene <- read.table(opt$gene, sep="\t", header = T, check.names = F)

Ctrl_name <- str_split(colnames(gene)[3], "[(|)]")[[1]][2]
Exp_name <- str_split(colnames(gene)[4], "[(|)]")[[1]][2]
xlab_name <- paste0(Exp_name, "/", Ctrl_name)

Up_regulated <- gene %>% filter(log2FoldChange > opt$fold & pvalue < opt$pvalue)
Down_regulated <- gene %>% filter(log2FoldChange < -(opt$fold) & pvalue < opt$pvalue)
Not_regulated <- gene %>% filter(abs(log2FoldChange) <= opt$fold | pvalue >= opt$pvalue)
up_number <- nrow(Up_regulated)
down_number <- nrow(Down_regulated)

label_data_up <- Up_regulated %>% top_n(10, log2FoldChange)
label_data_down <- Down_regulated %>% top_n(10, -log2FoldChange)

gene <- gene %>%
  mutate(regulation = case_when(
    log2FoldChange > opt$fold & pvalue < opt$pvalue ~ "Up_regulated",
    log2FoldChange < -(opt$fold) & pvalue < opt$pvalue ~ "Down_regulated",
    abs(log2FoldChange) <= opt$fold | pvalue >= opt$pvalue | padj == 'NA' ~ "Not_regulated"
  )) %>%
  mutate_at(vars(regulation), as.factor)
gene <- subset(gene, regulation != "NA")

#xlim <- ceiling(max(abs(gene$log2FoldChange)))
#ylim <- ceiling(max(-log10(as.numeric(gene$pvalue))))

xlim_left <- round(min(gene$log2FoldChange)*1.1, 0)
xlim_right <- round(max(gene$log2FoldChange)*1.1, 0)

logpvalue <- -log10(gene$pvalue)
logpvalue <- na.omit(logpvalue)
ylim <- round(max(logpvalue)*1.1, 0)

if (!is.null(opt$xlim)){
  xlim_left <- -opt$xlim
  xlim_right <- opt$xlim
}

if (!is.null(opt$ylim)){
  ylim <- opt$ylim
}
if(ylim > 200|ylim == Inf){ylim = 100}
#print("xlim_left: ", xlim_left)
#print("xlim_right: ", xlim_right)
#print("ylim: ", ylim)
print(ylim)
table(gene$regulation)
p <- ggplot()+
  geom_point_rast(data = gene, aes(as.numeric(log2FoldChange), -log10(as.numeric(pvalue)), color = regulation), size = 1, alpha= 0.6)+
  scale_color_manual(values = c("#57C3F3","grey", "#E95C59"))+
  geom_text_repel(data = label_data_up, aes(log2FoldChange, -log10(pvalue)), label = label_data_up$SYMBOL, color = "#625D9E", size = 4, max.overlaps = Inf)+
  geom_point(data = label_data_up, aes(log2FoldChange, -log10(pvalue)), size = 1, alpha= 0.6, color = "#625D9E")+
  
  geom_text_repel(data = label_data_down, aes(log2FoldChange, -log10(pvalue)), label = label_data_down$SYMBOL, color = "#68A180", size = 4, max.overlaps = Inf)+
  geom_point(data = label_data_down, aes(log2FoldChange, -log10(pvalue)), size = 1, alpha= 0.6, color = "#68A180")+
  
  annotate("text", x = xlim_left*0.7 , y = ylim*0.9, label=paste0("Down = ", as.character(down_number)), size= 4, color = "blue")+
  annotate("text", x = xlim_right*0.7, y = ylim*0.9, label=paste0("Up = ", as.character(up_number)), size= 4,  color = "red")+
  
  #annotate("text", x = xlim_left*0.7 , y = opt$ylim*0.9, label=paste0("Down = ", as.character(down_number)), size= 5, color = "blue")+
  #annotate("text", x = xlim_right*0.7, y = opt$ylim*0.9, label=paste0("Up = ", as.character(up_number)), size= 5,  color = "red")+
  
  geom_vline(xintercept = c(-(opt$fold), opt$fold), lty = 2, size = 0.5) +
  geom_hline(yintercept = -log10(opt$pvalue), lty = 2, size = 0.5) +
  #theme_classic(base_size = 12)+
  theme(legend.position = c(0.85,0.7),
        #legend.key.size = unit(2, 'cm'),
        axis.text = element_text(color = "black",size = 10),
        legend.text = element_text(color = "black",size = 12),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
        axis.ticks = element_line(color = "black", linewidth = 0.5))+
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim))+ 
  scale_x_continuous(expand = c(0,0), limits = c(xlim_left, xlim_right))+ 
  
  #scale_y_continuous(expand = c(0,0), limits = c(0, opt$ylime))+ 
  #scale_x_continuous(expand = c(0,0), limits = c(-opt$xlim, opt$xlim))+ 
  
  labs(x = paste0("log2", "(", xlab_name, ")"), y = bquote(-log[10](Pvalue)), color = "")
  #coord_fixed()

dir.create(opt$output, recursive = TRUE)
ggsave(p, file=paste0(opt$output,"/", opt$name, ".pdf"))
ggsave(p, file=paste0(opt$output,"/", opt$name, ".png"))
