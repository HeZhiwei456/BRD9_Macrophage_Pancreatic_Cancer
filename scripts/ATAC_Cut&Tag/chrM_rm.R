library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(getopt)
#print("Usage: Rscript chrM_report.R [-D Dir]")
options <- c("Dir", "D", 2, "character", "log file dir")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

data <- read.csv(paste0(opt$Dir,"/","chrM.txt"), sep = "\t", head = T, check.names = F)
data$chrM <- data$All - data$chrom
data <- data[, c("sample","chrom", "chrM")] %>% as.data.frame()
data <- melt(data, id="sample")
data <- ddply(data,'sample',transform, percent = value/sum(value)*100)

df = NULL
for (hist in unique(data$sample)){
  SMP <- subset(data, sample == hist)
  percent <- SMP$percent
  y = cumsum(percent)-percent/2
  SMP$pos <- y
  df <- rbind(df,SMP)
}

df$pos <- as.numeric(df$pos)
df$percent <- round(as.numeric(df$percent),2)
df$variable <- factor(df$variable, level=rev(c("chrom","chrM")))
df$sample <- factor(df$sample, level=unique(df$sample))

p <- ggplot(data = df, aes(x = sample, y = percent, fill = variable)) +
  geom_bar(position = "stack",width = 0.7, lwd = 1, stat="identity", color="black") +
  scale_fill_manual(values=c("#20B2AA", "#9370DB")) +
  geom_text(data = df,aes(x = sample, y = pos),label = df$percent, size = 5)+
  xlab("") +
  ylab("Ratio %") +coord_cartesian(ylim=c(0,100))+
  theme_bw() +	
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color="black",linewidth = 1, linetype="solid"),
        plot.title = element_blank(),
        panel.background = element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=15))
ggsave(p, file=paste0(opt$Dir,"/chrM_report.pdf"))
ggsave(p, file=paste0(opt$Dir,"/chrM_report.png"))


