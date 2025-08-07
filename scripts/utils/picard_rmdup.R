library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(getopt)
#print("Usage: Rscript picard_report.R [-D Dir]")
options <- c("Dir", "D", 2, "character", "log file dir")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
file <- dir(opt$Dir)
file <- file[grep("_rmDup", file)]
name <- str_split_fixed(file,"_rmDup", n=2)
sampleList <- unique(name[,1])
ratio = NULL
for (hist in sampleList){
  data <- read.csv(paste0(opt$Dir,"/",hist,"_rmDup.txt"), skip = 7, sep = "\t",head =F)
  data <- round(as.numeric(data[1,9])*100, 2)
  ratio <- rbind(ratio, c(hist, 100 - data, data))
}
ratio <- as.data.frame(ratio)
#head(ratio)
colnames(ratio) <- c("sample","Uniquely reads","Duplication reads")
ratio <- melt(ratio, id ="sample")
ratio$value <- as.numeric(ratio$value)
df = NULL
lab = NULL
for (sam in unique(ratio$sample)){
  pos <- subset(ratio, sample == sam)
  value <- pos$value
  y = cumsum(value)-value/2
  pos$lab <- y
  df <- rbind(df,pos)
}

df$value <- as.numeric(df$value)
df$lab <- as.numeric(df$lab)
df$variable <- factor(df$variable, level=rev(c("Uniquely reads","Duplication reads")))
df$sample <- factor(df$sample, level=unique(df$sample))
head(df)
p <- ggplot(data = df, aes(x =sample, y = value, fill = variable)) +
  geom_bar(position = "stack",width = 0.7, lwd = 1, stat="identity", color="black") +
  geom_text(data = df,aes(y= lab, x= sample),label=df$value, size = 5)+
  xlab("") +
  ylab("duplication %") +coord_cartesian(ylim=c(0,100))+
  theme_bw() +	
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color="black",size = 1, linetype="solid"),
        plot.title = element_blank(),
        panel.background = element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=16))

ggsave(p, file=paste0(opt$Dir,"/picard_report.pdf"), device = cairo_pdf)
ggsave(p, file=paste0(opt$Dir,"/picard_report.png"))


