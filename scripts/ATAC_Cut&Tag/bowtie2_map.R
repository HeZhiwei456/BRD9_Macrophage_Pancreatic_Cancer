library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(getopt)
#print("Usage: Rscript bowtie2.R [-D Dir]")
options <- c("Dir", "D", 2, "character", "log file dir",
             "type","t", 2, "character", "single or pair: S or P")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
path = opt$Dir

#path="/home/zm/Analysis/RNA-seq/STAR"
file <- dir(path)
file <- file[grep("log", file)]
name <- str_split_fixed(file,"[.]", n=2)
sampleList <- unique(name[,1])

ratio = NULL
for (hist in sampleList){
  if(opt$type == "S"){
    data <- read.csv(paste0(path,"/",hist,".log"), skip = 2, sep = " ",head =F)
    data <- as.numeric(data$V5)[1:3]
    data <- round(data/sum(data),4)*100
    data <- c(hist, data)
    ratio <- rbind(ratio, data)
  } else if(opt$type == "P"){
    data <- read.csv(paste0(path,"/",hist,".log"), skip = 2, sep = " ",head =F)[c(1:3,6),]
    map1 <- as.numeric(data[1:3,]$V5)
    map2 <- as.numeric(data[4,]$V7)
    data <- c(map1[1] - map2, map1[2], map1[3] + map2)
    data <- round(data/sum(data),4)*100
    data <- c(hist, data)
    ratio <- rbind(ratio, data)
  }
}
ratio <- as.data.frame(ratio)
colnames(ratio) <- c("sample","Unmapped reads","Uniquely mapped reads","multiple mapped reads")
ratio <- ratio[,c("sample","Uniquely mapped reads","multiple mapped reads","Unmapped reads")]
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
df$variable <- factor(df$variable, level=c("Unmapped reads","multiple mapped reads","Uniquely mapped reads"))
df$sample <- factor(df$sample, level=unique(df$sample))
df$value <- as.numeric(df$value)
df$lab <- as.numeric(df$lab)

p <- ggplot(data = df, aes(x =sample, y = value, fill = variable)) +
  geom_bar(position = "stack",width = 0.7, lwd = 1, stat="identity", color="black") +
  geom_text(data = df,aes(y= lab, x= sample),label=df$value, size = 5)+
  xlab("") +
  ylab("bowtie2 mapping ratio %") +coord_cartesian(ylim=c(0,100))+
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
        legend.text = element_text(size=14))

ggsave(p, file=paste0(path,"/bowtie2_report.pdf"), device = cairo_pdf)
ggsave(p, file=paste0(path,"/bowtie2_report.png"))
