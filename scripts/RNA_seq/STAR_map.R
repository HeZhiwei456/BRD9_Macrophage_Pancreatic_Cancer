library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(getopt)
#print("Usage: Rscript STAR_map.R [-D Dir]")
options <- c("Dir", "D", 2, "character", "json file dir",
             "output", "o", 2, "character", "output dir")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
path = opt$Dir
output = opt$output
#path="/home/zm/Analysis/RNA-seq/STAR"
file <- dir(path)
file <- file[grep("final.out", file)]
file_number <- length(file)
name <- str_split_fixed(file,"Log", n=2)
sampleList <- unique(name[,1])

ratio = NULL
for (hist in sampleList){
  data <- read.csv(paste0(path,"/",hist,"Log.final.out"), skip = 9, sep = "\t",head =F)
  data <- data[c(1,16,18,21,23,25,28),]
  data <- data[,2]
  data <- str_split_fixed(data,"%",2)[,1] %>% as.numeric()
  data <- c(hist, data)
  ratio <- rbind(ratio, data)
}
ratio <- as.data.frame(ratio)
ratio
colnames(ratio) <- c("sample","Uniquely mapped reads","multiple mapped reads", "too many loci reads","too many mismatches reads", "too short unmapped reads","unmapped other reads","chimeric reads")
ratio <- melt(ratio, id ="sample")
head(ratio)
ratio$variable <- factor(ratio$variable, level=rev(unique(ratio$variable)))

p <- ggplot(data = ratio, aes(x = sample, y = as.numeric(value), fill= variable)) +
  geom_bar(position = "stack",width = 0.7, lwd = 0.5, stat="identity", color="black") +
  #scale_fill_manual(values = c("hotpink3","darkolivegreen2"))+
  xlab("sample") +
  ylab("mapping ratio %") +
  scale_y_continuous(expand = c(0,0), limits = c(0,101))+
  #theme_bw() +	
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color="black",size = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        plot.title = element_blank(),
        panel.background = element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust=1),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(color = "black", size=10))
ggsave(p, file=paste0(output,"/STAR.pdf"))
ggsave(p, file=paste0(output,"/STAR.png"))
