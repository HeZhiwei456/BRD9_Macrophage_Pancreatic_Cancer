library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(getopt)
#print("Usage: Rscript picard_report.R [-D Dir]")
options <- c("dir", "d", 2, "character", "file dir",
             "suffix", "s", 2, "character", "file suffix",
             "output", "o", 2, "character", "output dir")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
#opt$dir = "/home/zm/Data/CutTag/ZZY/peak"
#opt$suffix = "narrowPeak"
file <- dir(opt$dir)
file <- file[grep(opt$suffix, file)]

name <- str_split_fixed(file,"_peaks", n=2)
sampleList <- unique(name[,1])
Number = NULL
for (hist in sampleList){
  data <- read.csv(paste0(opt$dir,"/",hist,"_peaks", ".", opt$suffix), sep = "\t",head =F)
  num <- nrow(data)
  Number <- rbind(Number, c(hist, num))
}
Number <- as.data.frame(Number)
colnames(Number) <- c("sample","peak_numbers")
Number$peak_numbers <- as.numeric(Number$peak_numbers)
head(Number)

p <- ggplot(Number, aes(x =sample, y = peak_numbers))+
  geom_col(position="dodge",width=0.8, fill="#D55E00", color = "black")+
  #scale_y_continuous(expand = c(0,0))+ 
  theme_bw()+ 
  theme_classic()+
  xlab("") + ylab("Peak Number") +
  theme(axis.text = element_text(color="black",size=12),
        axis.title.y = element_text(color = "black",size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(color="black",size=10, angle = 60, hjust=1))+
  geom_text(aes(label = peak_numbers), vjust = -1, hjust=0.5, colour = "black",position=position_dodge(1))

ggsave(p, file=paste0(opt$output,"/peak_numbers.pdf"), device = cairo_pdf)
ggsave(p, file=paste0(opt$output,"/peak_numbers.png"))


