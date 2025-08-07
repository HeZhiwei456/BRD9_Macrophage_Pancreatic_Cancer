library(DiffBind)
library(edgeR)
library(dplyr)
library(getopt)
library(ggplot2)
library(tibble)
library(tidyverse)
library(ggrepel)
library(DESeq2)
#print("Usage: Rscript reads_quality.R [-D Dir]")
options <- c("control","C", 2, "character","control",
             "experiment","E", 2, "character", "experiment",
             "control_name","c", 2, "character", "control name",
             "experiment_name", "e", 2, "character", "experiment name",
             "Input", "I", 2, "character", "input as control",
             #"input_name", "i", 2, "character", "input name",
             "bam", "b", 2, "character", "bam file dir",
             "Broad", "B", 2, "character", "Broad peak",
             "peak", "p", 2, "character", "peak file dir",
             "output","o", 2, "character","output dir",
             "Pvalue","P", 2, "numeric", "p.value",
             "method", "m", 2, "character", "diffpeak method",
             "diff", "d", 2, "character", "pvalue or FDR",
             "Fold", "F", 2, "numeric", "fold change")
spec <- matrix(options, byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

#-----------------------------------------------------
#make sample sheet file
cidr <- getwd()
dir.create(file.path(cidr, opt$output), recursive = TRUE) #, showWarnings = FALSE

ctr <- strsplit(opt$control, ",", fixed = TRUE)[[1]]
exp <- strsplit(opt$experiment, ",", fixed = TRUE)[[1]]
ctr_num <- length(ctr)
exp_num <- length(exp)

ctr_bam <- paste0(opt$bam, "/", ctr, ".bam")
exp_bam <- paste0(opt$bam, "/", exp, ".bam")
if(opt$Broad == "broad"){
  ctr_peak <- paste0(opt$peak, "/", ctr, "_peaks.broadPeak")
  exp_peak <- paste0(opt$peak, "/", exp, "_peaks.broadPeak")
  PeakType = "broad"
}else{
  ctr_peak <- paste0(opt$peak, "/", ctr, "_peaks.narrowPeak")
  exp_peak <- paste0(opt$peak, "/", exp, "_peaks.narrowPeak")
  PeakType = "narrow"
}

if (opt$Input == "NA" | opt$Input == "F") {
  df <- data.frame(SampleID=c(ctr, exp), 
                   Tissue=c(rep(NA, ctr_num+exp_num)), 
                   Factor=c(rep(opt$control_name, ctr_num),rep(opt$experiment_name, exp_num)),
                   Condition=c(rep(NA, ctr_num+exp_num)), 
                   Treatment=c(rep(NA, ctr_num+exp_num)), 
                   Replicate=c(1:ctr_num, 1:exp_num),
                   bamReads=c(ctr_bam,exp_bam), 
                   ControlID=c(rep(NA, ctr_num+exp_num)), 
                   bamControl=c(rep(NA, ctr_num+exp_num)),
                   Peaks=c(ctr_peak, exp_peak), 
                   PeakCaller=c(rep("narrow", ctr_num+exp_num)))
}else{
  input_bam <- paste0(opt$bam, "/", opt$Input, ".bam")
  df <- data.frame(SampleID=c(ctr, exp), 
                   Tissue=c(rep(NA, ctr_num+exp_num)), 
                   Factor=c(rep(opt$control_name, ctr_num),rep(opt$experiment_name, exp_num)),
                   Condition=c(rep(NA, ctr_num+exp_num)), 
                   Treatment=c(rep(NA, ctr_num+exp_num)), 
                   Replicate=c(1:ctr_num, 1:exp_num),
                   bamReads=c(ctr_bam,exp_bam), 
                   ControlID=c(rep(opt$Input, ctr_num+exp_num)), 
                   bamControl=c(rep(input_bam, ctr_num+exp_num)),
                   Peaks=c(ctr_peak, exp_peak), 
                   PeakCaller=c(rep("narrow", ctr_num+exp_num)))
}

write.table(df, file=paste0(opt$output,"/sample_sheet.csv"), sep = ",", row.names = F, quote =F)

#-----------------------------------------------------
#diffbind
#opt$output="/home/biology/Analysis/combine/20250227_昆泽_重庆/CutTag/diffpeak"
dbObj <- dba(sampleSheet=paste0(opt$output,"/sample_sheet.csv"))
dbObj

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

cairo_pdf(paste0(opt$output,"/sample_PCA.pdf"), family="DejaVu Sans")
#pdf(paste0(opt$output,"/sample_PCA.pdf"))
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

cairo_pdf(paste0(opt$output,"/sample_heatmap.pdf"), family="DejaVu Sans")
#pdf(paste0(opt$output,"/sample_heatmap.pdf"))
plot(dbObj)
dev.off()

png(paste0(opt$output,"/sample_heatmap.png"))
plot(dbObj)
dev.off()

dbObj <- dba.contrast(dbObj, categories = DBA_FACTOR, minMembers = 2)

dbObj <- dba.analyze(dbObj, method = DBA_ALL_METHODS)
dba.show(dbObj, bContrasts = T)

#opt$method = "deseq2"
if(opt$method == "deseq2"){
  dbObj_DB <- dba.report(dbObj, method = DBA_DESEQ2, contrast = 1, th=1) #DBA_DESEQ2
}else if(opt$method == "edger"){
  dbObj_DB <- dba.report(dbObj, method = DBA_EDGER, contrast = 1, th=1) #DBA_DESEQ2
}

out <- as.data.frame(dbObj_DB)
write.table(out, file = paste0(opt$output,"/All_peak_", opt$method,".xls"), sep = "\t", quote = F, row.names = F)

xlab = paste0("log2FC(", opt$experiment_name, "/", opt$control_name, ")")
if(opt$diff == "FDR"){
  diffpeak <- subset(out, FDR < opt$Pvalue & abs(Fold) > opt$Fold)
  up <- subset(diffpeak,  Fold > 0)
  down <- subset(diffpeak, Fold < 0)
  commonpeak <- subset(out, FDR >= opt$Pvalue | abs(Fold) <= opt$Fold)
  write.table(diffpeak, file = paste0(opt$output, "/diffpeak.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(up, file = paste0(opt$output, "/diffpeak_up.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(down, file = paste0(opt$output, "/diffpeak_down.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(commonpeak, file = paste0(opt$output, "/diffpeak_NoSig.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  out <- out %>% mutate(regulation = case_when(
    Fold > opt$Fold & FDR < opt$Pvalue ~ "Up_regulated",
    Fold < -(opt$Fold) & FDR < opt$Pvalue ~ "Down_regulated",
    abs(Fold) <= opt$Fold | FDR >= opt$Pvalue | FDR == 'NA' ~ "Not_regulated"
  )) %>%
    mutate_at(vars(regulation), as.factor)
  out <- subset(out, regulation != "NA")
  
  down_number <- nrow(subset(out, out$regulation == "Down_regulated"))
  up_number <- nrow(subset(out, out$regulation == "Up_regulated"))
  xlim_left <- round(min(out$Fold)*1.1, 0)
  xlim_right <- round(max(out$Fold)*1.1, 0)
  ylim <- round(max(-log10(out$FDR))*1.1, 0)
  
  p <- ggplot()+
    geom_point(data = out, aes(as.numeric(Fold), -log10(as.numeric(FDR)), color = regulation), size = 1, alpha= 0.6)+
    scale_color_manual(values = c("#57C3F3","grey", "#E95C59"))+
    
    annotate("text", x = xlim_left*0.7 , y = ylim*0.9, label=paste0("Down = ", as.character(down_number)), size= 4, color = "blue")+
    annotate("text", x = xlim_right*0.7, y = ylim*0.9, label=paste0("Up = ", as.character(up_number)), size= 4,  color = "red")+
    
    geom_vline(xintercept = c(-(opt$Fold), opt$Fold), lty = 2) +
    geom_hline(yintercept = -log10(opt$Pvalue), lty = 2) +
    theme_classic(base_size = 12)+
    theme(legend.position = "bottom",
          axis.text = element_text(color = "black",size = 10),
          legend.text = element_text(color = "black",size = 10))+
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim))+ 
    scale_x_continuous(expand = c(0,0), limits = c(xlim_left, xlim_right))+ 
    labs(x = xlab, y = "-log10(FDR)", color = "")
  #coord_fixed()
  ggsave(p, file=paste0(opt$output,"/sample_Volcano.pdf"), device = cairo_pdf)
  ggsave(p, file=paste0(opt$output,"/sample_Volcano.png"))
  
} else if(opt$diff == "pvalue"){
  diffpeak <- subset(out, p.value < opt$Pvalue & abs(Fold) > opt$Fold)
  head(diffpeak)
  up <- subset(diffpeak,  Fold > 0)
  down <- subset(diffpeak, Fold < 0)
  commonpeak <- subset(out, p.value >= opt$Pvalue | abs(Fold) <= opt$Fold)
  write.table(diffpeak, file = paste0(opt$output, "/diffpeak.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(up, file = paste0(opt$output, "/diffpeak_up.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(down, file = paste0(opt$output, "/diffpeak_down.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(commonpeak, file = paste0(opt$output, "/diffpeak_NoSig.xls"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  out <- out %>% mutate(regulation = case_when(
    Fold > opt$Fold & p.value < opt$Pvalue ~ "Up_regulated",
    Fold < -(opt$Fold) & p.value < opt$Pvalue ~ "Down_regulated",
    abs(Fold) <= opt$Fold | p.value >= opt$Pvalue | p.value == 'NA' ~ "Not_regulated"
  )) %>%
    mutate_at(vars(regulation), as.factor)
  out <- subset(out, regulation != "NA")
  
  down_number <- nrow(subset(out, out$regulation == "Down_regulated"))
  up_number <- nrow(subset(out, out$regulation == "Up_regulated"))
  xlim_left <- round(min(out$Fold)*1.1, 0)
  xlim_right <- round(max(out$Fold)*1.1, 0)
  ylim <- round(max(-log10(out$p.value))*1.1, 0)
  
  p <- ggplot()+
    geom_point(data = out, aes(as.numeric(Fold), -log10(as.numeric(p.value)), color = regulation), size = 1, alpha= 0.6)+
    scale_color_manual(values = c("#57C3F3","grey", "#E95C59"))+
    
    annotate("text", x = xlim_left*0.7 , y = ylim*0.9, label=paste0("Down = ", as.character(down_number)), size= 4, color = "blue")+
    annotate("text", x = xlim_right*0.7, y = ylim*0.9, label=paste0("Up = ", as.character(up_number)), size= 4,  color = "red")+
    
    geom_vline(xintercept = c(-(opt$Fold), opt$Fold), lty = 2) +
    geom_hline(yintercept = -log10(opt$Pvalue), lty = 2) +
    theme_classic(base_size = 12)+
    theme(legend.position = "bottom",
          axis.text = element_text(color = "black",size = 10),
          legend.text = element_text(color = "black",size = 10))+
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim))+ 
    scale_x_continuous(expand = c(0,0), limits = c(xlim_left, xlim_right))+ 
    labs(x = xlab, y = "-log10(pvalue)", color = "")
  #coord_fixed()
  ggsave(p, file=paste0(opt$output,"/sample_Volcano.pdf"), device = cairo_pdf)
  ggsave(p, file=paste0(opt$output,"/sample_Volcano.png"))
  
}

samples <- dbObj$samples$SampleID
df <- NULL
for(i in 1:length(samples)){
  rpkm <- dbObj$peaks[[i]][,c(1:3,5)]
  colnames(rpkm) <- c("seqnames", "start", "end", samples[i])
  if(is.null(df)){
    df <- rpkm
  }else{
    df <- merge(df, rpkm, by = c("seqnames", "start", "end"))
  }
}
df <- as.data.frame(df)
write.table(df, paste0(opt$output,"/sample_All_rpkm.xls"), row.names = F, col.names = T, sep = "\t", quote = F)

df_diff <- merge(diffpeak, df, by = c("seqnames", "start", "end"))
write.table(df_diff, paste0(opt$output,"/sample_diffpeak_rpkm.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
