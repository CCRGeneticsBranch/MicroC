suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

plotMultiGG <- function(ps, filename, nrow=3, ncol=3) {
  nplots <- nrow * ncol
  npages <- length(ps) %/% (nplots)
  if (length(ps) %% (nplots) > 0 )
    npages <- npages + 1
  pdf(filename)
  for (i in c(1:npages)) {
    sidx <- (nplots*(i-1)+1)
    eidx <- nplots*i
    if (eidx > length(ps))
      eidx <- length(ps)
    #print(sidx)
    #print(eidx)
    do.call(grid.arrange, c(ps[sidx:eidx], ncol=ncol, nrow=nrow))
    #grid.newpage()
  }  
  dev.off()
}

args <- commandArgs(trailingOnly=T)
treatment_ev_file <- args[1]
ctrl_ev_file <- args[2]
out_file <- args[3]

#setwd("Y:/projects/pipeline_dev_hc/dovetail_tools/out_yong2")
#treatment_ev_file <- 'RH4_PFI90_24h_merged/compartment/ev_corrected.txt'
#ctrl_ev_file <- 'RH4_DMSO_24h_merged/compartment/ev_corrected.txt'
treatment_ev <- read.table(treatment_ev_file, head=F, sep="\t")
colnames(treatment_ev) <- c("chromosome", "start", "end", "compartment", "loading", "dot")
ctrl_ev <- read.table(ctrl_ev_file, head=F, sep="\t")
colnames(ctrl_ev) <- c("chromosome", "start", "end", "compartment", "loading", "dot")

chrs <- unique(ctrl_ev$chromosome)
joined <- ctrl_ev %>% dplyr::inner_join(treatment_ev, by=c("chromosome"="chromosome", "start"="start","end"="end"))
colnames(joined)[5] <- "DMSO"
colnames(joined)[8] <- "PFI90"
ps <- list()
for (c in chrs) {
  chr_ev <- joined %>% dplyr::filter(chromosome==c)
  chr_ev$start <- chr_ev$start/1000000
  p <- ggplot(chr_ev, aes(x=start)) + 
    geom_line(aes(y = DMSO, color = "DMSO"), size=1.2) + 
    geom_line(aes(y = PFI90, color="PFI90"), size=1.2) + xlab(paste0(c,"(MB)")) + ylab("Loading")+ 
    scale_colour_manual("", values = c("DMSO"="red","PFI90"="blue")) +
    theme_light()+ theme(legend.position="right")
  ps[[length(ps)+1]] <- p
}
plotMultiGG(ps, out_file,3,1)



ggplot(chr_ev, aes(x=start)) + 
geom_line(aes(y = DMSO, color = "DMSO"), size=1.2) + 
geom_line(aes(y = PFI90, color="PFI90"), size=1.2) + xlab(c) + ylab("Loading")+ 
  scale_colour_manual("", values = c("DMSO"="red","PFI90"="blue")) +
  theme_light()+ theme(legend.position="right")