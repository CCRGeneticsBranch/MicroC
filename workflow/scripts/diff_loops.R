library(hictoolsr)
library(dbscan)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

run_diff_loops <- function (wt_loops, dox_loops, hicFiles, name="MicroC", res=10000, buffer=10, max_count=20) {
  out_prefix <- paste0("../../../analysis/", name, "_")
  res_dir <- paste0(res/1000, "k") 
  loops <- as_ginteractions(mergeBedpe(bedpeFiles=c(wt_loops, dox_loops), res=res, selectCol = 8))
  #colnames(loops) <- c("seqnames1", "start1", "end1","seqnames2", "start2", "end2", "value","source","cluster")
  #loops <- as_ginteractions(read.table(paste0(res_dir, "/H82N_H82NDox.merged.", res, ".bedpe")))
  
  loopCounts <- extractCounts(bedpe = loops,
                              hic = hicFiles,
                              chroms = paste0("chr", c(1:22, "X")),
                              res = res,
                              norm = 'NONE',
                              matrix = 'observed')
  #colnames(mcols(loopCounts)) <- gsub("H82N", "", gsub("H82N_", "WT_", gsub("\\.hic", "", gsub("^S._", "", colnames(mcols(loopCounts))))))
  #colnames(mcols(loopCounts)) <- gsub("H82N", "", gsub("H82N_", "WT_", gsub("\\.hic", "", gsub("hichip._", "", colnames(mcols(loopCounts))))))
  colnames(mcols(loopCounts))[(ncol(mcols(loopCounts))-3):ncol(mcols(loopCounts))] <- c("WT_Rep1","WT_Rep2","Dox_Rep1","Dox_Rep2" )
  loopsCountsDF <- as.data.frame(loopCounts)
  rn <- paste(loopsCountsDF$seqnames1,loopsCountsDF$start1,loopsCountsDF$end1,loopsCountsDF$seqnames2,loopsCountsDF$start2,loopsCountsDF$end2,sep="_")
  cnts <- 
    mcols(loopCounts)[grep("WT|Dox", colnames(mcols(loopCounts)))] |>
    as.matrix()
  rownames(cnts) <- rn
  colData <- 
    do.call(rbind, strsplit(x = colnames(cnts), split = "_")) |>
    as.data.frame(stringsAsFactors = TRUE) |>
    `colnames<-`(value = c("condition", "biorep"))
  
  dds <- 
    DESeqDataSetFromMatrix(countData = cnts,
                           colData = colData,
                           design = ~condition)
  deseq_res <-
    DESeq(dds) |>
    lfcShrink(coef = "condition_WT_vs_Dox", type="apeglm")
  
  padj_cutoff <- 0.25
  lfc_cutoff <- 0
  
  pdf(paste0(out_prefix, "Volcano.", res_dir, ".pdf"))
  p <- EnhancedVolcano(deseq_res,lab = rownames(deseq_res), ylab="-Log10 AdjP", x = 'log2FoldChange',y = 'padj', pCutoff = padj_cutoff,
                       FCcutoff = lfc_cutoff, title = name)
  print(p)
  dev.off()
  
  pdf(paste0(out_prefix, "MA.", res_dir, ".pdf"))
  plotMA(deseq_res)
  dev.off()
  
  pdf(paste0(out_prefix, "PCA.", res_dir, ".pdf"))
  p<-plotPCA(vst(dds), intgroup = "condition") + ggplot2::theme(aspect.ratio = 1)
  print(p)
  dev.off()
  
  mcols(loopCounts) <- cbind(mcols(loopCounts), deseq_res)
  loopCountsDF <- as.data.frame(loopCounts)
  
  genes_df <- read.table("../../../analysis/gencode.v39.annotation.coding_genes.bed", head=F)
  colnames(genes_df) <- c("seqnames", "start", "end", "gene_name")
  genes <- makeGRangesFromDataFrame(genes_df, ignore.strand=T)
  x <- makeGRangesFromDataFrame(loopCountsDF, 
                                seqnames.field="seqnames1",
                                start.field="start1",
                                end.field="end1",
                                strand.field="strand1")
  nearest_genes <- nearest(x, genes)
  nearest_gene_names <- genes_df$gene_name[nearest_genes]
  loopCountsDF$gene1 <- nearest_gene_names
  x <- makeGRangesFromDataFrame(loopCountsDF, 
                                seqnames.field="seqnames2",
                                start.field="start2",
                                end.field="end2",
                                strand.field="strand2")
  nearest_genes <- nearest(x, genes)
  nearest_gene_names <- genes_df$gene_name[nearest_genes]
  loopCountsDF$gene2 <- nearest_gene_names
  
  write.table(loopCountsDF, paste0(out_prefix, "DiffLoops.", res_dir, ".txt"), row.names=T, col.names=NA, sep="\t", quote=F)  
  
  wtLoops <- loopCounts[!is.na(loopCounts$padj) & loopCounts$padj <= padj_cutoff &
                          loopCounts$log2FoldChange > lfc_cutoff]
  doxLoops <- loopCounts[!is.na(loopCounts$padj) & loopCounts$padj <= padj_cutoff &
                           loopCounts$log2FoldChange < lfc_cutoff]
  
  #wtLoops <- loopCounts[loopCounts$log2FoldChange > 0]
  #doxLoops <- loopCounts[loopCounts$log2FoldChange < 0]
  
  
  wtLoops <- binBedpe(wtLoops, res, "center", "center")
  doxLoops <- binBedpe(doxLoops, res, "center", "center")
  loopList <- 
    list(allLoops = loopCounts,
         wtLoops = wtLoops,
         doxLoops = doxLoops)
  
  filteredLoops <- 
      lapply(X = loopList,
             FUN = filterBedpe,
             res = res,
             buffer = buffer) |>
      `names<-`(value = names(loopList))
   
  
  lapply(filteredLoops, summary)
  
  wtHicPath <- paste0(name, "_H82N/", name, "_H82N.hic")
  doxHicPath <- paste0(name, "_H82N_Dox/", name, "_H82N_Dox.hic")
  
  loopApaWtHic <-
    lapply(X = filteredLoops,
           FUN = calcApa,
           hic = wtHicPath,
           norm = "KR",
           buffer = buffer)
  loopApaDoxHic <-
    lapply(X = filteredLoops,
           FUN = calcApa,
           hic = doxHicPath,
           norm = "KR",
           buffer = buffer)
  nLoops <- lapply(filteredLoops, length)
  
  ## Divide each matrix by nLoops
  loopApaWtHic <- Map("/", loopApaWtHic, nLoops)
  loopApaDoxHic <- Map("/", loopApaDoxHic, nLoops)
  
  long <- 
    loopApaWtHic$allLoops |>
    as.table() |>
    as.data.frame() |>
    setNames(c('rows', 'cols', 'counts'))
  #ggplot(data = long,
  #       mapping = aes(x = rows, y = cols, fill = counts)) +
  #  scale_fill_distiller(palette = 'YlGnBu', direction = 1) + 
  #  geom_tile() +
  #  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust=1))
  
  toLong <- \(x) {
    x |>
      as.table() |>
      as.data.frame() |>
      setNames(c('rows', 'cols', 'counts'))
  }
  
  ## Apply function to convert all matrices to long format
  apas <- list(WT = lapply(loopApaWtHic, toLong),
               Dox = lapply(loopApaDoxHic, toLong))
  
  ## Add loopType to each data.frame and combine
  apas <- lapply(apas, \(x) do.call(rbind, Map(cbind, x, loopType = names(x))))
  
  ## Add hicMap to each data.frame and combine
  apas <- do.call(rbind, Map(cbind, apas, hicMap = names(apas)))
  
  ## Reorder factors
  apas$loopType <- factor(x = apas$loopType,
                          levels = c("allLoops", "wtLoops", "doxLoops"))
  apas$hicMap <- factor(x = apas$hicMap,
                        levels = c("WT", "Dox"))
  
  ## Visualize with ggplot2
  pdf(paste0(out_prefix, "APA.", res_dir, ".pdf"))
  apas$counts[which(apas$counts > max_count)] <- max_count
  p <- ggplot(data = apas,
         mapping = aes(x = rows, y = cols, fill = counts)) +
    scale_fill_distiller(palette = 'YlGnBu', direction = 1) + 
    facet_grid(hicMap~loopType, scales = "free") +
    geom_tile() +
    theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust=1))
  print(p)
  dev.off()
}
#setwd("/data/khanlab3/hsienchao/MicroC/processed_DATA/results/hic")
setwd("X:/hsienchao/MicroC/processed_DATA/results/hic")
wt_loops <- paste0(res_dir, "/H82N.merged.", res, ".bedpe")
dox_loops <-paste0(res_dir, "/H82NDox.merged.", res, ".bedpe")
hicFiles <- c("S4_H82N_Rep1.hic", "S3_H82N_Rep2.hic", "S2_H82NDox_Rep1.hic", "S1_H82NDox_Rep2.hic")
run_diff_loops(wt_loops, dox_loops, hicFiles, "MicroC", 10000, 10) 
wt_loops <- paste0(res_dir, "/hichip1_H82N.merged.", res, ".bedpe")
dox_loops <-paste0(res_dir, "/hichip2_H82NDox.merged.", res, ".bedpe")
hicFiles <- c("hichip1_H82N_Rep1.hic", "hichip1_H82N_Rep2.hic", "hichip2_H82NDox_Rep1.hic", "hichip2_H82NDox_Rep2.hic")
run_diff_loops(wt_loops, dox_loops, hicFiles, "HiChIP", 10000, 8) 

