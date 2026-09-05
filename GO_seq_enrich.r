library(tidyverse)
library(goseq)
library(GenomicFeatures)
library(msigdbr)

# gene lengths from gtf 
txdb <- makeTxDbFromGFF("~/gencode.v36.annotation.gtf.gz", format = "gtf")

txLengths <- transcriptLengths(txdb)
txLengths$gene_id <- gsub("\\..*", "", txLengths$gene_id)

medianLengths <- txLengths %>%
  group_by(gene_id) %>%
  summarise(medianLength = median(tx_len))

message("gene lengths done")

# hallmark gene sets using ensembl ids
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$ensembl_gene, f = .$gs_name)

# kegg
kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  split(x = .$ensembl_gene, f = .$gs_name)

message("gene sets loaded")

# enrichment function
run_goseq <- function(res, cancer_type, analysis, out_dir) {
  
  res <- res %>% drop_na()
  res$GeneID <- gsub("\\..*", "", rownames(res))
  
  # match gene lengths
  geneLengths <- medianLengths$medianLength[match(res$GeneID, medianLengths$gene_id)]
  
  if (sum(is.na(geneLengths)) == length(geneLengths)) {
    return(NULL)
  }
  
  prefix <- paste0(cancer_type, "_", toupper(analysis))
  
  for (direction in c("up", "down")) {
    
    sig <- as.integer(!is.na(res$padj) & res$padj < 0.05 &
                        if (direction == "up") res$log2FoldChange > 0 else res$log2FoldChange < 0)
    names(sig) <- res$GeneID
    
    n_sig <- sum(sig)
    
    if (n_sig == 0) next
    
    pwf <- tryCatch(
      nullp(sig, "hg38", "ensGene", bias.data = geneLengths, plot.fit = FALSE),
      error = function(e) { message("  nullp failed: ", e$message); NULL }
    )
    if (is.null(pwf)) next
    
    dir_prefix <- paste0(prefix, "_", toupper(direction))
    
    # hallmark
    res_hallmark <- tryCatch(
      goseq(pwf, "hg38", "ensGene", gene2cat = hallmark),
      error = function(e) { message("  hallmark failed: ", e$message); NULL }
    )
    if (!is.null(res_hallmark)) {
      res_hallmark$padj <- p.adjust(res_hallmark$over_represented_pvalue, method = "BH")
      write.csv(res_hallmark, file.path(out_dir, paste0(dir_prefix, "_HALLMARK.csv")), row.names = FALSE)
    }
    
    # kegg
    res_kegg <- tryCatch(
      goseq(pwf, "hg38", "ensGene", gene2cat = kegg),
      error = function(e) { message("  kegg failed: ", e$message); NULL }
    )
    if (!is.null(res_kegg)) {
      res_kegg$padj <- p.adjust(res_kegg$over_represented_pvalue, method = "BH")
      write.csv(res_kegg, file.path(out_dir, paste0(dir_prefix, "_KEGG.csv")), row.names = FALSE)
    }
  }
}

# loop
degs_dir <- "~/thesis/scripts/star_counts/degs/base"
out_dir  <- "~/thesis/scripts/star_counts/figures/enrichment/base"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(degs_dir, pattern = "_res_shrunk\\.rds$", full.names = TRUE)

for (f in files) {
  
  fname       <- basename(f)
  cancer_type <- gsub("_(bmi|age|gender)_res_shrunk\\.rds", "", fname)
  analysis    <- gsub(paste0(cancer_type, "_"), "", gsub("_res_shrunk\\.rds", "", fname))
  
  
  res <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(res)) next
  
  run_goseq(res, cancer_type, analysis, out_dir)
}

message("done")
