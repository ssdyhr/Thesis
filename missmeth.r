# libraries
library(tidyverse)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(msigdbr)

# annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# symbols to entrez
message("loading gene sets...")

symbols_to_entrez <- function(gene_list) {
  lapply(gene_list, function(symbols) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = symbols,
      column    = "ENTREZID",
      keytype   = "SYMBOL",
      multiVals = "first"
    )
    as.character(mapped[!is.na(mapped)])
  })
}

msig_h          <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_entrez <- symbols_to_entrez(split(msig_h$gene_symbol, msig_h$gs_name))
message(length(hallmark_entrez), " hallmark sets loaded")

# enrichment function 
run_enrichment <- function(sig_cpg, all_cpg, out_prefix, out_dir) {
  
  # filter to probes in annotation
  valid_probes  <- rownames(ann)
  sig_cpg_clean <- sig_cpg[sig_cpg %in% valid_probes]
  all_cpg_clean <- all_cpg[all_cpg %in% valid_probes]
  
  message("  sig: ", length(sig_cpg), " -> ", length(sig_cpg_clean),
          " | all: ", length(all_cpg), " -> ", length(all_cpg_clean))
  
  if (length(sig_cpg_clean) == 0) {
    message("  nothing left after filtering, skipping")
    return(NULL)
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # KEGG
  message("  kegg...")
  gst_kegg <- tryCatch({
    gometh(
      sig.cpg    = sig_cpg_clean,
      all.cpg    = all_cpg_clean,
      anno       = ann,
      sig.genes  = TRUE,
      collection = "KEGG",
      array.type = "450K"
    )
  }, error = function(e) { message("  kegg failed: ", e$message); NULL })
  
  if (!is.null(gst_kegg) && nrow(gst_kegg) > 0) {
    write.csv(gst_kegg, file.path(out_dir, paste0(out_prefix, "_KEGG.csv")), row.names = TRUE)
    message("  sig kegg: ", sum(gst_kegg$FDR < 0.05, na.rm = TRUE))
  }
  
  # hallmark
  message("  hallmark...")
  gst_hallmark <- tryCatch({
    gsameth(
      sig.cpg    = sig_cpg_clean,
      all.cpg    = all_cpg_clean,
      collection = hallmark_entrez,
      anno       = ann,
      array.type = "450K"
    )
  }, error = function(e) { message("  hallmark failed: ", e$message); NULL })
  
  if (!is.null(gst_hallmark) && nrow(gst_hallmark) > 0) {
    write.csv(gst_hallmark, file.path(out_dir, paste0(out_prefix, "_HALLMARK.csv")), row.names = TRUE)
    message("  sig hallmark: ", sum(gst_hallmark$FDR < 0.05, na.rm = TRUE))
  }
}

# list cancer types
results_dir <- "/home/animaldyhr/results"
files       <- list.files(results_dir, pattern = "_methylation_results_2\\.rds$",
                          full.names = TRUE)

for (f in files) {
  
  cancer_type <- gsub("_methylation_results_2\\.rds", "", basename(f))
  message("\n--- ", cancer_type, " ---")
  
  res <- readRDS(f)
  
  for (model in c("base", "tss_only", "subtype", "subtype_tss")) {
    for (analysis in c("bmi", "age", "gender")) {
      
      key     <- paste0(analysis, "_", model)
      res_obj <- res[[key]]
      
      if (is.null(res_obj)) next
      
      all_cpg <- rownames(res_obj)
      
      # split hypo/hyper
      sig_hyper <- rownames(res_obj[
        !is.na(res_obj$adj.P.Val) &
          res_obj$adj.P.Val < 0.05  &
          res_obj$logFC > 0, ])
      
      sig_hypo <- rownames(res_obj[
        !is.na(res_obj$adj.P.Val) &
          res_obj$adj.P.Val < 0.05  &
          res_obj$logFC < 0, ])
      
      out_dir <- file.path("/home/animaldyhr/figures/enrichment", model)
      
      if (length(sig_hyper) > 0) {
        run_enrichment(sig_hyper, all_cpg,
                       paste0(cancer_type, "_", toupper(analysis), "_HYPER"), out_dir)
      }
      
      if (length(sig_hypo) > 0) {
        run_enrichment(sig_hypo, all_cpg,
                       paste0(cancer_type, "_", toupper(analysis), "_HYPO"), out_dir)
      }
    }
  }
}
