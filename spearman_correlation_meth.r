library(tidyverse)
library(data.table)
library(maxprobes)
library(minfi)
library(matrixStats)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)

TOP_GENES <- 5000  # number of most variable genes to use — set to Inf for all
# annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# functional probes with gene annotation
functional_probes <- rownames(ann)[
  ann$UCSC_RefGene_Name != "" &
    !is.na(ann$UCSC_RefGene_Name) &
    ann$UCSC_RefGene_Group %in% c("TSS200", "TSS1500", "1stExon", "5'UTR", "Body", "3'UTR")
]
message("Functional probes in annotation: ", length(functional_probes))

# clinical
clinical_dataset <- readRDS("TCGA_clinical/clinical_dataset.rds")
clinical_dataset$BMI <- clinical_dataset$weight / (clinical_dataset$height / 100)^2
clinical_dataset <- clinical_dataset %>%
  distinct(sample, .keep_all = TRUE)

# probe map
probe_map <- fread("methlylation_data/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy")

# gene sets
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$ensembl_gene, f = .$gs_name)
message(length(hallmark), " Hallmark gene sets")

# SNP probe map
snps <- getSnpInfo(IlluminaHumanMethylation450kanno.ilmn12.hg19)
snps_to_remove <- rownames(snps)[
  (!is.na(snps$CpG_rs) & snps$CpG_maf > 0.01) |
    (!is.na(snps$SBE_rs) & snps$SBE_maf > 0.01)
]

# cross reactive probes
xloci <- maxprobes::xreactive_probes(array_type = "450K")

# output
out_dir <- "figures/spearman_methyl"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# fast vectorized Spearman
fast_spearman <- function(mat, var) {
  var_ranked <- rank(var)
  mat_ranked <- t(apply(mat, 1, rank))
  n          <- length(var)
  rho        <- cor(t(mat_ranked), var_ranked)[, 1]
  tstat      <- rho * sqrt((n - 2) / (1 - rho^2))
  pval       <- 2 * pt(abs(tstat), df = n - 2, lower.tail = FALSE)
  data.frame(
    gene = rownames(mat),
    rho  = rho,
    pval = pval,
    padj = p.adjust(pval, method = "BH")
  ) %>% arrange(desc(abs(rho)))
}

# spearman + fgsea 
run_spearman_fgsea <- function(gene_mat, clin_var, var_name, out_prefix) {
  
  keep <- !is.na(clin_var)
  mat  <- gene_mat[, keep]
  var  <- clin_var[keep]
  
  if (sum(keep) < 10) {
    message("  Skipping ", var_name, ": too few samples")
    return(NULL)
  }
    
  cor_df <- fast_spearman(mat, var)
  
  message("  significant genes : ", sum(cor_df$padj < 0.05, na.rm = TRUE))
  
  write.csv(cor_df,
            paste0(out_prefix, "_", var_name, "_spearman.csv"),
            row.names = FALSE)
  
  # rank by rho for fgsea
  rankings <- cor_df$rho
  names(rankings) <- cor_df$gene
  rankings <- sort(rankings, decreasing = TRUE)
  
  message("  Running fgsea Hallmark for ", var_name, "...")
  gseares <- tryCatch({
    fgsea(
      pathways  = hallmark,
      stats     = rankings,
      scoreType = "std",
      minSize   = 10,
      maxSize   = 500,
      nproc     = 1
    )
  }, error = function(e) {
    message("  fgsea error: ", e$message)
    NULL
  })
  
  if (!is.null(gseares)) {
    gseares_df <- as.data.frame(gseares)
    saveRDS(gseares_df,
            paste0(out_prefix, "_", var_name, "_fgsea_hallmark.rds"))
    message("  Significant Hallmark pathways (padj<0.05): ",
            sum(gseares_df$padj < 0.05, na.rm = TRUE))
  }
  
  return(cor_df)
}

data_dir <- "data/methlylation_data"
files    <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
files    <- files[!grepl("probeMap", files)]

for (input_file in files) {
  
  cancer_type <- gsub("TCGA-|\\.rds", "", basename(input_file))
  message("\n", strrep("=", 60))
  message("Processing: ", cancer_type, " | Top genes: ", TOP_GENES)
  message(strrep("=", 60))
  
  df <- tryCatch(readRDS(input_file), error = function(e) {
    message("Failed to load: ", input_file)
    return(NULL)
  })
  if (is.null(df)) next
  
  # QC
  df <- as.data.frame(df)
  df <- df[!duplicated(df$sample), ]
  rownames(df) <- df$sample
  df <- df[, -1]
  
  clin <- clinical_dataset %>% filter(grepl("-01", sample))
  
  common_samples <- intersect(colnames(df), clin$sample)
  if (length(common_samples) < 10) {
    message("Skipping: too few samples (", length(common_samples), ")")
    next
  }
  
  df_filtered <- df[, common_samples]
  clin        <- clin[match(common_samples, clin$sample), ]
  
  # standard probe QC
  df_filtered <- df_filtered[!(rownames(df_filtered) %in% xloci), ]
  df_filtered <- df_filtered[!(rownames(df_filtered) %in% snps_to_remove), ]
  
  autosomal_probes <- probe_map$`#id`[!(probe_map$chrom %in% c("chrX", "chrY", ""))]
  df_filtered      <- df_filtered[rownames(df_filtered) %in% autosomal_probes, ]
  
  df_filtered <- df_filtered[rowSums(is.na(df_filtered)) < 0.1 * ncol(df_filtered), ]
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  df_filtered <- log2(as.matrix(df_filtered) / (1 - as.matrix(df_filtered)))
  df_filtered <- df_filtered[is.finite(rowSums(df_filtered)), ]
  
  # keep only functional probes with gene annotation
  df_filtered <- df_filtered[rownames(df_filtered) %in% functional_probes, ]
  
  message("  Samples: ", ncol(df_filtered), " | Functional probes: ", nrow(df_filtered))
  
  # map probes to genes and average 
  probe_gene <- data.frame(
    probe = rownames(ann),
    gene  = ann$UCSC_RefGene_Name,
    stringsAsFactors = FALSE
  ) %>%
    filter(probe %in% rownames(df_filtered)) %>%
    filter(gene != "") %>%
    mutate(gene = strsplit(gene, ";")) %>%
    unnest(gene) %>%
    distinct(probe, gene)
  
  gene_matrix <- df_filtered[probe_gene$probe, ] %>%
    as.data.frame() %>%
    mutate(gene = probe_gene$gene[match(rownames(.), probe_gene$probe)]) %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
    filter(!is.na(gene)) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  message("  Genes after averaging: ", nrow(gene_matrix))
  
  # convert gene symbols to ensembl_ids
  symbol_to_ensembl <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = rownames(gene_matrix),
    column    = "ENSEMBL",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
  
  keep_ens        <- !is.na(symbol_to_ensembl)
  gene_matrix_ens <- gene_matrix[keep_ens, ]
  rownames(gene_matrix_ens) <- symbol_to_ensembl[keep_ens]
  
  # remove duplicate ensemble_ids
  gene_vars       <- rowVars(gene_matrix_ens)
  gene_matrix_ens <- gene_matrix_ens[order(gene_vars, decreasing = TRUE), ]
  gene_matrix_ens <- gene_matrix_ens[!duplicated(rownames(gene_matrix_ens)), ]
  
  message("  Genes after Ensembl mapping and dedup: ", nrow(gene_matrix_ens))
  
  # select top N most variable genes
  if (is.finite(TOP_GENES) && nrow(gene_matrix_ens) > TOP_GENES) {
    gene_vars_final <- rowVars(gene_matrix_ens)
    top_idx         <- order(gene_vars_final, decreasing = TRUE)[1:TOP_GENES]
    gene_matrix_ens <- gene_matrix_ens[top_idx, ]
  }
  
  out_prefix <- file.path(out_dir, cancer_type)
  
  # bmi
  if (sum(!is.na(clin$BMI)) > 10) {
    run_spearman_fgsea(
      gene_mat   = gene_matrix_ens,
      clin_var   = clin$BMI,
      var_name   = "BMI",
      out_prefix = out_prefix
    )
  } else {
    message("  skipping bmi")
  }
  
  # age
  if (sum(!is.na(clin$age_at_initial_pathologic_diagnosis)) > 10) {
    run_spearman_fgsea(
      gene_mat   = gene_matrix_ens,
      clin_var   = clin$age_at_initial_pathologic_diagnosis,
      var_name   = "AGE",
      out_prefix = out_prefix
    )
  } else {
    message("  skipping age")
  }
  
  rm(df, df_filtered, clin, gene_matrix, gene_matrix_ens, probe_gene)
  gc()
}

message("\done!")
