# ORFeome – Analyzing ORFeome screening data
# Copyright (C) 2025 Anushka Shome and Franck Dumetz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see https://www.gnu.org/licenses/.

# This script identifies differentially expressed genes from ORFeome screening data and saves fold-change results as CSV files.

library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
fc_input <- as.integer(args[1])

counts <- read.csv("counts.csv", header = TRUE, row.names = 1, check.names = FALSE)
gene_desc <- counts$Description
names(gene_desc) <- rownames(counts)
counts <- as.matrix(counts[, !(colnames(counts) %in% c("Description"))])

coldata <- read.csv("treatments.csv", row.names = 1, stringsAsFactors = FALSE)
coldata$Condition <- factor(coldata$Condition)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

stopifnot(all(rownames(coldata) == colnames(counts)))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~Condition)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)

replicate_groups <- split(colnames(norm_counts), coldata$Condition)

get_fc_tables <- function(treated, untreated, norm_counts, dds, replicates) {
  contrast <- c("Condition", treated, untreated)
  res <- results(dds, contrast = contrast)
  
  fc <- data.frame()
  
  res <- res[!is.na(res$padj) & res$padj < 0.05, ]
  
  for (gene in rownames(res)) {
    untreated_vals <- norm_counts[gene, replicates[[untreated]]]
    treated_vals   <- norm_counts[gene, replicates[[treated]]]
    
    valid <- untreated_vals != 0
    
    fold_changes <- treated_vals[valid] / untreated_vals[valid]
    
    #fold_changes <- treated_vals / untreated_vals
    
    name_and_func <- strsplit(gene, ";")[[1]]
    nf <- strsplit(name_and_func, "=")
    named <- setNames(sapply(nf, `[`, 2), sapply(nf, `[`, 1))
    desc <- gene_desc[gene]
    gene_info <- data.frame(
      Gene = gene,
      Description = desc,
      Phenotype = "overrepresented",
      stringsAsFactors = FALSE
    )
    
    treated_cols <- setNames(as.list(treated_vals), paste0("Treated", seq_along(treated_vals)))
    untreated_cols <- setNames(as.list(untreated_vals), paste0("Untreated", seq_along(untreated_vals)))
    
    rep_names <- unlist(mapply(
      function(t, u) c(t, u),
      names(treated_cols),
      names(untreated_cols),
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    ))
    rep_values <- unlist(mapply(
      function(t, u) list(treated_cols[[t]], untreated_cols[[u]]),
      names(treated_cols),
      names(untreated_cols),
      SIMPLIFY = FALSE
    ))
    
    #rep_df <- as.data.frame(setNames(rep_values, rep_names))
    rep_df <- as.data.frame(as.list(rep_values))
    names(rep_df) <- rep_names
    
    if(all(valid)) {
      stats <- data.frame(
        TreatedAvg = mean(treated_vals[valid]),
        UntreatedAvg = mean(untreated_vals[valid]),
        FoldChange = mean(fold_changes),
        padj = res[gene, "padj"],
        stringsAsFactors = FALSE
      )
    } else {
      stats <- data.frame(
        TreatedAvg = mean(treated_vals[valid]),
        UntreatedAvg = mean(untreated_vals[valid]),
        FoldChange = mean(fold_changes),
        padj = res[gene, "padj"],
        stringsAsFactors = FALSE
      )
    }
    gene_info <- cbind(gene_info, rep_df, stats)
    
    if (max(untreated_vals) < 5) next
    if (min(untreated_vals) == 0) next
    
    if (all(fold_changes > fc_input, na.rm = TRUE)) fc <- rbind(fc, gene_info)
  }
  
  return(list(fc=fc))
}

no_counts <- function(untreated, replicates){
  zero <- data.frame(Gene = character())
  untreated <- unlist(replicates[untreated])
  for (gene in rownames(counts)) {
    untreated_vals <- counts[gene, untreated]
    
    name_and_func <- strsplit(gene, ";")[[1]]
    nf <- strsplit(name_and_func, "=")
    named <- setNames(sapply(nf, `[`, 2), sapply(nf, `[`, 1))
    gene_info <- data.frame(
      Gene = gene,
      Description = gene_desc[gene],
      stringsAsFactors = FALSE
    )
    
    if (max(untreated_vals) == 0) zero <- rbind(zero, gene_info)
  }
  return(zero)
}

treated_m1    <- grep("^treated.*_m1$", names(replicate_groups), value = TRUE)
untreated_m1  <- grep("^untreated.*_m1$", names(replicate_groups), value = TRUE)
treated_m10   <- grep("^treated.*_m10$", names(replicate_groups), value = TRUE)
untreated_m10 <- grep("^untreated.*_m10$", names(replicate_groups), value = TRUE)
untreated_both <- grep("^untreated", names(replicate_groups), value = TRUE)
treated_both <- grep("^treated", names(replicate_groups), value = TRUE )

zero <- no_counts(untreated_both, replicate_groups)
write.csv(zero, "results/no-counts-untreated.csv", row.names = FALSE)
zero_t <- no_counts(treated_both, replicate_groups)
write.csv(zero_t, "results/no-counts-treated.csv", row.names = FALSE)
common <- merge(zero, zero_t, by = "Gene")
common$Description.y <- NULL
names(common)[names(common) == "Description.x"] <- "Description"

write.csv(common, "results/no-counts-both.csv", row.names = FALSE)


fc_list <- list()

for (treated in c(treated_m1, treated_m10)) {
  untreated_group <- if (grepl("_m1$", treated)) untreated_m1 else untreated_m10
  for (untreated in untreated_group) {
    treated_name <- gsub("^treated_", "", treated)
    untreated_name <- gsub("^untreated_", "", untreated)
    name <- paste(treated_name, "vs", untreated_name, sep = "_")
    cat("Running:", name, "\n")
    
    fc_tables <- get_fc_tables(treated, untreated, norm_counts, dds, replicate_groups)
    
    fc_list[[name]] <- fc_tables$fc
    
    if (nrow(fc_tables$fc) > 0) {
      write.csv(fc_tables$fc, paste0("results/foldchange_", fc_input, "_", name, ".csv"), row.names = FALSE)
    }
  }
}

# Save common genes to CSV
if (length(fc_list) > 0) {
  common_genes <- Reduce(intersect, lapply(fc_list, function(x) x$Gene))
  if (length(common_genes) > 0) {
    fc_all <- do.call(rbind, lapply(fc_list, as.data.frame))
    common_df <- fc_all[fc_all$Gene %in% common_genes, ]
    common_df <- unique(common_df[, c("Gene", "Phenotype")])
    write.csv(common_df, paste0("results/foldchange_", fc_input, "_CommonHits.csv"), row.names = FALSE)
  }
}
