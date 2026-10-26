# ==============================================================================
# RECOVER PROBE IDS AND RUN ANNOTATION
# ==============================================================================

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(tidyverse)

project_folder <- "F:/Malnutrition_Epigenomics_Project"
setwd(project_folder)

cat("=== RECOVERING PROBE IDS ===\n")

# Re-run limma analysis (fast since data is already in memory)
pheno_clean$Group <- factor(pheno_clean$Group, levels = c("NESAM", "ESAM"))
design <- model.matrix(~ Group + Immune_Proportion, data = pheno_clean)

fit <- lmFit(beta_clean, design)
fit_eb <- eBayes(fit)

# Get results WITH probe IDs as row names
results <- topTable(fit_eb, coef = "GroupESAM", number = Inf, adjust.method = "BH")

# CRITICAL: Convert row names to a column BEFORE saving
results$Probe_ID <- rownames(results)

# Save with probe IDs as a proper column
write.csv(results, "All_Probes_IDAT.csv", row.names = FALSE)
cat("✓ Re-saved All_Probes_IDAT.csv with probe IDs\n")
cat("✓ First 5 probe IDs:", head(results$Probe_ID, 5), "\n\n")

# Now extract significant DMPs
sig_dmps <- results %>% 
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.05)

cat("Found", nrow(sig_dmps), "significant DMPs\n\n")

# ==============================================================================
# GENE ANNOTATION
# ==============================================================================

cat("=== GETTING GENOMIC ANNOTATION ===\n")
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_df <- as.data.frame(ann)
ann_df$Probe_ID <- rownames(ann_df)

# Match probe IDs
matched_probes <- intersect(sig_dmps$Probe_ID, ann_df$Probe_ID)
cat("Found", length(matched_probes), "matching probes\n\n")

# Select correct columns and rename
ann_subset <- ann_df %>%
  filter(Probe_ID %in% matched_probes) %>%
  select(Probe_ID, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_Island) %>%
  rename(Chr = chr,
         Position = pos,
         Gene_Name = UCSC_RefGene_Name,
         Gene_Group = UCSC_RefGene_Group,
         Island_Relation = Relation_to_Island)

# Merge
annotated_dmps <- sig_dmps %>%
  left_join(ann_subset, by = "Probe_ID")

# Save results
write.csv(annotated_dmps, "Annotated_DMPs_IDAT.csv", row.names = FALSE)
cat("✓ Saved: Annotated_DMPs_IDAT.csv\n\n")

# Summary statistics
cat("=== SUMMARY STATISTICS ===\n")
cat("Total significant DMPs:", nrow(annotated_dmps), "\n")

has_gene <- sum(!is.na(annotated_dmps$Gene_Name) & annotated_dmps$Gene_Name != "", na.rm = TRUE)
cat("DMPs with gene associations:", has_gene, "\n")
cat("DMPs intergenic:", nrow(annotated_dmps) - has_gene, "\n\n")

cat("Genomic context:\n")
print(table(annotated_dmps$Gene_Group, useNA = "ifany"))

cat("\nCpG island relation:\n")
print(table(annotated_dmps$Island_Relation, useNA = "ifany"))

# Extract gene list for pathway analysis
cat("\n=== PREPARING GENE LIST ===\n")
all_genes <- unlist(strsplit(annotated_dmps$Gene_Name, split = ";"))
all_genes <- trimws(all_genes)
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
gene_list <- unique(all_genes)

cat("Unique genes identified:", length(gene_list), "\n")
if (length(gene_list) > 0) {
  cat("Top 20 genes:\n")
  print(head(gene_list, 20))
  
  writeLines(gene_list, "Gene_List_for_Pathway_Analysis.txt")
  cat("\n✓ Saved: Gene_List_for_Pathway_Analysis.txt\n")
  cat("  (Use this with g:Profiler for pathway analysis)\n")
}

# Create clean summary table
summary_table <- annotated_dmps %>%
  select(Probe_ID, logFC, adj.P.Val, Chr, Position, Gene_Name, Gene_Group, Island_Relation) %>%
  arrange(adj.P.Val)

write.csv(summary_table, "DMPs_Summary_Table.csv", row.names = FALSE)
cat("\n✓ Saved: DMPs_Summary_Table.csv\n\n")

cat("🎉 COMPLETE! Check your folder for the annotated files.\n")