library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)

dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path("data", "processed", "tcga_counts_matrix.tsv")
metadata_path <- file.path("data", "processed", "tcga_coldata.tsv")

counts_matrix <- as.matrix(read.table(counts_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
colData_df <- read.table(metadata_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

colData_df$condition <- factor(colData_df$condition, levels = c("Normal", "Tumor"))
rownames(counts_matrix) <- sub("\\..*$", "", rownames(counts_matrix))

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData_df,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res <- results(dds)
res$padj <- p.adjust(res$pvalue, method = "BH")

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

write.table(res_df, file.path("results", "tables", "tcga_deseq2_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

target_gene <- rownames(res_df)[1]

png(file.path("results", "figures", "Figure5_TCGA_Volcano_PNPLA3.png"), width = 1800, height = 1200, res = 150)
EnhancedVolcano(
  res_df,
  lab = res_df$gene,
  x = "log2FoldChange",
  y = "padj",
  selectLab = target_gene,
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.5,
  labSize = 4.5,
  title = "TCGA-LIHC Differential Expression",
  subtitle = "Tumor vs Normal Liver Tissue",
  caption = "PNPLA3 highlighted",
  colAlpha = 0.9
)
dev.off()

saveRDS(dds, file.path("data", "processed", "dds_object.rds"))
