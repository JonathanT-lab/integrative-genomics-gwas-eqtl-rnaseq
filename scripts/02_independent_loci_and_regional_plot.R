library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

gwas <- fread(file.path("data", "processed", "gwas_cleaned.tsv"))
eqtl <- fread(file.path("data", "raw", "GTEx_liver_eqtl.csv"))

sig <- gwas[P < 5e-8][order(P)]
window <- 500000

independent_snps <- data.frame()

while (nrow(sig) > 0) {
  lead_snp <- sig[1, ]
  independent_snps <- rbind(independent_snps, lead_snp)
  
  sig <- sig %>%
    filter(!(CHR == lead_snp$CHR &
               BP >= (lead_snp$BP - window) &
               BP <= (lead_snp$BP + window)))
}

fwrite(independent_snps, file.path("results", "tables", "independent_loci.tsv"), sep = "\t")

lead <- independent_snps[1, ]
lead_chr <- lead$CHR
lead_pos <- lead$BP

locus_gwas <- gwas %>%
  filter(CHR == lead_chr,
         BP >= (lead_pos - window),
         BP <= (lead_pos + window))

eqtl <- as.data.frame(eqtl)

eqtl <- eqtl %>%
  mutate(
    parts = stringr::str_split_fixed(`Variant Id`, "_", 5),
    CHR = as.numeric(stringr::str_remove(parts[, 1], "chr")),
    BP = as.numeric(parts[, 2])
  ) %>%
  select(-parts)

locus_eqtl <- eqtl %>%
  filter(
    CHR == lead_chr,
    BP >= (lead_pos - window),
    BP <= (lead_pos + window),
    `Gene Symbol` == "PNPLA3"
  ) %>%
  rename(
    SNP = `SNP Id`,
    eQTL_P = `P-Value`
  )

plot_df <- locus_gwas %>%
  select(SNP, BP, GWAS_P = P) %>%
  left_join(
    locus_eqtl %>% select(SNP, eQTL_P),
    by = "SNP"
  )

write.table(plot_df, file.path("results", "tables", "regional_plot_data.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p <- ggplot(plot_df, aes(x = BP)) +
  geom_point(aes(y = -log10(GWAS_P)), color = "black", alpha = 0.7) +
  geom_point(aes(y = -log10(eQTL_P)), color = "red", alpha = 0.6, na.rm = TRUE) +
  geom_vline(xintercept = lead_pos, linetype = "dashed") +
  labs(
    x = paste0("Genomic Position (chr", lead_chr, ")"),
    y = "-log10(P-value)",
    title = "Regional Association Plot: GWAS (black) vs liver eQTL (red)"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path("results", "figures", "Figure3_Regional_Association_Plot.png"),
  plot = p, width = 12, height = 5, dpi = 150
)

cat("Number of independent loci:", nrow(independent_snps), "\n")
cat("Lead SNP:", lead$SNP, "chr", lead_chr, "position", lead_pos, "\n")

