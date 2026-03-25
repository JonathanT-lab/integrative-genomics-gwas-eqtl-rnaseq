library(data.table)
library(dplyr)
library(coloc)
library(tidyr)
library(ggplot2)

dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

gwas <- fread(file.path("data", "processed", "gwas_cleaned.tsv"))
eqtl <- fread(file.path("data", "raw", "GTEx_liver_eqtl.csv"))
lead_table <- fread(file.path("results", "tables", "independent_loci.tsv"))

lead <- lead_table[1, ]
lead_chr <- lead$CHR
lead_pos <- lead$BP
window <- 500000

locus_gwas <- gwas %>%
  filter(CHR == lead_chr,
         BP >= (lead_pos - window),
         BP <= (lead_pos + window))

eqtl <- as.data.frame(eqtl) %>%
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

locus_gwas <- as.data.table(locus_gwas)
locus_eqtl <- as.data.table(locus_eqtl)

gwas_coloc <- locus_gwas[, .(
  SNP,
  beta = BETA,
  varbeta = SE^2,
  MAF = EAF,
  pvalues = P
)]

eqtl_coloc <- locus_eqtl[, .(
  SNP,
  beta = NES,
  varbeta = (NES / qnorm(1 - eQTL_P / 2))^2,
  MAF = 0.1,
  pvalues = eQTL_P
)]

eqtl_coloc <- eqtl_coloc[order(pvalues)]
eqtl_coloc <- eqtl_coloc[!duplicated(SNP)]

common_snps <- intersect(gwas_coloc$SNP, eqtl_coloc$SNP)

gwas_coloc <- gwas_coloc[gwas_coloc$SNP %in% common_snps, ]
eqtl_coloc <- eqtl_coloc[eqtl_coloc$SNP %in% common_snps, ]

N_total <- 6540 + 2096759
s_prop <- 6540 / N_total

coloc_res <- coloc.abf(
  dataset1 = list(
    beta = gwas_coloc$beta,
    varbeta = gwas_coloc$varbeta,
    type = "cc",
    snp = gwas_coloc$SNP,
    MAF = gwas_coloc$MAF,
    N = N_total,
    s = s_prop
  ),
  dataset2 = list(
    beta = eqtl_coloc$beta,
    varbeta = eqtl_coloc$varbeta,
    type = "quant",
    snp = eqtl_coloc$SNP,
    MAF = eqtl_coloc$MAF,
    N = 208
  )
)

pp_df <- as.data.frame(t(coloc_res$summary[1:5]))
colnames(pp_df) <- c("PP0", "PP1", "PP2", "PP3", "PP4")
pp_df_long <- pivot_longer(pp_df, everything(), names_to = "Hypothesis", values_to = "Posterior")

write.table(pp_df_long, file.path("results", "tables", "coloc_posterior_probabilities.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p <- ggplot(pp_df_long, aes(x = Hypothesis, y = Posterior)) +
  geom_bar(stat = "identity") +
  ylim(0, 1) +
  labs(
    title = "Colocalisation Posterior Probabilities",
    y = "Posterior Probability"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path("results", "figures", "Figure4_Colocalisation_PP_Plot.png"),
  plot = p, width = 8, height = 5, dpi = 150
)

print(coloc_res$summary)

