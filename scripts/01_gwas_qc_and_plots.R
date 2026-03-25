library(data.table)
library(qqman)

dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

gwas_path <- file.path("data", "raw", "GCST90809296.tsv")
gwas <- fread(gwas_path, sep = "\t", header = TRUE)

colnames(gwas) <- c(
  "CHR", "BP", "EA", "OA", "BETA", "SE", "EAF", "P",
  "SNP", "DIR", "HETDF", "HET_PVAL", "N"
)

gwas[, CHR := as.numeric(CHR)]
gwas[, BP := as.numeric(BP)]
gwas[, P := as.numeric(P)]
gwas[, BETA := as.numeric(BETA)]
gwas[, SE := as.numeric(SE)]
gwas[, EAF := as.numeric(EAF)]

gwas <- gwas[!is.na(P)]
gwas <- gwas[P > 0 & P <= 1]
gwas <- gwas[CHR %in% 1:22]
gwas <- gwas[!is.na(BP) & !is.na(SNP)]

fwrite(gwas, file.path("data", "processed", "gwas_cleaned.tsv"), sep = "\t")

png(file.path("results", "figures", "Figure1_Manhattan_plot.png"), width = 1800, height = 900, res = 150)
manhattan(
  gwas,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5),
  col = c("steelblue4", "skyblue")
)
dev.off()

png(file.path("results", "figures", "Figure2_QQ_plot.png"), width = 1800, height = 900, res = 150)
qq(gwas$P, main = "QQ Plot")
dev.off()

cat("Total variants after cleaning:", nrow(gwas), "\n")
cat("Genome-wide significant variants:", sum(gwas$P < 5e-8, na.rm = TRUE), "\n")
