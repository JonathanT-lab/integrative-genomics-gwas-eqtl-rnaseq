library(DESeq2)
library(ggplot2)

dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(file.path("data", "processed", "dds_object.rds"))

rownames(dds) <- sub("\\..*$", "", rownames(dds))
pnpla3_id <- rownames(dds)[1]

norm_counts <- counts(dds, normalized = TRUE)
pnpla3_expr <- norm_counts[pnpla3_id, ]

pnpla3_expr_df <- data.frame(
  Expression = pnpla3_expr,
  Condition = colData(dds)$condition
)

p <- ggplot(pnpla3_expr_df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, color = "black") +
  scale_y_log10() +
  labs(
    title = "PNPLA3 Expression in TCGA-LIHC",
    x = "",
    y = "Normalized Counts (log scale)"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path("results", "figures", "Figure6_PNPLA3_Boxplot.png"),
  plot = p, width = 7, height = 5, dpi = 150
)

