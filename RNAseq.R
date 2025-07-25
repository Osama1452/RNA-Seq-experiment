if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RnaSeqGeneEdgeRQL")
browseVignettes("RnaSeqGeneEdgeRQL")
library(RnaSeqGeneEdgeRQL)
targetsFile <- system.file("extdata", "targets.txt",
                           package="RnaSeqGeneEdgeRQL")
targets <- read.delim(targetsFile, stringsAsFactors=FALSE)
targets
group <- paste(targets$CellType, targets$Status, sep=".")
group <- factor(group)
table(group)


if( !file.exists("GSE60450_Lactation-GenewiseCounts.txt.gz") ) {
  FileURL <- paste(
    "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
    "format=file",
    "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
    sep="&")
  download.file(FileURL, method="libcurl", "GSE60450_Lactation-GenewiseCounts.txt.gz")
}

GenewiseCounts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz",
                             row.names="EntrezGeneID")
colnames(GenewiseCounts) <- substring(colnames(GenewiseCounts),1,7)
dim(GenewiseCounts)
head(GenewiseCounts)
GenewiseCounts

library(edgeR)
y <- DGEList(GenewiseCounts[,-1], group=group,
             genes=GenewiseCounts[,1,drop=FALSE])
options(digits=3)
y$samples

write.csv(y$samples, file = "samples_data.csv")
write.csv(y$counts, file = "gene_counts.csv")

library(org.Mm.eg.db)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENTREZID", column="SYMBOL")
head(y$genes)
y$genes
write.csv(y$genes, file = "gene_mapping.csv")

y <- y[!is.na(y$genes$Symbol), ]
dim(y)
write.csv(y[!is.na(y$genes$Symbol), ], file = "filtered_genes.csv")

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
write.csv(design, file = "design_matrix.csv")

keep <- filterByExpr(y, design)
table(keep)
write.csv(y[keep, ], file = "filtered_expression_data.csv")

y <- y[keep, , keep.lib.sizes=FALSE]

AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)

y <- calcNormFactors(y)
y$samples
write.csv(y$samples, file = "normalized_samples.csv")

pch <- c(0,1,2,15,16,17)
colors <- rep(c("gold", "red", "navy"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)


plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

## ----mdplot11, fig.width=7, fig.height=7, fig.cap="Mean-difference plot of log2-expression in sample 11 versus the average log2-expression across all other samples. The plot shows a number of genes that are both highly expressed and highly up-regulated."----
plotMD(y, column=11)
abline(h=0, col="red", lty=2, lwd=2)

par(mfrow=c(3,4)) # Arrange plots in a 3x4 grid for 12 samples

# Define sample names based on provided data
sample_names <- c("MCL1_DG", "MCL1_DH", "MCL1_DI", "MCL1_DJ", "MCL1_DK", "MCL1_DL", 
                  "MCL1_LA", "MCL1_LB", "MCL1_LC", "MCL1_LD", "MCL1_LE", "MCL1_LF")

# Loop through each of the 12 samples
for (i in 1:12) {
  # Generate MD plot for the i-th sample
  plotMD(y, column=i)
  # Add a horizontal line at y=0
  abline(h=0, col="red", lty=2, lwd=2)
  # Add a title with the actual sample name and increased size for clarity
  title(main=sample_names[i], cex.main=0.1, font.main=2)
}
dev.off()
dev.new()

y <- estimateDisp(y, design, robust=TRUE)
topTags(glmQLFTest(glmQLFit(y, design)), n=100)
write.csv(topTags(glmQLFTest(glmQLFit(y, design)), n=Inf)$table, file = "all_DE_results.csv")

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
fit$coefficients
write.csv(fit$coefficients, file = "glm_coefficients.csv")

plotQLDisp(fit)

title(main="Quarter-Root Dispersion vs Average Expression", 
      xlab="Average Log2 CPM", ylab="Quarter-Root QL Dispersion")

summary(fit$df.prior)
fit$df.prior
write.csv(data.frame(df.prior = fit$df.prior), file = "df_prior_values.csv", row.names = FALSE)


B.LvsP <- makeContrasts(B.lactating-B.pregnant, levels=design)
contrast.matrix <- makeContrasts(
  B.LvsP = B.lactating - B.pregnant,
  B.PvsV = B.pregnant - B.virgin,
  B.LvsV = B.lactating - B.virgin,
  L.LvsP = L.lactating - L.pregnant,
  L.PvsV = L.pregnant - L.virgin,
  L.LvsV = L.lactating - L.virgin,
  levels = design
)


res <- glmQLFTest(fit, contrast=B.LvsP)
results_list <- list()  
for (i in 1:ncol(contrast.matrix)) {
  contrast_name <- colnames(contrast.matrix)[i]
  res <- glmQLFTest(fit, contrast = contrast.matrix[, i])
  results_list[[contrast_name]] <- topTags(res, n = Inf)$table
  write.csv(results_list[[contrast_name]], 
            file = paste0(contrast_name, "_DE_results.csv"),
            row.names = FALSE)
}

topTags(res)
res
write.csv(topTags(res, n = Inf)$table, file = "B_LvsP_DE_results.csv", row.names = FALSE)

res_table <- topTags(res, n = Inf)$table
logFC_cutoff <- 1
FDR_cutoff <- 0.05
res_table$Significant <- with(res_table, 
                              ifelse(FDR < FDR_cutoff & abs(logFC) > logFC_cutoff, "DE", "Not_DE"))
plot(res_table$logFC, -log10(res_table$FDR),
     col = ifelse(res_table$Significant == "DE", "red", "gray"),
     pch = 19, cex = 0.6,
     xlab = "log2 Fold Change",
     ylab = "-log10(FDR)",
     main = "Volcano Plot for B.LvsP")
abline(v = c(-logFC_cutoff, logFC_cutoff), col = "blue", lty = 2)
abline(h = -log10(FDR_cutoff), col = "green", lty = 2)

res_list <- list()  
for (i in 1:ncol(contrast.matrix)) {
  contrast_name <- colnames(contrast.matrix)[i]
  res <- glmQLFTest(fit, contrast = contrast.matrix[, i])
  res_list[[contrast_name]] <- res
  result_table <- topTags(res, n = Inf)$table
  write.csv(result_table, file = paste0(contrast_name, "_DE_results.csv"), row.names = FALSE)
}

dir.create("volcano_plots", showWarnings = FALSE)
logFC_cutoff <- 1
FDR_cutoff <- 0.05

for (contrast_name in names(res_list)) {
  res_table <- topTags(res_list[[contrast_name]], n = Inf)$table
  res_table$Significant <- with(res_table,
                                ifelse(FDR < FDR_cutoff & abs(logFC) > logFC_cutoff, "DE", "Not_DE")
  )
  png(filename = paste0("volcano_plots/", contrast_name, "_volcano.png"),
      width = 800, height = 600)
  plot(res_table$logFC, -log10(res_table$FDR),
       col = ifelse(res_table$Significant == "DE", "red", "navy"),
       pch = 19, cex = 0.6,
       xlab = "log2 Fold Change",
       ylab = "-log10(FDR)",
       main = paste("Volcano Plot:", contrast_name))
  abline(v = c(-logFC_cutoff, logFC_cutoff), col = "orange", lty = 2)
  abline(h = -log10(FDR_cutoff), col = "green", lty = 2)
  
  dev.off()
}


is.de <- decideTests(res)
summary(is.de)
is.de
write.csv(is.de, file = "Differential_Significance_Matrix.csv", row.names = TRUE)

plotMD(res, status=is.de)

logFC_cutoff <- 1
FDR_cutoff <- 0.05

dir.create("md_plots_annotated", showWarnings = FALSE)
for (contrast_name in names(res_list)) {
  res_obj <- res_list[[contrast_name]]
  res_table <- topTags(res_obj, n = Inf)$table
  
  is_significant <- with(res_table, FDR < FDR_cutoff & abs(logFC) > logFC_cutoff)
  
  status <- ifelse(res_table$FDR < FDR_cutoff & res_table$logFC > logFC_cutoff, 1,
                   ifelse(res_table$FDR < FDR_cutoff & res_table$logFC < -logFC_cutoff, -1, 0))
  
  png(filename = paste0("md_plots_annotated/", contrast_name, "_MD_plot.png"),
      width = 2500, height = 1500, res = 300)
  
  plotMD(res_obj, status = status,
         xlab = "Average LogCPM", ylab = "log2 Fold Change",
         main = paste("MD Plot -", contrast_name))
  
  abline(h = c(-logFC_cutoff, logFC_cutoff), col = "blue", lty = 2)
  
  top_genes <- rownames(res_table[order(res_table$FDR), ])[1:50]  # أول 5 جينات حسب FDR
  top_coords <- res_table[top_genes, ]
  with(top_coords, text(logCPM, logFC, labels = top_genes, pos = 3, cex = 0.7, col = "purple"))
  
  dev.off()
}


tr <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
topTags(tr)
tr

lfc_threshold <- log2(1.5)
logFC_cutoff <- 1
FDR_cutoff <- 0.05

dir.create("treat_results", showWarnings = FALSE)
dir.create("treat_volcano", showWarnings = FALSE)
treat_list <- list()
for (contrast_name in colnames(contrast.matrix)) {
    tr <- glmTreat(fit, contrast = contrast.matrix[, contrast_name], lfc = lfc_threshold)
    treat_list[[contrast_name]] <- tr
    result_table <- topTags(tr, n = Inf)$table
    write.csv(result_table, file = paste0("treat_results/", contrast_name, "_Treat_DE.csv"), row.names = FALSE)
  
  result_table$Significant <- with(result_table,
                                   ifelse(FDR < FDR_cutoff & abs(logFC) > logFC_cutoff, "DE", "Not_DE")
  )
  
  png(filename = paste0("treat_volcano/", contrast_name, "_Treat_volcano.png"),
      width = 3500, height = 3000, res = 300)
  
  plot(result_table$logFC, -log10(result_table$FDR),
       col = ifelse(result_table$Significant == "DE", "red", "skyblue"),
       pch = 19, cex = 0.6,
       xlab = "log2 Fold Change",
       ylab = "-log10(FDR)",
       main = paste("Treat Volcano Plot:", contrast_name))
  
  abline(v = c(-logFC_cutoff, logFC_cutoff), col = "blue", lty = 2)
  abline(h = -log10(FDR_cutoff), col = "green", lty = 2)
  
  top_genes <- rownames(result_table[order(result_table$FDR), ])[1:50]
  top_coords <- result_table[top_genes, ]
  with(top_coords, text(logFC, -log10(FDR), labels = top_genes, pos = 3, cex = 0.7, col = "purple"))
  
  dev.off()
}

is.de <- decideTests(tr)
summary(is.de)
plotMD(tr, status=is.de)

dir.create("md_treat_plots", showWarnings = FALSE)
decide_list <- list()
for (contrast_name in names(treat_list)) {
    tr <- treat_list[[contrast_name]]
  is.de <- decideTests(tr)
  
  decide_list[[contrast_name]] <- is.de
  
  cat("\n--- Summary for", contrast_name, "---\n")
  print(summary(is.de))
  
  png(filename = paste0("md_treat_plots/", contrast_name, "_MD_Treat.png"),
      width  = 3000, height = 2600, res = 300)
  
  plotMD(tr, status = is.de,
         xlab = "Average LogCPM",
         ylab = "log2 Fold Change (Treat)",
         main = paste("MD plot (Treat):", contrast_name))
  
  abline(h = 0, col = "red", lty = 2)
  abline(h = c(-log2(1.5), log2(1.5)), col = "blue", lty = 3)
  top50 <- rownames(topTags(tr, n=50)$table)
  coords <- topTags(tr, n=5)$table
  with(coords, text(logCPM, logFC, labels = top50, pos=3, cex=0.7, col="maroon"))
  
  dev.off()
}


dir.create("treat_decide_tests", showWarnings = FALSE)
decide_list <- list()

for (contrast_name in names(treat_list)) {
  
  tr <- treat_list[[contrast_name]]
  is.de <- decideTests(tr)
  
  decide_list[[contrast_name]] <- is.de
  
  write.csv(as.data.frame(is.de), 
            file = paste0("treat_decide_tests/", contrast_name, "_decideTests.csv"),
            row.names = TRUE)
}


logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")

o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:30],]

coolmap(logCPM, margins=c(7,7), lhei=c(1,6), lwid=c(1,3))

annotation_col <- data.frame(Group = y$samples$group)
rownames(annotation_col) <- colnames(logCPM)
pheatmap(logCPM,
         scale = "row",  
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         fontsize_row = 7,
         fontsize_col = 9,
         main = "Top 30 DE Genes: B.lactating vs B.pregnant",
         border_color = NA)

logCPM <- cpm(y, prior.count = 2, log = TRUE)


rownames(logCPM) <- rownames(y$counts)
colnames(logCPM) <- paste(y$samples$group, seq_len(ncol(logCPM)), sep = "-")

head(rownames(logCPM))

top_genes_list <- lapply(res_list, function(res) {
  rownames(topTags(res, n = 30)$table)
})

all_top30 <- unique(unlist(top_genes_list))
length(all_top30)       

common_genes <- intersect(all_top30, rownames(logCPM))
cat("Number of shared genes:", length(common_genes), "\n")

if (length(common_genes) < 2) {
  stop("Number of shared genes is less than 2, provide n in topTags or check IDs")
  }
heatmap_data <- logCPM[common_genes, , drop = FALSE]

heatmap_scaled <- t(scale(t(heatmap_data)))

annotation_col <- data.frame(Group = y$samples$group)
rownames(annotation_col) <- colnames(heatmap_scaled)

library(pheatmap)
dir.create("heatmap_all_contrasts", showWarnings = FALSE)
pheatmap(
  heatmap_scaled,
  annotation_col       = annotation_col,
  color                = colorRampPalette(c("navy","white","firebrick3"))(50),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method    = "complete",
  fontsize_row         = 6,
  fontsize_col         = 8,
  main                 = "DE Genes Across All Comparisons",
  filename             = "heatmap_all_contrasts/DE_genes_all_contrasts.png",
  width                = 25,
  height               = 50,
  dpi                  = 300
)


library(heatmaply)
heatmaply(
  heatmap_scaled,
  k_row = 3,                  # clustering rows
  k_col = 5,                  # clustering cols
  show_dendrogram = c(TRUE, TRUE),
  hide_colorbar = FALSE,
  main = "Interactive DE Heatmap"
)

library(heatmaply)
my_colors <- colorRampPalette(c("navy","white","firebrick3"))(256)

heatmaply(
  heatmap_scaled,
  
  colors        = my_colors,
  
  width         = 3000,
  height        = 1500,

  fontsize_row  = 7,
  fontsize_col  = 12,
  
  showticklabels = c(TRUE, TRUE),
  
  main          = "DE Genes Across All Comparisons",
  
  column_side_colors = annotation_col
)


con <- makeContrasts(
  L.PvsL = L.pregnant - L.lactating,
  L.VvsL = L.virgin - L.lactating,
  L.VvsP = L.virgin - L.pregnant, levels=design)
con

groups <- levels(y$samples$group)
cmb <- combn(groups, 2, simplify = FALSE)
contrast_list <- lapply(cmb, function(x) {
  nm <- paste0(x[1], "_vs_", x[2])
  expr <- paste(x[1], "-", x[2])
  setNames(expr, nm)
})
contrast_list <- do.call(c, contrast_list)
contrast_list

res <- glmQLFTest(fit, contrast=con)
topTags(res)

dir.create("anovaQLF_results", showWarnings = FALSE)
for (coef in colnames(con)) {
  
  res_i <- glmQLFTest(fit, contrast = con[, coef])
  
  table_i <- topTags(res_i, n = Inf)$table
  
  write.csv(table_i,
            file      = file.path("anovaQLF_results",
                                  paste0(coef, "_QLFtest_results.csv")),
            row.names = TRUE)
}

all_results <- do.call(rbind, lapply(colnames(con), function(coef) {
  ti <- read.csv(file.path("anovaQLF_results",
                           paste0(coef, "_QLFtest_results.csv")),
                 row.names = 1)
  ti$Comparison <- coef
  ti
}))
write.csv(all_results,
          file      = "anovaQLF_results/all_contrasts_QLFtest_results.csv",
          row.names = TRUE)


con <- makeContrasts(
  (L.lactating-L.pregnant)-(B.lactating-B.pregnant), 
  levels=design)

library(limma)

# 1) Define your six within-group contrasts as named character vectors
L_con <- c(
  L.LvsP = "L.lactating - L.pregnant",
  L.PvsV = "L.pregnant  - L.virgin",
  L.LvsV = "L.lactating - L.virgin"
)

B_con <- c(
  B.LvsP = "B.lactating - B.pregnant",
  B.PvsV = "B.pregnant  - B.virgin",
  B.LvsV = "B.lactating - B.virgin"
)

# 2) Build a list of all “(L-diff) – (B-diff)” expressions
diff_list <- list()
for(lname in names(L_con)){
  for(bname in names(B_con)){
    # name each contrast e.g. "L.LvsP_vs_B.LvsP"
    name <- paste0(lname, "_vs_", bname)
    expr <- paste0("(", L_con[lname], ") - (", B_con[bname], ")")
    diff_list[[name]] <- expr
  }
}

# 3) Call makeContrasts() with your design matrix
#    (assumes you already have `design` defined)
con_dd <- do.call(
  makeContrasts,
  c(diff_list, list(levels = design))
)

# inspect
head(con_dd)
dim(con_dd)        # should be 9 columns, one per diff-of-diff

# 4) Save to CSV for later reference
dir.create("contrast_matrices", showWarnings = FALSE)
write.csv(
  con_dd,
  file       = "contrast_matrices/L_vs_B_diff_of_diffs.csv",
  row.names  = TRUE
)


res <- glmQLFTest(fit, contrast=con)
topTags(res)
write.csv(topTags(res, n = Inf), file = "results.csv")

contrast.names <- colnames(contrast.matrix)

for (name in contrast.names) {
  res <- glmQLFTest(fit, contrast = contrast.matrix[, name])
  
  tags <- topTags(res, n = Inf)
  
  write.csv(tags, file = paste0(name, "_results.csv"))
}


go <- goana(tr, species="Mm")
topGO(go, n=15)
go
write.csv(topGO(go, n = Inf), file = "go_results.csv")

dir.create("GO_analysis", showWarnings = FALSE)

go_results    <- list()
topgo_results <- list()

for (cmp in names(treat_list)) {
  tr <- treat_list[[cmp]]
  
  go <- goana(tr, species = "Mm")
  go_results[[cmp]] <- go
  
  tg <- topGO(go, n = 20)
  topgo_results[[cmp]] <- tg
  
  write.csv(as.data.frame(go),
            file = file.path("GO_analysis", paste0(cmp, "_goana.csv")),
            row.names = TRUE)
  
  write.csv(tg,
            file = file.path("GO_analysis", paste0(cmp, "_topGO20.csv")),
            row.names = FALSE)
}

lapply(topgo_results, head)
all_go   <- do.call(rbind, lapply(names(go_results), function(cmp) {
  df <- as.data.frame(go_results[[cmp]])
  df$Comparison <- cmp
  df
}))
write.csv(all_go, "GO_analysis/all_goana_combined.csv", row.names = TRUE)


library(dplyr)
library(tidyr)
library(ggplot2)

files <- list.files("GO_analysis", pattern="_topGO20.csv", full.names=TRUE)
go_list <- lapply(files, read.csv, stringsAsFactors=FALSE)
names(go_list) <- gsub("_topGO20.csv","", basename(files))

go_long <- bind_rows(lapply(names(go_list), function(cmp) {
  df <- go_list[[cmp]]
  
  df_sub <- df %>% select(Term, P.Up, P.Down)
  
  df_long <- df_sub %>%
    pivot_longer(cols = c(P.Up, P.Down),
                 names_to  = "Type",
                 values_to = "P") %>%
    group_by(Comparison = cmp, Type) %>%
    arrange(P) %>%
    slice(1:20) %>%         
    ungroup() %>%
    mutate(
      term     = factor(Term, levels = rev(unique(Term))),
      negLogP  = -log10(P)
    )
  
  df_long
}), .id = NULL)


p <- ggplot(go_long, aes(x = term, y = negLogP, fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  facet_wrap(~ Type, ncol = 1, scales = "free_y") +
  labs(
    x     = "",
    y     = "-log10(P-value)",
    title = "Top 10 GO Terms for Everyone Comparison for P.Up vs P.Down"
  ) +
  theme_minimal() +
  theme(
    legend.position   = "bottom",
    strip.text        = element_text(face = "bold"),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 8)
  )
print(p)
ggsave("GO_analysis/Top 20 GO Terms for Everyone Comparison for P.Up vs P.Down.png", plot = p,
       width = 40, height = 25, dpi = 300)


go_long <- bind_rows(lapply(names(go_list), function(cmp) {
  df <- go_list[[cmp]]
  
  df_sub <- df %>% select(Term, P.Up, P.Down)
  
  df_long <- df_sub %>%
    pivot_longer(cols = c(P.Up, P.Down),
                 names_to  = "Type",
                 values_to = "P") %>%
    group_by(Comparison = cmp, Type) %>%
    arrange(P) %>%
    slice(1:5) %>%           
    ungroup() %>%
    mutate(
      term     = factor(Term, levels = rev(unique(Term))),
      negLogP  = -log10(P)
    )
  
  df_long
}), .id = NULL)


p <- ggplot(go_long, aes(x = term, y = negLogP, fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  facet_wrap(~ Type, ncol = 1, scales = "free_y") +
  labs(
    x     = "",
    y     = "-log10(P-value)",
    title = "Top 5 GO Terms for Everyone Comparison for P.Up vs P.Down"
  ) +
  theme_minimal() +
  theme(
    legend.position   = "bottom",
    strip.text        = element_text(face = "bold"),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 8)
  )

print(p)
ggsave("GO_analysis/Top 5 GO Terms for Everyone Comparison for P.Up vs P.Down.png", plot = p,
       width = 30, height = 20, dpi = 300)



library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
tr <- treat_list[["B.LvsP"]]
res_table <- topTags(tr, n = Inf)$table

geneList <- res_table$logFC
names(geneList) <- rownames(res_table)    

geneList <- sort(geneList, decreasing = TRUE)
gsea_res <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",   
  ont           = "BP",        
  nPerm         = 1000,         
  minGSSize     = 10,           
  maxGSSize     = 500,          
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

str(gsea_res)
head(as.data.frame(gsea_res), n = 6)
dir.create("GO_analysis", showWarnings = FALSE)

write.csv(
  as.data.frame(gsea_res),
  file      = "GO_analysis/B.LvsP_gseaGO_results.csv",
  row.names = FALSE
)
ridgeplot(gsea_res, showCategory = 10) +
  ggtitle("Ridge Plot: Top 10 Enriched GO Terms (B.LvsP)")

firstID <- gsea_res@result$ID[1]
firstDesc <- gsea_res@result$Description[1]

gseaplot2(
  gsea_res,
  geneSetID = firstID,
  title     = firstDesc
)

dotplot(gsea_res, showCategory = 50) +
  ggtitle("Dotplot: Top 50 Enriched GO Terms (B.LvsP)")

gsea_list <- list()

for(cmp in names(treat_list)) {
  tr <- treat_list[[cmp]]
  df <- topTags(tr, n = Inf)$table
  gl <- sort(setNames(df$logFC, rownames(df)), decreasing = TRUE)
  
  gsea_list[[cmp]] <- gseGO(
    geneList      = gl,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    nPerm         = 1000,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  write.csv(
    as.data.frame(gsea_list[[cmp]]),
    file      = file.path("GO_analysis", paste0(cmp, "_gseaGO.csv")),
    row.names = FALSE
  )
}

lapply(gsea_list, function(x) head(as.data.frame(x), 3))


library(RColorBrewer)
cols <- brewer.pal(11, "Set3")  

plot_gsea_ridge <- function(gsea_obj, title_text) {
  ridgeplot(
    gsea_obj,
    showCategory = 10,
    fill = "p.adjust"
  ) +
    scale_fill_gradient(
      low = "skyblue",
      high = "navy"
    ) +
    ggtitle(title_text) +
    theme_minimal()
}

plot_gsea_dot <- function(gsea_obj, title_text) {
  dotplot(
    gsea_obj,
    showCategory = 20,
    color = "p.adjust"
  ) +
    scale_color_gradient(
      low = "gold",
      high = "red"
    ) +
    ggtitle(title_text) +
    theme_classic()
}

plot_gsea_detail <- function(gsea_obj, top_n = 50) {
  ids <- gsea_obj@result$ID[1:top_n]
  desc <- gsea_obj@result$Description[1:top_n]
  plots <- lapply(seq_along(ids), function(i) {
    gseaplot2(
      gsea_obj,
      geneSetID = ids[i],
      title     = desc[i]
    )
  })
  return(plots)
}


library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(RColorBrewer)
library(edgeR)

contrast_list <- c(
  "B.LvsP", "B.PvsV", "B.LvsV",
  "L.LvsP", "L.PvsV", "L.LvsV"
)

dir.create("GO_analysis", showWarnings = FALSE)

ridge_cols <- colorRampPalette(brewer.pal(9, "Blues"))(100)
dot_cols   <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

# 4.1 Colored ridgeplot (top 10)
plot_gsea_ridge <- function(gsea_obj, title_text) {
  ridgeplot(gsea_obj, showCategory = 50, fill = "p.adjust") +
    scale_fill_gradientn(colors = ridge_cols) +
    ggtitle(title_text) +
    theme_minimal()
}

# 4.2 Colored dotplot (top 20)
plot_gsea_dot <- function(gsea_obj, title_text) {
  dotplot(gsea_obj, showCategory = 50, color = "p.adjust") +
    scale_color_gradientn(colors = dot_cols) +
    ggtitle(title_text) +
    theme_classic()
}

# 4.3 Detailed GSEA plots for top N terms
plot_gsea_detail <- function(gsea_obj, top_n = 50) {
  res_df <- as.data.frame(gsea_obj)
  if (nrow(res_df) < 1) {
    warning("No enriched pathways returned; skipping detailed plots.")
    return(list())
  }
  ids  <- res_df$ID[1:min(top_n, nrow(res_df))]
  desc <- res_df$Description[1:length(ids)]
  lapply(seq_along(ids), function(i) {
    gseaplot2(
      gsea_obj,
      geneSetID = ids[i],
      title     = desc[i]
    )
  })
}

for (cmp in contrast_list) {
  
  tr <- treat_list[[cmp]]
  if (is.null(tr)) {
    warning("No object found in treat_list for contrast ", cmp)
    next
  }
  
  res_tab <- edgeR::topTags(tr, n = Inf)$table
  
  # Ensure ENTREZID rownames
  if (!all(rownames(res_tab) %in% keys(org.Mm.eg.db, keytype = "ENTREZID"))) {
    stop("rownames(res_tab) are not valid ENTREZIDs for ", cmp)
  }
  
  geneList <- sort(
    setNames(res_tab$logFC, rownames(res_tab)),
    decreasing = TRUE
  )
  
  gsea_res <- clusterProfiler::gseGO(
    geneList      = geneList,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  # Check for empty result
  gsea_df <- as.data.frame(gsea_res)
  if (nrow(gsea_df) == 0) {
    warning("GSEA returned zero pathways for ", cmp, 
            ". Try lowering minGSSize or raising pvalueCutoff.")
    next
  }
  
  ridge_plot <- plot_gsea_ridge(gsea_res, paste("Ridge Plot:", cmp))
  ggsave(
    filename = file.path("GO_analysis", paste0(cmp, "_ridge.pdf")),
    plot     = ridge_plot,
    width    = 20,
    height   = 15
  )
  
  dot_plot <- plot_gsea_dot(gsea_res, paste("Dotplot:", cmp))
  ggsave(
    filename = file.path("GO_analysis", paste0(cmp, "_dot.pdf")),
    plot     = dot_plot,
    width    = 20,
    height   = 15
  )
  
  detail_plots <- plot_gsea_detail(gsea_res, top_n = 50)
  if (length(detail_plots) > 0) {
    pdf(file.path("GO_analysis", paste0(cmp, "_detailed.pdf")),
        width = 20, height = 15)
    for (p in detail_plots) print(p)
    dev.off()
  }
  
  message("Completed GSEA for ", cmp)
}

message("All contrasts processed. Check the GO_analysis/ folder for PDF outputs.")


library(impute)
library(preprocessCore)
library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)

enableWGCNAThreads()
expr <- cpm(y, log=TRUE, prior.count=1)

varGenes <- apply(expr, 1, var)
keep     <- varGenes > quantile(varGenes, 0.5)
datExpr  <- expr[keep, ]
dim(datExpr)  # rows = genes, cols = samples

rownames(datExpr) <- rownames(expr)[keep]
colnames(datExpr) <- colnames(expr)
enableWGCNAThreads()
powers <- c(1:10, seq(12, 20, 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType  = "signed"  # أو "unsigned"
)

par(mfrow=c(1,2))
plot(
  sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit (R^2)",
  type="b", main="Power Selection Curve"
)
abline(h=0.8, col="red")

plot(
  powers, sft$fitIndices[,5],
  xlab="Soft Threshold (power)",
  ylab="Mean Connectivity",
  type="b", main="Average Connection vs. Power"
)

outdir <- "WGCNA_results"
dir.create(outdir, showWarnings = FALSE)

write.csv(
  sft$fitIndices,
  file      = file.path(outdir, "pickSoftThreshold_results.csv"),
  row.names = TRUE
)

pdf(file.path(outdir, "pickSoftThreshold_plots.pdf"), width = 7, height = 5)

plot(
  sft$fitIndices[, "Power"],
  -sign(sft$fitIndices[, "slope"]) * sft$fitIndices[, "SFT.R.sq"],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit (signed R^2)",
  type = "b",
  main = "Power Selection – Scale-Free Benchmark"
)
abline(h = 0.8, col = "red", lty = 2)

plot(
  sft$fitIndices[, "Power"],
  sft$fitIndices[, "mean.k."],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "b",
  main = "Average Connection vs. Power"
)

dev.off()

png(file.path(outdir, "scaleFreeTopology_fit.png"), width = 800, height = 600)
plot(
  sft$fitIndices[, "Power"],
  -sign(sft$fitIndices[, "slope"]) * sft$fitIndices[, "SFT.R.sq"],
  type = "b",
  xlab = "Soft Threshold (power)",
  ylab = "signed R^2",
  main = "Scale-Free Topology fit"
)
abline(h = 0.8, col = "red", lty = 2)
dev.off()

png(file.path(outdir, "meanConnectivity.png"), width = 800, height = 600)
plot(
  sft$fitIndices[, "Power"],
  sft$fitIndices[, "mean.k."],
  type = "b",
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  main = "Mean Connectivity vs Power"
)
dev.off()



options(timeout = 700)
keg <- kegga(tr, species="Mm")
topKEGG(keg, n=20, truncate=34)
outdir <- "KEGG_analysis"
dir.create(outdir, showWarnings = FALSE)

comparisons <- c(
  "B.LvsP",  # B.lactating - B.pregnant
  "B.PvsV",  # B.pregnant  - B.virgin
  "B.LvsV",  # B.lactating - B.virgin
  "L.LvsP",  # L.lactating - L.pregnant
  "L.PvsV",  # L.pregnant  - L.virgin
  "L.LvsV"   # L.lactating - L.virgin
)

kegg_results    <- list()
topkegg_results <- list()


options(timeout = 700)
library(RCurl)
curlSetOpt(timeout = 700)

library(KEGGREST)

for(cmp in comparisons) {
  tr <- treat_list[[cmp]]
  if (is.null(tr)) {
    warning("I did not find a result for ", cmp)
    next
  }
  
  keg <- tryCatch({
    kegga(tr, species = "Mm")
  }, error = function(e) {
    warning("kegga failed to compare ", cmp, ":\n", e$message)
    return(NULL)
  })
  if (is.null(keg)) next
  
  topkeg <- topKEGG(keg, n = 20, truncate = 34)
  
  write.csv(
    as.data.frame(keg),
    file      = file.path(outdir, paste0(cmp, "_kegga.csv")),
    row.names = TRUE
  )
  
  write.csv(
    topkeg,
    file      = file.path(outdir, paste0(cmp, "_topKEGG20.csv")),
    row.names = FALSE
  )
  
  message("✅ KEGG saved to ", cmp)
}



all_kegga <- do.call(rbind, lapply(names(kegg_results), function(cmp) {
  df <- as.data.frame(kegg_results[[cmp]])
  df$Comparison <- cmp
  df
}))
write.csv(
  all_kegga,
  file      = file.path(outdir, "all_kegga_combined.csv"),
  row.names = TRUE
)


all_topKEGG <- do.call(rbind, lapply(names(topkegg_results), function(cmp) {
  df <- topkegg_results[[cmp]]
  df$Comparison <- cmp
  df
}))
write.csv(
  all_topKEGG,
  file      = file.path(outdir, "all_topKEGG20_combined.csv"),
  row.names = FALSE
)

lapply(topkegg_results, head, n = 3)


library(clusterProfiler)
library(org.Mm.eg.db)
df     <- topTags(tr, n = Inf)$table
geneUp <- rownames(df)[df$logFC > 1 & df$FDR < 0.02]
geneDn <- rownames(df)[df$logFC < -1 & df$FDR < 0.02]

ekegg_up <- enrichKEGG(
  gene         = geneUp,
  organism     = "mmu",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05
)

ekegg_dn <- enrichKEGG(
  gene         = geneDn,
  organism     = "mmu",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05
)

head(ekegg_up)
ekegg_up
ekegg_dn
write.csv(as.data.frame(ekegg_up), "KEGG_enrich_up.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_dn), "KEGG_enrich_dn.csv", row.names = FALSE)

geneList <- df$logFC
names(geneList) <- rownames(df)
geneList <- sort(geneList, decreasing = TRUE)

gkegg <- gseKEGG(
  geneList      = geneList,
  organism      = "mmu",
  keyType       = "ncbi-geneid",
  minGSSize     = 10,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

head(as.data.frame(gkegg), 6)
ridgeplot(gkegg, showCategory=10) + ggtitle("GSEA KEGG: Ridge Plot")
dotplot(gkegg, showCategory=15) + ggtitle("GSEA KEGG: Dotplot")

p1 <- ridgeplot(gkegg, showCategory=10) +
  ggtitle("GSEA KEGG: Ridge Plot")
ggsave(
  filename = "GSEA_KEGG_ridge.png",
  plot     = p1,
  width    = 15, 
  height   = 10, 
  dpi      = 300
)

p2 <- dotplot(gkegg, showCategory=15) +
  ggtitle("GSEA KEGG: Dotplot")
ggsave(
  filename = "GSEA_KEGG_dot.png",
  plot     = p2,
  width    = 15, 
  height   = 10, 
  dpi      = 300
)



geneLists <- lapply(treat_list, function(tr){
  df <- topTags(tr, n = Inf)$table
  rownames(df)[df$FDR < 0.05 & abs(df$logFC) > 1]
})

cc_KEGG <- compareCluster(
  geneCluster  = geneLists,
  fun          = "enrichKEGG",
  organism     = "mmu",
  keyType      = "ncbi-geneid"
)

dotplot(cc_KEGG, showCategory=10) + ggtitle("CompareCluster KEGG")



comparisons <- c("B.LvsP","B.PvsV","B.LvsV","L.LvsP","L.PvsV","L.LvsV")
outdir <- "KEGG_ORA_plots"
dir.create(outdir, showWarnings = FALSE)

ekegg_list <- list()
geneLists  <- list()

for(cmp in comparisons){
  df  <- topTags(treat_list[[cmp]], n=Inf)$table
  sig <- subset(df, FDR<0.05 & abs(logFC)>1)
  genes <- rownames(sig)  # ENTREZIDs expected
  
  # 2.2 enrichKEGG ORA
  ekegg <- enrichKEGG(
    gene         = genes,
    organism     = "mmu",
    keyType      = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
  ekegg_list[[cmp]] <- ekegg
  geneLists[[cmp]]  <- genes
  
  # 2.3 barplot
  p_bar <- barplot(ekegg, showCategory=10, x="Count", color="p.adjust") +
    scale_fill_viridis_c(option="C", direction=-1) +
    ggtitle(paste(cmp,"KEGG Barplot")) + theme_minimal()
  ggsave(file.path(outdir,paste0(cmp,"_barplot.png")), p_bar, width=10, height=7, dpi=300)
  
  # 2.4 dotplot
  p_dot <- dotplot(ekegg, showCategory=10, color="p.adjust") +
    scale_color_gradient(low="gold", high="red") +
    ggtitle(paste(cmp,"KEGG Dotplot")) + theme_classic()
  ggsave(file.path(outdir,paste0(cmp,"_dotplot.png")), p_dot, width=10, height=7, dpi=300)
}


library(UpSetR)
upset(fromList(geneLists), order.by="freq")
ggsave("KEGG_upset_UpSetR.png", width=10, height=8, dpi=300)


sig_entrez <- rownames(res_table)[
  res_table$FDR < 0.05 &
    abs(res_table$logFC) > 1
]
sig_entrez <- as.character(sig_entrez)

gene_list <- res_table$logFC
names(gene_list) <- rownames(res_table)
gene_list <- sort(gene_list, decreasing = TRUE)

kk_ora <- enrichKEGG(
  gene         = sig_entrez,
  organism     = "mmu",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10
)

kk_gsea <- gseKEGG(
  geneList     = gene_list,
  organism     = "mmu",
  nPerm        = 10000,
  minGSSize    = 10,
  pvalueCutoff = 0.05
)


dotplot(kk_ora, showCategory=20) +
  ggtitle("Top 20 KEGG Pathways (ORA)") +
  theme_minimal(base_size=14)


firstID <- kk_gsea@result$ID[1]
firstDesc <- kk_gsea@result$Description[1]
gseaplot2(kk_gsea,
          geneSetID = firstID,
          title     = paste("GSEA:", firstDesc),
          base_size = 12)

pathview(
  gene.data  = gene_list,
  pathway.id = firstID,
  species    = "mmu",
  out.suffix = "pathview"
)

p1 <- dotplot(kk_ora, showCategory=15)
ggplotly(p1)


write.csv(as.data.frame(kk_ora), file="KEGG_ORA_results.csv", row.names=FALSE)
write.csv(as.data.frame(kk_gsea@result), file="KEGG_GSEA_results.csv", row.names=FALSE)
ggsave("KEGG_ORA_dotplot.png", width=8, height=6, dpi=300)
ggsave("GSEA_running_score.png", width=8, height=6, dpi=300)

library(GO.db)
cyt.go <- c("GO:0032465", "GO:0000281")
term <- select(GO.db, keys=cyt.go, columns="TERM")
term
message("Available keytypes: ", paste(keytypes(GO.db), collapse = ", "))
message("Available columns: ", paste(columns(GO.db), collapse = ", "))
cyt.go <- c("GO:0032465", "GO:0000281")

go_terms <- AnnotationDbi::select(
  x       = GO.db,
  keys    = cyt.go,
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
)

print(go_terms)


cyt.go <- c("GO:0032465", "GO:0000281")
message("select() belongs to package: ", find("select")[1])
showMethods("select", classes = "OrgDb")
go2entrez <- AnnotationDbi::select(
  x       = org.Mm.eg.db,
  keys    = cyt.go,
  columns = c("ENTREZID"),
  keytype = "GO"
)
print(go2entrez)
write.csv(df_terms,
          file      = "go_terms_results.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")

meta <- y$samples
meta$Group <- factor(meta$group,
                     levels = c("B.virgin","B.pregnant","B.lactating",
                                "L.virgin","L.pregnant","L.lactating"))


design <- model.matrix(~0 + Group, data = meta)
colnames(design) <- levels(meta$Group)
design

cont.matrix <- makeContrasts(
  B.LvsP = B.lactating - B.pregnant,
  B.PvsV = B.pregnant  - B.virgin,
  B.LvsV = B.lactating - B.virgin,
  L.LvsP = L.lactating - L.pregnant,
  L.PvsV = L.pregnant  - L.virgin,
  L.LvsV = L.lactating - L.virgin,
  levels = design
)

library(GO.db)
columns(GO.db)
#> [1] "GOID"     "EVIDENCE" "ONTOLOGY" "TERM"     "DEFINITION"

library(org.Mm.eg.db)
columns(org.Mm.eg.db)        
keytypes(org.Mm.eg.db)       
cyt.go <- c("GO:0032465","GO:0000281")
df <- AnnotationDbi::select(
  x       = org.Mm.eg.db,
  keys    = cyt.go,
  columns = "ENTREZID",
  keytype = "GO"
)
head(df)
cyt.go.genes <- split(df$ENTREZID, f = df$GO)
str(cyt.go.genes)
logCPM <- cpm(y, log=TRUE)
fry_results <- lapply(colnames(contrast.matrix), function(ct) {
  fr <- fry(y      = logCPM,
            index  = cyt.go.genes,
            design = design,
            contrast = contrast.matrix[, ct])
  fr$Contrast <- ct
  fr
})
fry_df <- do.call(rbind, fry_results)

head(fry_df)
fry_df$FDR <- p.adjust(fry_df$PValue, method = "BH")

write.csv(fry_df,
          file      = "GO_enrichment_fry_results.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
library(ggplot2)
fry_df$negLogFDR <- -log10(fry_df$FDR + 1e-10)

ggplot(fry_df, aes(x = Contrast, y = negLogFDR, fill = Contrast)) +
  geom_bar(stat="identity") +
  facet_wrap(~rownames(fry_df), scales="free_y") +
  theme_minimal() +
  labs(y = "-log10(FDR)", x = "")


B.VvsL <- makeContrasts(B.virgin-B.lactating, levels=design)
fry(y, index=cyt.go.genes, design=design, contrast=B.VvsL)

y <- calcNormFactors(y)

meta <- y$samples
meta$Group <- factor(meta$group,
                     levels = c("B.virgin","B.pregnant","B.lactating",
                                "L.virgin","L.pregnant","L.lactating"))
design <- model.matrix(~0 + Group, data = meta)
colnames(design) <- levels(meta$Group)
contrast.matrix <- makeContrasts(
  B.LvsP = B.lactating - B.pregnant,
  B.PvsV = B.pregnant  - B.virgin,
  B.LvsV = B.lactating - B.virgin,
  L.LvsP = L.lactating - L.pregnant,
  L.PvsV = L.pregnant  - L.virgin,
  L.LvsV = L.lactating - L.virgin,
  levels = design
)

logCPM <- cpm(y, log=TRUE)

results_list <- lapply(colnames(contrast.matrix), function(ct) {
  fr <- fry(
    y        = logCPM,
    index    = cyt.go.genes,      
    design   = design,
    contrast = contrast.matrix[, ct]
  )
  fr$Contrast <- ct
  fr
})
names(results_list) <- colnames(contrast.matrix)

all_results <- do.call(rbind, results_list)
print(all_results)
write.csv(all_results, "GO_fry_all_contrasts.csv", row.names = FALSE)


res <- glmQLFTest(fit, contrast=B.VvsL)
index <- rownames(fit) %in% cyt.go.genes[[1]]
barcodeplot(res$table$logFC, index=index, labels=c("B.lactating","B.virgin"), 
            main=cyt.go[1])

y <- calcNormFactors(y)
labels_list <- list(
  B.LvsP = c("B.pregnant",   "B.lactating"),
  B.PvsV = c("B.virgin",     "B.pregnant"),
  B.LvsV = c("B.virgin",     "B.lactating"),
  L.LvsP = c("L.pregnant",   "L.lactating"),
  L.PvsV = c("L.virgin",     "L.pregnant"),
  L.LvsV = c("L.virgin",     "L.lactating")
)
contrast_names <- colnames(contrast.matrix)

for (ct in contrast_names) {
  res <- glmQLFTest(fit, contrast = contrast.matrix[, ct])
  
  pdf(file = paste0("barcodeplot_", ct, ".pdf"),
      width = 8, height = 5, onefile = TRUE)
  
  for (i in seq_along(cyt.go)) {
    go_id   <- cyt.go[i]
    genes_i <- cyt.go.genes[[i]]
    index   <- rownames(res$table) %in% genes_i
    
    barcodeplot(
      res$table$logFC,
      index   = index,
      labels  = labels_list[[ct]],
      main    = paste(ct, go_id, sep = " | ")
    )
  }
  
  dev.off()
  
  message("The barcodeplot for ", ct, " has been saved to file: ",
          paste0("barcodeplot_", ct, ".pdf"))
}
dev.off()

Mm.c2 <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.all.v7.1.entrez.rds"))
Mm.c2
idx <- ids2indices(Mm.c2,id=row.names(y))
BvsL.v <- makeContrasts(B.virgin - L.virgin, levels=design)
cam <- camera(y, idx, design, contrast=BvsL.v, inter.gene.cor=0.01)
options(digits=2)
head(cam,14)
cam

contrast.matrix <- makeContrasts(
  B.LvsP = B.lactating - B.pregnant,
  B.PvsV = B.pregnant  - B.virgin,
  B.LvsV = B.lactating - B.virgin,
  L.LvsP = L.lactating - L.pregnant,
  L.PvsV = L.pregnant  - L.virgin,
  L.LvsV = L.lactating - L.virgin,
  levels = design
)

camera_results <- lapply(colnames(contrast.matrix), function(ct) {
  contrast_vec <- contrast.matrix[, ct]
  
  cam <- camera(y,      
                idx,   
                design,
                contrast = contrast_vec,
                inter.gene.cor = 0.01)
  
  cam$Contrast <- ct
  cam
})

names(camera_results) <- colnames(contrast.matrix)

all_cam <- do.call(rbind, camera_results)
all_cam$Contrast <- rownames(all_cam)
print(head(all_cam))
write.csv(all_cam, "camera_all_contrasts.csv", row.names = FALSE)



pdf("barcodeplots_all_contrasts.pdf", width=8, height=6.4)
for (ct in colnames(contrast.matrix)) {
  
  res <- glmQLFTest(fit, contrast=contrast.matrix[, ct])
  
  logFC_vec <- res$table$logFC
  
  barcodeplot(
    logFC_vec,
    index   = idx[["LIM_MAMMARY_STEM_CELL_UP"]],
    index2  = idx[["LIM_MAMMARY_STEM_CELL_DN"]],
    labels  = c("Up genes","Down genes"),
    main    = paste0(ct, " – Mammary stem cell signature"),
    alpha   = 1
  )
}
Cairo::CairoPDF("barcodeplots.pdf", width=10, height=8)
# ... plotting code ...
dev.off()


dev.off()

res <- glmQLFTest(fit, contrast=BvsL.v)
barcodeplot(res$table$logFC,
            index=idx[["LIM_MAMMARY_STEM_CELL_UP"]],
            index2=idx[["LIM_MAMMARY_STEM_CELL_DN"]],
            labels=c("L.virgin","B.virgin"),
            main="LIM_MAMMARY_STEM_CELL",
            alpha=1)


upSet   <- idx[["LIM_MAMMARY_STEM_CELL_UP"]]
downSet <- idx[["LIM_MAMMARY_STEM_CELL_DN"]]

if (!dir.exists("plots")) dir.create("plots")

for (ct in colnames(contrast.matrix)) {
  
  safe_name <- gsub("[^A-Za-z0-9_]", "_", ct)
  fname     <- file.path("plots", paste0("barcodeplot_", safe_name, ".pdf"))
  
  pdf(file   = fname,
      width  = 8,
      height = 6)
  
  res <- glmQLFTest(fit, contrast = contrast.matrix[, ct])
  
  barcodeplot(
    res$table$logFC,
    index   = upSet,
    index2  = downSet,
    labels  = c("Up genes", "Down genes"),
    main    = paste0("Barcode plot – ", ct),
    alpha   = 1
  )
  
  dev.off()
  
  message("Saved: ", fname)
}

