# libraries
library(ggpubr)
library(DESeq2)
library(PCAtools)

# define path
path_to_project_file <- "~/Dropbox (MIT)/metabolism_Joslin/Jan_Willem/Jan-Willem Human Study_RNAseq"

# run data_preprocessing.R
source(paste0(path_to_project_file, "/Data Analysis R Files/data_preprocessing.R"))

# run data_analysis_helper.R
source(paste0(path_to_project_file, "/Data Analysis R Files/data_analysis_helper.R"))

# run plotting.R
source(paste0(path_to_project_file, "/Data Analysis R Files/plotting.R"))

# set working directory
setwd(path_to_project_file)

# PCA_plot function
PCA_plot <- function(vsd, coldata, color_grouping, shape_grouping, label_grouping){
  rv <- matrixStats::rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(2000, length(rv)))] # number of top genes to use for PCs, selected by highest row variance
  pca <- prcomp(t(assay(vsd)[select,]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  pca_anno <- merge(pca$x, coldata, by.x = "row.names", by.y = "row.names")

  plot <- ggplot(pca_anno, aes(PC1, PC2, label = Row.names)) +
    geom_point(aes(color = get(color_grouping), shape = get(shape_grouping)), size = 2) +
    # geom_line(col = "gray") +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    scale_colour_manual(name = color_grouping, values = custom.colours) +
    scale_shape_manual(name = shape_grouping, values = custom.shapes) +
    geom_text_repel(aes(label = get(label_grouping)), min.segment.length = 0) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(list(plot, pca))
}

# set basic variables
bulk_deconvolution <- c("bulk")
wilcoxon_DESeq2 <- c("wilcoxon")
counts_tpm <- c("tpm")

if(counts_tpm == "cts"){
  print("--> COUNTS")
  data <- cts
}else if(counts_tpm == "tpm"){
  print("--> TPM")
  data <- tpm
}else{
  stop("Input 'counts_tpm' should be either 'cts' or 'tpm")
}


# Identify sex-specific genes at baseline using TPM

## Set what and how to analyze
baseline_training <- c("baseline")
comparison <- c("m", "f")

## CLEAN data and metadata
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
baseline_data <- res[[1]]
baseline_coldata <- res[[2]]
## WILCOXON RANK SUM TEST => DEGs
sex_specific_data <- wilcoxon_stats(baseline_data, baseline_coldata, comparison, baseline_training)
sex_specific_data <- drop_rows_with_zeros(sex_specific_data, baseline_coldata) # exclude lowly expressed genes
select_genes <- sex_specific_data[sex_specific_data$p.value > 0.01, "gene.name"] # exclude sex-specific genes with p-values <0.01
sex_specific_gene <- sex_specific_data[sex_specific_data$p.value < 0.01, "gene.name"] # 189 genes


# PCA PLOTS

## Set what and how to analyze
baseline_training <- c("baseline")
comparison <- c("lean", "obese")

## Baseline samples
data <- cts # use counts to generate PCA
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
baseline_data <- res[[1]]
baseline_coldata <- res[[2]]
### Without filtering out sex-specific genes
dds <- DESeqDataSetFromMatrix(countData = baseline_data, colData = baseline_coldata, design = ~large_group)
vsd <- vst(dds, blind = FALSE) # normalize counts
shape_grouping <- "Sex"
color_grouping <- "large_group"
plot <- PCA_plot(vsd, baseline_coldata, color_grouping, shape_grouping)
ggsave(paste0(path_to_project_file, "/Plots/pca/baseline_PCA_plot_all_genes.pdf"), plot[[1]], width=4.5, height=3)

## PCAtools package (use this comprehensive PCA analysis package to double check/further explore existing results)
p <- pca(assay(vsd), metadata = baseline_coldata, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, colby = "large_group", hline = 0, vline = 0, legendPosition = "right", ellipse = TRUE, ellipseType = 't', ellipseLevel = 0.95, ellipseFill = TRUE, ellipseAlpha = 1/4, ellipseLineSize = 1.0,ylim = c(-30, 30))
biplot(p, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p, colby = "sex_large_group")
plotloadings(p, labSize = 3)
eigencorplot(p, components = getComponents(p, 1:29), metavars = c('Age','Sex','group_name','V1_Weight_kg','V1_Body_Fat_percent','V1_Waist_Circ_CM','V1_Fat_Mass_kg','V1_FFM','V1_TBW','V1_BMI','V1_SBP','V1_DBP','V1_RHR_bpm','V3_VO2_Peak_LMIN','V3_Peak_Watts','V1_Albumin','V1_ALK_PHOS','V1_ALT','V1_AST','V1_Cholesterol','V1_Creatinine','V1_HDL^','V1_LDL^','V1_Triglycerides','V1_Total_Bili','V1_Total_Protein','V1_HbA1C_percent','V1_FFA','V1_Insulin','V1_AdipoIR'))

pca <- plot[[2]]
baseline_coldata$group_three <- baseline_coldata$group_name
baseline_coldata[grepl("lean", baseline_coldata$group_three), "group_three"] <- "lean"
color_grouping <- "group_three"
plot <- PCA_plot(vsd, baseline_coldata, color_grouping, shape_grouping)
ggsave(paste0(path_to_project_file, "/Plots/pca/baseline_lean_obese_PCA_plot_all_genes.pdf"), plot[[1]], width=5, height=3)

### Filtering out sex-specific genes (did not make much of a difference)
dds_filt <- DESeqDataSetFromMatrix(countData = baseline_data[rownames(baseline_data) %in% select_genes,], colData = baseline_coldata, design = ~large_group)
vsd_filt <- vst(dds_filt, blind = FALSE) # normalize counts
shape_grouping <- "Sex"
color_grouping <- "large_group"
plot_filt <- PCA_plot(vsd_filt, baseline_coldata, color_grouping, shape_grouping)
ggsave(paste0(path_to_project_file, "/Plots/pca/baseline_PCA_plot_filt.pdf"), plot_filt[[1]], width=4.5, height=3)
pca_filt <- plot_filt[[2]]
color_grouping <- "group_name"
plot_filt <- PCA_plot(vsd_filt, baseline_coldata, color_grouping, shape_grouping)
ggsave(paste0(path_to_project_file, "/Plots/pca/baseline_lean_obese_PCA_plot_filt.pdf"), plot_filt[[1]], width=5, height=3)

## All samples (map the samples after training into the baseline PCA space)
baseline_training <- c("training")
comparison <- c("A", "B")
data <- cts # use counts to generate PCA
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
all_data <- res[[1]]
all_coldata <- res[[2]]
all_data <- all_data[, grepl("B", colnames(all_data))]
all_coldata <- all_coldata[all_coldata$comparison == "B",]
### Without filtering out sex-specific genes
dds_all <- DESeqDataSetFromMatrix(countData = all_data, colData = all_coldata, design = ~1)
vsd_all <- vst(dds_all, blind = FALSE) # normalize counts
B_pca_coord <- scale(t(assay(vsd_all)[names(pca$center),]), pca$center, pca$scale) %*% pca$rotation
pca_coord <- rbind(pca$x, B_pca_coord)
pca_coord <- cbind(pca_coord, coldata[, c("large_group", "group_name", "Sex", "time", "subject")])
plot_all <- ggplot(pca_coord, aes(PC1, PC2)) +
  geom_point(aes(color = large_group, shape = Sex, alpha = time), size = 2) +
  geom_path(aes(x = PC1, y = PC2, group = subject), 
            arrow = arrow(length = unit(0.1, "cm"))) +
  scale_colour_manual(values = custom.colours) +
  scale_shape_manual(values = custom.shapes) +
  scale_alpha_manual(values = c("A" = 1, "B" = 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(paste0(path_to_project_file, "/Plots/pca/all_PCA_plot_all_genes.pdf"), plot_all, width=4.5, height=3)
### Filtering out sex-specific genes
dds_all_filt <- DESeqDataSetFromMatrix(countData = all_data[rownames(all_data) %in% select_genes,], colData = all_coldata, design = ~1)
vsd_all_filt <- vst(dds_all_filt, blind = FALSE) # normalize counts
B_pca_coord_filt <- scale(t(assay(vsd_all_filt)[names(pca_filt$center),]), pca_filt$center, pca_filt$scale) %*% pca_filt$rotation
pca_coord_filt <- rbind(pca_filt$x, B_pca_coord_filt)
pca_coord_filt <- cbind(pca_coord_filt, coldata[, c("large_group", "group_name", "Sex", "time", "subject")])
plot_all_filt <- ggplot(pca_coord_filt, aes(PC1, PC2)) +
  geom_point(aes(color = large_group, shape = Sex, alpha = time), size = 2) +
  geom_path(aes(x = PC1, y = PC2, group = subject), 
            arrow = arrow(length = unit(0.1, "cm"))) +
  scale_colour_manual(values = custom.colours) +
  scale_shape_manual(values = c("f" = 17, "m" = 16)) +
  scale_alpha_manual(values = c("A" = 1, "B" = 0.5)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(paste0(path_to_project_file, "/Plots/pca/all_PCA_plot_filt.pdf"), plot_all_filt, width=4.5, height=3)

# Correlate PCA-space distance and clinical measurement changes
## Calculate euclidean distances in the 29-dimensional space or the 2-dimensional spaces
euclidean <- function(a, b) sqrt(sum((a - b)^2))
for (temp in c("all", "filt")) {
  temp.pca_coord_dist <- data.frame()
  if (temp == "all") {
    temp.pca_coord <- pca_coord
  } else {
    temp.pca_coord <- pca_coord_filt
  }
  for (temp.dim in c(29, 2)) {
    for (temp.sub in unique(temp.pca_coord$subject)) {
      temp.dist <- euclidean(t(temp.pca_coord[temp.pca_coord$subject == temp.sub,][1, 1:temp.dim]), t(temp.pca_coord[temp.pca_coord$subject == temp.sub,][2, 1:temp.dim]))
      temp.pca_coord_dist[temp.sub, paste0("dist_", temp.dim)] <- temp.dist
    }
  }
  temp.pca_coord_dist <- merge(temp.pca_coord_dist, coldata[coldata$time == "B", -4], by.x = "row.names", by.y = "subject")
  assign(paste0("pca_coord_dist_", temp), temp.pca_coord_dist)
  temp.pca_coord_dist_p <- data.frame()
  for (temp.dim in c(29, 2)) { 
    for (temp.col in colnames(coldata)[6:length(colnames(coldata))]) {
      if (temp.col %in% c("large_group", "Sex")) {
        temp.p <- t.test(temp.pca_coord_dist[, paste0("dist_", temp.dim)] ~ temp.pca_coord_dist[, temp.col])$p.value
        temp.pca_coord_dist_p[temp.col, paste0("p_value_", temp.dim)] <- temp.p
      } else if (temp.col %in% c("group_name", "sex_large_group")) {
        for (temp.subgroup in unique(temp.pca_coord_dist$large_group)) {
          temp.data <- temp.pca_coord_dist[temp.pca_coord_dist$large_group == temp.subgroup,]
          temp.p <- t.test(temp.data[, paste0("dist_", temp.dim)] ~ temp.data[, temp.col])$p.value
          temp.pca_coord_dist_p[paste0(temp.col, "_", temp.subgroup), paste0("p_value_", temp.dim)] <- temp.p
        }
      } else {
        temp.p <- cor.test(temp.pca_coord_dist[, paste0("dist_", temp.dim)], temp.pca_coord_dist[, temp.col])$p.value
        temp.pca_coord_dist_p[temp.col, paste0("p_value_", temp.dim)] <- temp.p
      }
    }
    temp.pca_coord_dist_p[, paste0("p_value_adj_", temp.dim)] <- p.adjust(temp.pca_coord_dist_p[, paste0("p_value_", temp.dim)], method = "BH")
    for (temp.var in c(rownames(temp.pca_coord_dist_p[temp.pca_coord_dist_p[, paste0("p_value_", temp.dim)] < 0.05,]), "group_name_obese")) {
      if (temp.var %in% c("group_name_obese", "group_name_lean", "sex_large_group_obese", "sex_large_group_lean")) {
        temp.var.plot <- gsub("_obese", "", gsub("_lean", "", temp.var))
        plot <- ggplot(temp.pca_coord_dist, aes(get(temp.var.plot), get(paste0("dist_", temp.dim)), color = large_group)) +
          geom_boxplot(outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
          geom_point(aes(shape=Sex), position = position_jitterdodge(dodge.width = 1), size = 0.8) +
          scale_colour_manual(values = custom.colours) +
          scale_shape_manual(values = custom.shapes) +
          ylab(paste0("dist_", temp.dim)) +
          xlab(temp.var) +
          stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", paired = FALSE, hide.ns = TRUE) +
          theme_bw() +
          theme(panel.grid = element_blank())
      } else {
        temp.var.plot <- temp.var
        plot <- ggplot(temp.pca_coord_dist, aes(get(paste0("dist_", temp.dim)), get(temp.var.plot), color = large_group)) +
          geom_point(aes(shape = Sex)) +
          scale_colour_manual(values = custom.colours) +
          scale_shape_manual(values = custom.shapes) +
          xlab(paste0("dist_", temp.dim)) +
          ylab(temp.var) +
          geom_smooth(method=lm, se=FALSE) +
          stat_cor(method = "pearson", label.x = 3) +
          theme_bw() +
          theme(panel.grid = element_blank())
      }
      ggsave(paste0(path_to_project_file, paste0("/Plots/pca/all_PCA_dist_", temp.dim, "_", temp.var, "_scatterplot_", temp, ".pdf")), plot, width=4, height=3)
    }
  }
  assign(paste0("pca_coord_dist_p_", temp), temp.pca_coord_dist_p)
}

# Add sex-specific gene column to the DEG excel spreadsheets
temp <- as.data.frame(read_excel(paste0(path_to_project_file, "Data Analysis/training/bulk/wilcoxon/DEGs.xlsx"), sheet = 7, col_names = TRUE))
colnames(temp)[1] <- "gene.name"
temp <- temp[, c("gene.name", "p.value", "mean.B", "mean.A", "FC", "logFC", "min.log.pval", "p.adj")]
temp[temp$gene.name %in% sex_specific_gene, "baseline_sex_specific"] <- "Y"
temp.wb <- loadWorkbook(paste0(path_to_project_file, "supp_tables/TableS5_training_DEGs.xlsx"))
temp.sheets <- getSheets(temp.wb)
# removeSheet(temp.wb, sheetName="Lean")
temp.yourSheet <- createSheet(temp.wb, sheetName="Lean male")
addDataFrame(temp, temp.yourSheet)
saveWorkbook(temp.wb, paste0(path_to_project_file, "supp_tables/TableS5_training_DEGs.xlsx"))
training_obese_sig <- as.data.frame(read_excel(paste0(path_to_project_file, "supp_tables/TableS5_training_DEGs.xlsx"), sheet = 1, col_names = TRUE))
training_obese_sig <- training_obese_sig[abs(training_obese_sig$logFC) > 0.5,]
colnames(training_obese_sig) <- paste0("to_", colnames(training_obese_sig))
training_lean_sig <- as.data.frame(read_excel(paste0(path_to_project_file, "supp_tables/TableS5_training_DEGs.xlsx"), sheet = 2, col_names = TRUE))
training_lean_sig <- training_lean_sig[abs(training_lean_sig$logFC) > 0.5,]
colnames(training_lean_sig) <- paste0("tl_", colnames(training_lean_sig))
baseline_obese_vs_lean_sig <- as.data.frame(read_excel(paste0(path_to_project_file, "supp_tables/TableS2_baseline_DEGs.xlsx"), sheet = 1, col_names = TRUE))
baseline_obese_vs_lean_sig <- baseline_obese_vs_lean_sig[abs(baseline_obese_vs_lean_sig$logFC) > 0.5,]
colnames(baseline_obese_vs_lean_sig) <- paste0("ol_", colnames(baseline_obese_vs_lean_sig))
global_ol <- merge(training_obese_sig, baseline_obese_vs_lean_sig, by.x = "to_gene.name", by.y = "ol_gene.name")
global_ol$dir <- ifelse(sign(global_ol$to_logFC) == -sign(global_ol$ol_logFC), "oppo", "same")
global_ol_oppo <- global_ol[global_ol$dir == "oppo", "to_gene.name"]

# Correlation between gene-level change and clinical measurement changes
baseline_training <- c("training")
comparison <- c("A", "B")
data <- cts
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
all_data <- res[[1]]
all_coldata <- res[[2]]
dds <- DESeqDataSetFromMatrix(countData = all_data, colData = all_coldata, design = ~1)
vsd <- vst(dds, blind = FALSE) # normalize counts
candidate_gene_list <- training_obese_sig$gene.name[training_obese_sig$gene.name %in% c(training_lean_sig$gene.name, baseline_obese_vs_lean_sig$gene.name)]
candidate_gene_mat <- assay(vsd)[candidate_gene_list,]
candidate_gene_delta <- candidate_gene_mat[, 30:58] - candidate_gene_mat[, 1:29]
colnames(candidate_gene_delta) <- c(1,2,seq(4,30))
pca_coord_dist_all <- merge(pca_coord_dist_all, t(candidate_gene_delta), by.x = "Row.names", by.y = "row.names")
for (temp.gene in candidate_gene_list) {
  for (temp.col in colnames(pca_coord_dist_all)[c(2,3,100:108)]) {
    if (temp.col %in% c("large_group", "Sex")) {
      temp.p <- t.test(pca_coord_dist_all[, temp.gene] ~ pca_coord_dist_all[, temp.col])$p.value
      pca_coord_dist_p_all[temp.col, paste0("p_value_", temp.gene)] <- temp.p
    } else if (temp.col %in% c("group_name", "sex_large_group")) {
      for (temp.subgroup in unique(pca_coord_dist_all$large_group)) {
        temp.data <- pca_coord_dist_all[pca_coord_dist_all$large_group == temp.subgroup,]
        temp.p <- t.test(temp.data[, temp.gene] ~ temp.data[, temp.col])$p.value
        pca_coord_dist_p_all[paste0(temp.col, "_", temp.subgroup), paste0("p_value_", temp.gene)] <- temp.p
      }
    } else {
      temp.p <- cor.test(pca_coord_dist_all[, temp.gene], pca_coord_dist_all[, temp.col])$p.value
      pca_coord_dist_p_all[temp.col, paste0("p_value_", temp.gene)] <- temp.p
    }
  }
  pca_coord_dist_p_all[, paste0("p_value_adj_", temp.gene)] <- p.adjust(pca_coord_dist_p_all[, paste0("p_value_", temp.gene)], method = "BH")
  if (sum(pca_coord_dist_p_all[, paste0("p_value_", temp.gene)] < 0.005) > 0) {
    for (temp.var in rownames(pca_coord_dist_p_all[pca_coord_dist_p_all[, paste0("p_value_", temp.gene)] < 0.005,])) {
      if (temp.var %in% c("large_group", "sex_large_group_lean", "sex_large_group_obese", "group_name_obese", "group_name_lean")) {
        temp.var.plot <- gsub("_obese", "", gsub("_lean", "", temp.var))
        plot <- ggplot(pca_coord_dist_all, aes(get(temp.var.plot), get(temp.gene), color = large_group)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(shape = Sex), width = 0.2) +
          scale_colour_manual(values = custom.colours) +
          scale_shape_manual(values = c("f" = 17, "m" = 16)) +
          ylab(paste0("Delta_", temp.gene)) +
          xlab(temp.var) +
          theme_bw() +
          theme(panel.grid = element_blank())
      } else if (temp.var == "Sex") {
        plot <- ggplot(pca_coord_dist_all, aes(get(temp.var), get(temp.gene))) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(aes(shape = Sex, color = large_group), width = 0.2) +
          scale_colour_manual(values = custom.colours) +
          scale_shape_manual(values = c("f" = 17, "m" = 16)) +
          ylab(paste0("Delta_", temp.gene)) +
          xlab(temp.var) +
          theme_bw() +
          theme(panel.grid = element_blank())
        } else {
        plot <- ggplot(pca_coord_dist_all, aes(get(temp.gene), get(temp.var), color = large_group)) +
          geom_point(aes(shape = Sex)) +
          scale_colour_manual(values = custom.colours) +
          scale_shape_manual(values = c("f" = 17, "m" = 16)) +
          xlab(paste0("Delta_", temp.gene)) +
          ylab(temp.var) +
          geom_smooth(method=lm, se=FALSE) +
          stat_cor(method = "pearson", label.x = 0) +
          theme_bw() +
          theme(panel.grid = element_blank())
      }
      ggsave(paste0(path_to_project_file, paste0("/Plots/pca/Delta_", temp.gene, "_", temp.var, "_scatterplot.pdf")), plot, width=4, height=3)
    }
  }
}
b3gnt10 <- as.data.frame(t(candidate_gene_mat["B3GNT10",, drop = F]))
b3gnt10$timepoint <- ifelse(grepl("A", rownames(b3gnt10)), "A", "B")
b3gnt10$subject <- gsub("A", "", gsub("B", "", rownames(b3gnt10)))
rownames(coldata) <- paste0(coldata$time, coldata$subject)
b3gnt10 <- merge(b3gnt10, coldata[, c("Sex", "large_group", "delta_Creatinine", "Days_between_last_ex_session")], by.x = "row.names", by.y = "row.names")
plot <- ggplot(b3gnt10, aes(timepoint, B3GNT10, color = large_group, group = subject)) +
  geom_point(aes(shape = Sex)) +
  geom_line() +
  scale_colour_manual(values = custom.colours) +
  scale_shape_manual(values = c("f" = 17, "m" = 16)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(paste0(path_to_project_file, paste0("/Plots/pca/B3GNT10_lineplot.pdf")), plot, width=4, height=3)
