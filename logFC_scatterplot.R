# libraries
library(ggpubr)

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



## TRAINING EFFECTS

# Set what and how to analyze
baseline_training <- c("training")
comparison <- c("B", "A")

# CLEAN data and metadata
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
training_data <- res[[1]]
training_coldata <- res[[2]]
for (temp.selection in c("obese", "lean", "obese female", "obese male", "lean female", "lean male", "non-diabetic obese", "T2D obese", "MIT lean", "HIIT lean")) {
  # SELECT specific individuals for the comparison
  temp.res <- select_individuals(temp.selection, training_data, training_coldata)
  temp.data <- temp.res[[1]]
  temp.coldata <- temp.res[[2]]
  # WILCOXON RANK SUM TEST => DEGs
  temp.data <- wilcoxon_stats(temp.data, temp.coldata, comparison, baseline_training)
  # significant genes with many expression levels equal to zero as are no longer significant
  temp.data <- drop_rows_with_zeros(temp.data, temp.coldata)
  temp.data <- temp.data[temp.data$p.value != 1,]
  temp.sig <- temp.data[temp.data$p.value < 0.05 & abs(temp.data$logFC) > 0.5,]
  assign(paste0(baseline_training, "_", gsub(" ", "_", temp.selection), "_data"), temp.data)
  assign(paste0(baseline_training, "_", gsub(" ", "_", temp.selection), "_sig"), temp.sig)
}


## BASELINE OBESE VS LEAN

# Set what and how to analyze
baseline_training <- c("baseline")
comparison <- c("obese", "lean")
selection <- c("all")

# CLEAN data and metadata
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
baseline_obese_vs_lean_data <- res[[1]]
baseline_obese_vs_lean_coldata <- res[[2]]
# SELECT specific individuals for the comparison
res <- select_individuals(selection, baseline_obese_vs_lean_data, baseline_obese_vs_lean_coldata)
baseline_obese_vs_lean_data <- res[[1]]
baseline_obese_vs_lean_coldata <- res[[2]]
# WILCOXON RANK SUM TEST => DEGs
baseline_obese_vs_lean_data <- wilcoxon_stats(baseline_obese_vs_lean_data, baseline_obese_vs_lean_coldata, comparison, baseline_training)
# significant genes with many expression levels equal to zero as are no longer significant
baseline_obese_vs_lean_data <- drop_rows_with_zeros(baseline_obese_vs_lean_data, baseline_obese_vs_lean_coldata)
baseline_obese_vs_lean_data <- baseline_obese_vs_lean_data[baseline_obese_vs_lean_data$p.value != 1,]
baseline_obese_vs_lean_sig <- baseline_obese_vs_lean_data[baseline_obese_vs_lean_data$p.value < 0.05 & abs(baseline_obese_vs_lean_data$logFC) > 0.5,]

## BASELINE OBESE T2D VS OBESE NORMAL GLYCEMIC CONTROL

# Set what and how to analyze
baseline_training <- c("baseline")
comparison <- c("T2D obese", "non-diabetic obese")
selection <- c("all")

# CLEAN data and metadata
res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
baseline_obese_t2d_vs_norm_data <- res[[1]]
baseline_obese_t2d_vs_norm_coldata <- res[[2]]
# SELECT specific individuals for the comparison
res <- select_individuals(selection, baseline_obese_t2d_vs_norm_data, baseline_obese_t2d_vs_norm_coldata)
baseline_obese_t2d_vs_norm_data <- res[[1]]
baseline_obese_t2d_vs_norm_coldata <- res[[2]]
# WILCOXON RANK SUM TEST => DEGs
baseline_obese_t2d_vs_norm_data <- wilcoxon_stats(baseline_obese_t2d_vs_norm_data, baseline_obese_t2d_vs_norm_coldata, comparison, baseline_training)
# significant genes with many expression levels equal to zero as are no longer significant
baseline_obese_t2d_vs_norm_data <- drop_rows_with_zeros(baseline_obese_t2d_vs_norm_data, baseline_obese_t2d_vs_norm_coldata)
baseline_obese_t2d_vs_norm_data <- baseline_obese_t2d_vs_norm_data[baseline_obese_t2d_vs_norm_data$p.value != 1,]
baseline_obese_t2d_vs_normn_sig <- baseline_obese_t2d_vs_norm_data[baseline_obese_t2d_vs_norm_data$p.value < 0.05 & abs(bbaseline_obese_t2d_vs_norm_data$logFC) > 0.5,]


## SCATTERPLOTS

# training effect obese AND training effect lean
scatterplot_axis <- c("training_obese", "training_lean")
sp_data <- list(training_obese_data, training_lean_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both", cap_threshold = 2)
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/obese-lean comparison/logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect obese AND baseline obese vs lean (baseline as x-axis, training as y-axis)
scatterplot_axis <- c("baseline_obese_vs_lean", "training_obese")
sp_data <- list(baseline_obese_vs_lean_data, training_obese_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both", cap_threshold = 2)
ggsave(paste0(path_to_project_file, "/Plots/baseline-training comparison/bulk/wilcoxon/obesity_rescue_logFC_scatterplot_with_R_P.pdf"), plot, width=4, height=4)

# training effect lean AND baseline obese vs lean
scatterplot_axis <- c("baseline_obese_vs_lean", "training_lean")
sp_data <- list(baseline_obese_vs_lean_data, training_lean_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/baseline-training comparison/bulk/wilcoxon/logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect lean MIT AND training effect lean HIIT
scatterplot_axis <- c("training_lean_MIT", "training_lean_HIIT")
sp_data <- list(training_lean_MIT_data, training_lean_HIIT_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/baseline-training comparison/bulk/wilcoxon/lean_MIT_HIIT_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect lean MIT AND baseline obese vs lean
scatterplot_axis <- c("training_lean_MIT", "baseline_obese_vs_lean")
sp_data <- list(training_lean_MIT_data, baseline_obese_vs_lean_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/baseline-training comparison/bulk/wilcoxon/lean_MIT_baseline_obese_lean_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect lean HIIT AND baseline obese vs lean
scatterplot_axis <- c("training_lean_HIIT", "baseline_obese_vs_lean")
sp_data <- list(training_lean_HIIT_data, baseline_obese_vs_lean_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/baseline-training comparison/bulk/wilcoxon/lean_HIIT_baseline_obese_lean_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training female effect obese AND training male effect obese
scatterplot_axis <- c("training_obese_female", "training_obese_male")
sp_data <- list(training_obese_female_data, training_obese_male_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis)
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_female_male_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect obese normal glycemic control AND training effect obese T2D
scatterplot_axis <- c("training_obese_norm", "training_obese_t2d")
sp_data <- list(training_obese_norm_data, training_obese_t2d_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_norm_t2d_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect obese T2D AND baseline obese T2D vs obese
scatterplot_axis <- c("baseline_obese_t2d_vs_norm", "training_obese_t2d")
sp_data <- list(baseline_obese_t2d_vs_norm_data, training_T2D_obese_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both", cap_threshold = 2)
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_t2d_baseline_t2d_norm_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect obese T2D AND baseline obese T2D vs lean
scatterplot_axis <- c("training_obese_t2d", "baseline_obese_t2d_lean")
sp_data <- list(training_obese_t2d_data, baseline_obese_t2d_lean_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_t2d_baseline_t2d_lean_logFC_scatterplot.pdf"), plot, width=4, height=4)

# training effect obese normal glycemic control AND baseline obese normal glycemic control vs lean
scatterplot_axis <- c("training_obese_norm", "baseline_obese_norm_lean")
sp_data <- list(training_obese_norm_data, baseline_obese_norm_lean_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, regression_line = "all", label = "both")
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_norm_baseline_norm_lean_logFC_scatterplot.pdf"), plot, width=4, height=4)

for (temp.tf in c("ARNTL", "AR", "PPARG", "HNF1A")) {
  plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, highlight = regulons[regulons$tf == temp.tf, "target"], regression_line = "highlight", label = "highlight")
  ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_female_male_logFC_scatterplot_", temp.tf, ".pdf"), plot, width=4, height=4)
}
for (i in 1:length(select_pathways)) {
  plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, highlight = get(paste0("pathway", i)), regression_line = "highlight", label = "highlight")
  ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/obese_female_male_logFC_scatterplot_", select_pathways[i], ".pdf"), plot, width=4, height=4)
}

# training female effect lean AND training male effect lean
scatterplot_axis <- c("training_lean_female", "training_lean_male")
sp_data <- list(training_lean_female_data, training_lean_male_data)
plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis)
ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/lean_female_male_logFC_scatterplot.pdf"), plot, width=4, height=4)
for (temp.tf in c("ARNTL", "AR", "PPARG", "HNF1A")) {
  plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, highlight = regulons[regulons$tf == temp.tf, "target"], regression_line = "highlight", label = "highlight")
  ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/lean_female_male_logFC_scatterplot_", temp.tf, ".pdf"), plot, width=4, height=4)
}
for (i in 1:length(select_pathways)) {
  plot <- compare_DEGs_scatterplot(sp_data, scatterplot_axis, highlight = get(paste0("pathway", i)), regression_line = "highlight", label = "highlight")
  ggsave(paste0(path_to_project_file, "/Plots/training/bulk/wilcoxon/comparison/lean_female_male_logFC_scatterplot_", select_pathways[i], ".pdf"), plot, width=4, height=4)
}

# interesting TF regulons (AR 138 targets, ARNTL 10 targets, TCF7L2 20 targets, HNF1A 27 targets, PPARG 41 targets)
data(dorothea_hs, package = "dorothea")
regulons <- as.data.frame(dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B")))
i <- 1
for (temp.sheet in c("obese", "lean")) {
  temp.data <- as.data.frame(read_excel(paste0(path_to_project_file, "supp_tables/TableS7_training_TF_activity_Dorothea.xlsx"), sheet = i, col_names = TRUE))
  temp.data <- temp.data[order(temp.data$NES, decreasing = T),]
  assign(paste0("tf_", temp.sheet), temp.data)
  i <- i+1
}
sum(tf_obese[1:20, "TF"] %in% tf_lean[1:20, "TF"]) # 9
sum(tf_obese[(nrow(tf_obese)-19):nrow(tf_obese), "TF"] %in% tf_lean[(nrow(tf_lean)-19):nrow(tf_lean), "TF"]) # 4
tf_list <- c("ARNTL", "PPARG", "AR")
for (i in 1:length(tf_list)) {
  temp.genes <- regulons[regulons$tf == tf_list[i], "target"]
  assign(paste0("tf", i), temp.genes)
}

# interesting biological pathways
i <- 1
for (temp.sheet in c("obese", "lean")) {
  temp.data <- as.data.frame(read_excel(paste0(path_to_project_file, "supp_tables/TableS6_training_Pathways.xlsx"), sheet = i, col_names = TRUE))
  assign(paste0("pathway_", temp.sheet), temp.data)
  i <- i+1
}
sum(pathway_obese$pathway %in% pathway_lean$pathway)
select_pathways <- c("NABA_ECM_REGULATORS", "WP_PPAR_SIGNALING_PATHWAY", "REACTOME_MOLECULES_ASSOCIATED_WITH_ELASTIC_FIBRES")
for (i in 1:length(select_pathways)) {
  temp.genes <- unlist(strsplit(pathway_obese[pathway_obese$pathway == select_pathways[i], "leadingEdge"], ", "))
  temp.genes <- c(temp.genes, unlist(strsplit(pathway_lean[pathway_lean$pathway == select_pathways[i], "leadingEdge"], ", ")))
  temp.genes <- unique(temp.genes)
  assign(paste0("pathway", i), temp.genes)
}


compare_DEGs_scatterplot <- function(scatterplot_data, axis, pval_threshold=0.05, logFC_threshold=0.5, cap_threshold=3, highlight=NULL, regression_line=c("all","highlight"), label=c("both", "highlight")){
  data.x <- scatterplot_data[[1]]
  data.y <- scatterplot_data[[2]]
  
  data.x <- data.x[, names(data.x) %in% c("p.value", "logFC")]
  data.y <- data.y[, names(data.y) %in% c("p.value", "logFC")]
  
  scatterplot_data <- merge(data.x, data.y, by="row.names")
  names(scatterplot_data)[names(scatterplot_data) == "Row.names"] <- "gene"

  scatterplot_data$significant <- "none"
  scatterplot_data$sig.x <- ifelse(scatterplot_data$p.value.x<pval_threshold & abs(scatterplot_data$logFC.x) > logFC_threshold, TRUE, FALSE)
  scatterplot_data$sig.y <- ifelse(scatterplot_data$p.value.y<pval_threshold & abs(scatterplot_data$logFC.y) > logFC_threshold, TRUE, FALSE)
  if (sum(scatterplot_data$sig.x & scatterplot_data$sig.y) > 0) {
    scatterplot_data[scatterplot_data$sig.x & scatterplot_data$sig.y,]$significant <- "both"
  }
  if (sum(scatterplot_data$sig.x & !scatterplot_data$sig.y) > 0) {
    scatterplot_data[scatterplot_data$sig.x & !scatterplot_data$sig.y,]$significant <- axis[1]
  }
  if (sum(!scatterplot_data$sig.x & scatterplot_data$sig.y) > 0) {
    scatterplot_data[!scatterplot_data$sig.x & scatterplot_data$sig.y,]$significant <- axis[2]
  }
  
  scatterplot_data$capped <- 0
  # cap if abs(logFC)>cap_threshold
  if(length(scatterplot_data[abs(scatterplot_data$logFC.x)>cap_threshold | abs(scatterplot_data$logFC.y)>cap_threshold,]$capped) > 0){
    scatterplot_data[abs(scatterplot_data$logFC.x)>cap_threshold | abs(scatterplot_data$logFC.y)>cap_threshold,]$capped <- 1
  }
  scatterplot_data[abs(scatterplot_data$logFC.x)>cap_threshold, "logFC.x"] <- sign(scatterplot_data[abs(scatterplot_data$logFC.x)>cap_threshold, "logFC.x"])*cap_threshold
  scatterplot_data[abs(scatterplot_data$logFC.y)>cap_threshold, "logFC.y"] <- sign(scatterplot_data[abs(scatterplot_data$logFC.y)>cap_threshold, "logFC.y"])*cap_threshold
  
  # PLOT
  plot <- ggplot(scatterplot_data, aes(x=logFC.x, y=logFC.y)) +
    scale_x_continuous(limits = c(-cap_threshold,cap_threshold), expand = expansion(add = 0.05)) +
    scale_y_continuous(limits = c(-cap_threshold,cap_threshold), expand = expansion(add = 0.05)) +
    # none and capped
    geom_point(data = scatterplot_data[scatterplot_data$significant == "none" & scatterplot_data$capped==1, ],
               aes(x=logFC.x, y=logFC.y), shape=25, fill="grey40", color="grey45", alpha = 0.5, size=0.7) +
    geom_point(data = scatterplot_data[scatterplot_data$significant == "none" & scatterplot_data$capped==0, ],
               aes(x=logFC.x, y=logFC.y), shape=21, fill="grey40", color="grey45", alpha = 0.5, size=0.5)
  if(grepl("female", axis[1])) {
    temp.shape.x <- 17
  } else if (grepl("male", axis[1])) {
    temp.shape.x <- 16
  } else {
    temp.shape.x <- 21
  }
  if(grepl("female", axis[2])) {
    temp.shape.y <- 17
  } else if (grepl("male", axis[2])) {
    temp.shape.y <- 16
  } else {
    temp.shape.y <- 21
  }
  for (temp.axis in c(1, 2)) {
    if(grepl("training_obese", axis[temp.axis])) {
      assign(paste0("temp.fill", temp.axis), custom.colours[["obese"]])
      assign(paste0("temp.color", temp.axis), "orangered2")
    } else if (grepl("training_lean", axis[temp.axis])) {
      assign(paste0("temp.fill", temp.axis), custom.colours[["lean"]])
      assign(paste0("temp.color", temp.axis), "green3")
    } else {
      assign(paste0("temp.fill", temp.axis), "#00B6EB")
      assign(paste0("temp.color", temp.axis), "royalblue1")
    } 
  }
  
  # axis[1]
  temp.axis1.label <- paste0("FC for ", axis[1])
  plot <- plot +
    geom_point(data = scatterplot_data[scatterplot_data$significant == axis[1] & scatterplot_data$capped==1, ],
               aes(x=logFC.x, y=logFC.y), shape=25, fill=temp.fill1, color=temp.color1, alpha = 0.9, size=1) +
    geom_point(data = scatterplot_data[scatterplot_data$significant == axis[1] & scatterplot_data$capped==0, ],
               aes(x=logFC.x, y=logFC.y), shape=temp.shape.x, fill=temp.fill1, color=temp.color1, alpha = 0.9, size=1) +
    labs(x = bquote(~log[2]~ .(temp.axis1.label)))

  # axis[2]
  temp.axis2.label <- paste0("FC for ", axis[2])
  plot <- plot +
    geom_point(data = scatterplot_data[scatterplot_data$significant == axis[2] & scatterplot_data$capped==1, ],
               aes(x=logFC.x, y=logFC.y), shape=25, fill=temp.fill2, color=temp.color2, alpha = 0.9, size=1) +
    geom_point(data = scatterplot_data[scatterplot_data$significant == axis[2] & scatterplot_data$capped==0, ],
               aes(x=logFC.x, y=logFC.y), shape=temp.shape.y, fill=temp.fill2, color=temp.color2, alpha = 0.9, size=1) +
    labs(y = bquote(~log[2]~ .(temp.axis2.label)))

  # both and capped 
  plot <- plot +
    geom_point(data = scatterplot_data[scatterplot_data$significant == "both" & scatterplot_data$capped=="1", ],
             aes(x=logFC.x, y=logFC.y), shape=25, fill="grey14", color="black", alpha = 0.9, size=1.7) +
    geom_point(data = scatterplot_data[scatterplot_data$significant == "both" & scatterplot_data$capped=="0", ],
               aes(x=logFC.x, y=logFC.y), shape=21, fill="grey14", color="black", alpha = 0.9, size=1.5)
  
  # candidate genes to highlight
  if (!is.null(highlight)) {
    plot <- plot +
      geom_point(data = scatterplot_data[scatterplot_data$gene %in% highlight & scatterplot_data$capped=="1", ],
                 aes(x=logFC.x, y=logFC.y), shape=25, fill="deeppink", color="deeppink4", alpha = 0.9, size=1.7) +
      geom_point(data = scatterplot_data[scatterplot_data$gene %in% highlight & scatterplot_data$capped=="0", ],
                 aes(x=logFC.x, y=logFC.y), shape=21, fill="deeppink", color="deeppink4", alpha = 0.9, size=1.5)
  }
    
  plot <- plot +
    geom_hline(yintercept = 0.5, size = 0.3, linetype = "dashed") +
    geom_hline(yintercept = -0.5, size = 0.3, linetype = "dashed") +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_vline(xintercept = 0.5, size = 0.3, linetype = "dashed") +
    geom_vline(xintercept = -0.5, size = 0.3, linetype = "dashed") +
    geom_vline(xintercept = 0, size = 0.3) +
    theme_classic()
  
  if (regression_line == "all") {
    plot <- plot +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x.npc = 0.7, label.y.npc = "top", size = 3) +
      geom_smooth(method=glm, se=FALSE, fullrange=TRUE, size=0.5, color = "black")
  } else if (regression_line == "highlight") {
    plot <- plot +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), data=scatterplot_data[scatterplot_data$gene %in% highlight,], label.x.npc = 0.7, label.y.npc = "top", size = 3) +
      geom_smooth(data=scatterplot_data[scatterplot_data$gene %in% highlight,], method=glm, se=FALSE, fullrange=TRUE, size=0.5, color = "deeppink4")
  }
    
  if (label == "both") {
    plot <- plot +
      geom_text_repel(data = scatterplot_data[scatterplot_data$significant == "both" & abs(scatterplot_data$logFC.x)>logFC_threshold & abs(scatterplot_data$logFC.y)>logFC_threshold, ],
                      aes(x=logFC.x, y=logFC.y, label=gene), size=3, segment.size = 0.25, color="black", fontface = 'bold', max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
  } else if (label == "highlight") {
    plot <- plot +
      geom_text_repel(data = scatterplot_data[scatterplot_data$gene %in% highlight, ],
                      aes(x=logFC.x, y=logFC.y, label=gene), size=3, segment.size = 0.25, color="deeppink4", fontface = 'bold', max.overlaps = getOption("ggrepel.max.overlaps", default = 30))
  }

  return(plot)
}
