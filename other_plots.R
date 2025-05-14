# libraries
library(ggpubr)
library(tidyr)

# define path
path_to_project_file <- "~/Dropbox (MIT)/metabolism_Joslin/Jan_Willem/Jan-Willem Human Study_RNAseq"

# LDSC enrichment of GWAS-DEG links
h2 <- read.table(paste0(path_to_project_file, "/human_exercise_obesity_benj/h2-nov27.tsv"), header = T)
h2_plot <- h2[c(1,17,20,22,24,26,28,30,32,34,36,38), c(5,7,8)]
h2_plot$peaks <- gsub("_linked_0_001", "", h2_plot$peaks)
h2_plot[c(7,8,11,12), "peaks"] <- paste0("all_", h2_plot[c(7,8,11,12), "peaks"])
h2_plot <- h2_plot %>% separate(peaks, into = c("sex", "group", "temp1", "temp2", "temp3", "pval"), sep = "_", remove = F)
h2_plot$sex <- factor(h2_plot$sex, levels = c("f", "m"))
h2_plot$group <- factor(h2_plot$group, levels = c("lean", "obese"))
h2_plot <- h2_plot[, -grep("temp", colnames(h2_plot))]
h2_plot$pval <- factor(h2_plot$pval, levels = c("01", "05"))
h2_plot <- h2_plot[order(h2_plot$group, h2_plot$sex),]
h2_plot$peaks <- factor(h2_plot$peaks, levels = h2_plot$peaks)
ggplot(h2_plot, aes(peaks, Enrichment)) +
  geom_point(aes(color = group, shape = sex), size = 3) +
  geom_hline(yintercept = h2_plot[h2_plot$peaks == "BSS01665_Tx", "Enrichment"]) +
  scale_color_manual(values = custom.colours, na.value = "grey50") +
  scale_shape_manual(values = custom.shapes, na.value = 18) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# GWAS TF enrichment
gwas_tf <- read.table(paste0(path_to_project_file, "/human_exercise_obesity_benj/tf_enrich_231128.tsv"), header = T)
