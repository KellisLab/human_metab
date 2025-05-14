library(reshape2)
library(ggpubr)

path_to_project_file <- "~/Dropbox (MIT)/metabolism_Joslin/Jan_Willem/Jan-Willem Human Study_RNAseq/"
path_to_intial_data <- paste0(path_to_project_file, "/Data")

# run data_preprocessing.R
source(paste0(path_to_project_file, "/Data Analysis R Files/data_preprocessing.R"))

prop <- read.csv(paste0(path_to_intial_data, "/CIBERSORTx_Job19_output/CIBERSORTxGEP_Job19_Fractions-Adjusted.txt"), sep="\t")

prop <- prop[, !names(prop) %in% c("P.value", "Correlation", "RMSE")]
prop <- merge(prop, coldata[, names(coldata) %in% c("subject", "time", "Sex", "large_group", "sex_large_group", "group_name", "V1_BMI")], by.x="Mixture", by.y="row.names")
prop <- reshape2::melt(prop, id.vars = c("Mixture", "subject", "time", "Sex", "large_group", "sex_large_group", "group_name", "V1_BMI"))
names(prop)[names(prop) == "variable"] <- "cell_type"
levels(prop$cell_type)[levels(prop$cell_type)=="ASC"] <- "Adipose Stem Cell"
levels(prop$cell_type)[levels(prop$cell_type)=="EC"] <- "Endothelial"
levels(prop$cell_type)[levels(prop$cell_type)=="SM"] <- "Smooth muscle"
levels(prop$cell_type)[levels(prop$cell_type)=="Lymphoid.cell"] <- "Lymphoid"
levels(prop$cell_type)[levels(prop$cell_type)=="Myeloid.cell"] <- "Myeloid"
prop$cell_type <- factor(prop$cell_type, levels = names(custom.colours))

subject_time <- c(1:30)
subject_time <- subject_time[subject_time!= 3]
a <- paste0("A", subject_time)
b <- paste0("B", subject_time)
#subject_time_factor <- c(sapply(seq_along(a), function(i) append(a[i], b[i], i)))
subject_time_factor <- c(a,b)

prop$Mixture <- factor(prop$Mixture,levels = subject_time_factor)

prop[prop$large_group == "lean", "large_group"] <- "lean"
prop[prop$large_group == "obese", "large_group"] <- "obese"
prop$large_group <- factor(prop$large_group, levels = c("lean", "obese"))
prop$overweight <- ifelse(prop$V1_BMI >= 25 & prop$V1_BMI < 30, "yes", NA)

# Baseline - lean vs obese (large group)
p1 <- ggplot(prop[prop$time == "A",], aes(x=cell_type, y=value, group = large_group)) + 
  geom_boxplot(aes(fill=large_group), outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
  geom_point(aes(shape=Sex, color=overweight), position = position_dodge(1)) +
  labs(x = "Cell type", y = "Relative proportion", fill = "Sex") +
  geom_text_repel(aes(label = V1_BMI), data = prop[prop$time == "A" & prop$overweight == "yes" & !is.na(prop$overweight),], min.segment.length = 0, position = position_dodge(1)) +
  scale_fill_manual(values = custom.colours) +
  scale_color_manual(values = c("red", "black")) +
  scale_shape_manual(values = custom.shapes) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = FALSE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/baseline/obese vs lean/deconvolution/proportion/proportion_overweight_label.pdf"), p1, width=10, height=4)

# Baseline - lean vs overweight vs obese (three group)
prop$three_group <- ifelse(!is.na(prop$overweight), "overweight", NA)
prop[is.na(prop$three_group), "three_group"] <- as.character(prop[is.na(prop$three_group), "large_group"])
prop$three_group <- factor(prop$three_group, levels = c("lean", "overweight", "obese"))
p1 <- ggplot(prop[prop$time == "A",], aes(x=cell_type, y=value, group = three_group)) + 
  geom_boxplot(aes(fill=three_group), outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
  geom_point(aes(shape=Sex), position = position_dodge(1)) +
  facet_grid(~cell_type, scales = "free_x") +
  labs(x = "Cell type", y = "Relative proportion", fill = "Sex") +
  scale_fill_manual(values = custom.colours) +
  scale_shape_manual(values = custom.shapes) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  geom_pwc(
    aes(group = three_group), tip.length = 0,
    method = "wilcox_test", label = "p.format"
  )

ggsave(paste0(path_to_project_file, "/Plots/baseline/obese vs lean/deconvolution/proportion/proportion_three_group.pdf"), p1, width=10, height=4)

p2 <- ggplot(prop[prop$time == "A",], aes(x=cell_type, y=value, group = large_group)) + 
  geom_boxplot(aes(fill = large_group), outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
  geom_point(aes(shape=Sex), position = position_jitterdodge(dodge.width = 1), size = 0.8) +
  facet_grid(~cell_type, scales = "free_x") +
  labs(x = "Cell type", y = "Relative proportion", fill = "Sex") +
  scale_fill_manual(values = custom.colours) +
  scale_shape_manual(values = custom.shapes) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = FALSE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/baseline/obese vs lean/deconvolution/proportion/proportion2.pdf"), p2, width=6, height=3)


# Baseline - male vs female
p1 <- ggplot(prop[prop$time == "A",], aes(x=cell_type, y=value, fill=Sex)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge2(padding = 0.2)) +
  geom_point(aes(shape=Sex), position = position_jitterdodge(jitter.width = 0.2), size = 0.8) +
  facet_grid(large_group~cell_type, scales = "free_x") +
  labs(x = "Cell type", y = "Relative proportion", fill = "Sex") +
  scale_fill_manual(values = custom.colours) +
  scale_shape_manual(values = custom.shapes) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = FALSE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/baseline/male vs female/deconvolution/proportion/proportion.pdf"), p1, width=6, height=3)


# Obese T2D vs obese
p1 <- ggplot(prop[prop$time == "A" & prop$large_group == "obese",], aes(x=cell_type, y=value, group=group_name)) + 
  geom_boxplot(aes(fill=group_name), outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
  geom_point(aes(shape=Sex, color=overweight), position = position_dodge(1)) +
  facet_grid(.~cell_type, scales = "free_x") +
  labs(x = "Cell type", y = "Relative proportion", fill = "Sex") +
  scale_fill_manual(values = custom.colours) +
  scale_shape_manual(values = custom.shapes) +
  scale_color_manual(values = c("red", "black")) +
  geom_text_repel(aes(label = V1_BMI), data = prop[prop$time == "A" & prop$large_group == "obese" & prop$overweight == "yes" & !is.na(prop$overweight),], min.segment.length = 0, position = position_dodge(1)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = FALSE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/baseline/T2D obese vs non-diabetic obese/deconvolution/proportion/proportion_overweight_label.pdf"), p1, width=8, height=3)

p2 <- ggplot(prop[prop$time == "A" & prop$large_group == "obese" & prop$V1_BMI > 27,], aes(x=cell_type, y=value, group=group_name)) + 
  geom_boxplot(aes(fill=group_name), outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
  geom_point(aes(shape=Sex), position = position_dodge(1)) +
  facet_grid(.~cell_type, scales = "free_x") +
  labs(x = "Cell type", y = "Relative proportion", fill = "Sex") +
  scale_fill_manual(values = custom.colours) +
  scale_shape_manual(values = custom.shapes) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = FALSE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/baseline/T2D obese vs non-diabetic obese/deconvolution/proportion/proportion_two_BMI_under_27_removed.pdf"), p2, width=8, height=3)

# post-training vs pre-training
p <- ggplot(prop, aes(x=time, y=value, fill=time)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge2(padding = 0.2)) +
  geom_point(aes(shape=Sex), size = 0.9, color = "red") +
  geom_line(aes(group=subject), size = 0.3) +
  facet_grid(large_group~cell_type) +
  labs(x = "Cell type", y = "Relative proportion", fill = "time") +
  scale_fill_manual(values = custom.colours) +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = TRUE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/training/deconvolution/proportion/proportion_with_lines.pdf"), p, width=7, height=3)

p1 <- ggplot(prop, aes(x=time, y=value, fill=time)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge(1)) +
  geom_point(aes(shape=Sex), size = 0.9, color = "red") +
  geom_line(aes(group=subject), size = 0.3) +
  facet_grid(large_group+Sex~cell_type) +
  labs(x = "Cell type", y = "Relative proportion", fill = "time") +
  scale_fill_manual(values = custom.colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = TRUE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/training/deconvolution/proportion/proportion_by_sex_with_lines.pdf"), p1, width=6, height=4)

p2 <- ggplot(prop, aes(x=time, y=value, fill=time)) + 
  geom_boxplot(outlier.shape = NA, lwd = 0.25, width = 0.8, position = position_dodge2(padding = 0.2)) +
  geom_point(aes(shape=Sex), size = 0.9, color = "red") +
  geom_line(aes(group=subject), size = 0.3) +
  facet_grid(group_name+Sex~cell_type) +
  labs(x = "Cell type", y = "Relative proportion", fill = "time") +
  scale_fill_manual(values = custom.colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", paired = TRUE, hide.ns = TRUE)

ggsave(paste0(path_to_project_file, "/Plots/training/deconvolution/proportion/proportion_by_small_groups_with_lines.pdf"), p2, width=6, height=6)
