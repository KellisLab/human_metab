library(fmsb)

# define path
# path_to_project_file <- "/Users/Danae/Documents/MIT Internship/Jan-Willem Human Study"
path_to_project_file <- "~/Dropbox (MIT)/metabolism_Joslin/Jan_Willem/Jan-Willem Human Study_RNAseq"

# run data_preprocessing.R
source(paste0(path_to_project_file, "/Data Analysis R Files/data_preprocessing.R"))

# run data_analysis_helper.R
source(paste0(path_to_project_file, "/Data Analysis R Files/data_analysis_helper.R"))

# run plotting.R
source(paste0(path_to_project_file, "/Data Analysis R Files/plotting.R"))

# set working directory
setwd(path_to_project_file)

# Set what and how to analyze
wilcoxon_DESeq2 <- c("wilcoxon") # wilcoxon or DESeq2
counts_tpm <- c("tpm")
comparison <- c("obese", "lean")
# comparison <- c("T2D obese", "non-diabetic obese")
# comparison <- c("A", "B")
selection <- c("all")
# selection <- c("obese", "lean")
baseline_training <- c("baseline")
# baseline_training <- c("training")

# Get the gene sets (this step is for baseline, the gene sets for training are in the "global_ol_oppo" object defined in PCA_plot and the "pathway1/2/3" and "tf1/2/3" objects)
bulk_deconvolution <- c("bulk") # bulk or deconvolution
# get path to excel file
if(baseline_training == "baseline"){
  path_to_folder <- paste0(path_to_project_file, "/Data Analysis", 
                           "/", baseline_training, 
                           "/", comparison[1], " vs ", comparison[2], 
                           "/", bulk_deconvolution, 
                           "/", wilcoxon_DESeq2)
}else{
  path_to_folder <- paste0(path_to_project_file, "/Data Analysis", 
                           "/", baseline_training, 
                           "/", bulk_deconvolution,
                           "/", wilcoxon_DESeq2)
}


# read excel file containing the pathways
if(length(selection) == 1){
  bulk_pathways <- as.data.frame(read_excel(paste0(path_to_folder, "/Pathways.xlsx"), sheet = selection, col_names = TRUE))
}else{
  bulk_pathways <- as.data.frame(read_excel(paste0(path_to_folder, "/Pathways.xlsx"), sheet = paste0(selection[1], " ", selection[2]), col_names = TRUE))
}



# select a group of pathways
# for Figure 1F
up_down = "up"
immune_group_name = "immune_related"
immune_select_pathways <- c("REACTOME_NEUTROPHIL_DEGRANULATION",
                     "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
                     "REACTOME_INTERFERON_SIGNALING",
                     "REACTOME_INTERFERON_GAMMA_SIGNALING")
ecm_group_name = "ECM_related"
ecm_select_pathways <- c("NABA_CORE_MATRISOME",
                     "WP_BURN_WOUND_HEALING",
                     "REACTOME_INTERFERON_SIGNALING",
                     "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                     "NABA_PROTEOGLYCANS"
                     )

up_down = "down"
tca_group_name = "TCA_ATP_respiratory_related"
tca_select_pathways <- c("REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS",
                     "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                     "KEGG_CITRATE_CYCLE_TCA_CYCLE"
                     )
fa_group_name = "fatty_acid_related"
fa_select_pathways <- c("KEGG_FATTY_ACID_METABOLISM",
                     "WP_FATTY_ACID_BIOSYNTHESIS",
                     "WP_FATTY_ACID_BETAOXIDATION",
                     "WP_OXIDATIVE_PHOSPHORYLATION"
)

# for Figure 2F (weak signal, we do not have this figure panel anymore)
up_down = "up"
group_name = "ribo_related"
select_pathways <- c("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
                     "WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS",
                     "KEGG_RIBOSOME",
                     "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
                     "REACTOME_SELENOAMINO_ACID_METABOLISM",
                     "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
                     "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
                     "REACTOME_INFLUENZA_INFECTION")
group_name = "ECM_related"
select_pathways <- c("NABA_CORE_MATRISOME",
                     "NABA_PROTEOGLYCANS"
)
group_name = "necrosis"
select_pathways <- c("REACTOME_REGULATED_NECROSIS")
group_name = "DNA_damage"
select_pathways <- c("REACTOME_DNA_DAMAGE_BYPASS")

up_down = "down"
group_name = "hormone"
select_pathways <- c("WP_CORTICOTROPINRELEASING_HORMONE_SIGNALING_PATHWAY")
                     
# get the leading edges for this group of pathways
for (temp in c("immune", "ecm", "tca", "fa")) {
  temp.bulk_genes <- bulk_pathways[bulk_pathways$pathway %in% get(paste0(temp, "_select_pathways")), "leadingEdge"]
  
  # from a list of genes to a single vector
  for(i in 1:length(temp.bulk_genes)){
    temp.bulk_genes[i] <- strsplit(temp.bulk_genes[[i]], ", ")
  }
  temp.bulk_genes <- unique(unlist(temp.bulk_genes))
  
  assign(paste0(temp, "_bulk_genes"), temp.bulk_genes)
}

# check the number of these genes that are significant in each cell type
bulk_deconvolution <- c("deconvolution")
cells <- c("Adipocyte", "Adipose Stem Cell", "Fibroblast", "Myeloid", "Lymphoid", "Endothelial", "Smooth muscle")
# create data frame for the spider web plot
spider_input <- function(path_to_project_file, baseline_training, comparison, bulk_deconvolution, wilcoxon_DESeq2, selection, gene_set, cells) {
  data <- as.data.frame(matrix(0, 1, length(cells)))
  colnames(data) <- cells
  for(cell in cells){
    print(cell)
    
    # get path to excel file
    if(baseline_training == "baseline"){
      path_to_folder <- paste0(path_to_project_file, "/Data Analysis", 
                               "/", baseline_training, 
                               "/", comparison[1], " vs ", comparison[2], 
                               "/", bulk_deconvolution, 
                               "/", wilcoxon_DESeq2,
                               "/", cell)
    }else{
      path_to_folder <- paste0(path_to_project_file, "/Data Analysis", 
                               "/", baseline_training, 
                               "/", bulk_deconvolution,
                               "/", wilcoxon_DESeq2,
                               "/", cell)
    }
  
    # read excel file containing the DEGs
    for (i in 1:length(selection)) {
      deconvoluted_DEGs <- as.data.frame(read_excel(paste0(path_to_folder, "/DEGs.xlsx"), sheet = selection[i], col_names = TRUE))
      # drop the first column "..1"
      deconvoluted_DEGs[1] <- NULL
      rownames(deconvoluted_DEGs) <- deconvoluted_DEGs$gene.name
      # number of genes in the gene set that are significant in this cell type
      nb_sig_genes <- nrow(deconvoluted_DEGs[deconvoluted_DEGs$gene.name %in% gene_set,])
      
      data[selection[i], cell] <- nb_sig_genes
    }
  }
  return(data)
}
for (temp in c("immune_bulk_genes", "ecm_bulk_genes", "tca_bulk_genes", "fa_bulk_genes")) {
  temp.input <- spider_input(path_to_project_file, baseline_training, comparison, bulk_deconvolution, wilcoxon_DESeq2, selection, get(temp), cells)
  assign(paste0(temp, "_spider_input"), temp.input)
}

data <- rbind(ecm_bulk_genes_spider_input["all",], immune_bulk_genes_spider_input["all",], fa_bulk_genes_spider_input["all",], tca_bulk_genes_spider_input["all",])
data$group_name <- c(ecm_group_name, immune_group_name, fa_group_name, tca_group_name)
data$color <- c("red", "red", "blue", "blue")
data_long <- reshape2::melt(data)
if(up_down == "up"){
  spider_color = custom.colours[["Up"]]
}else if(up_down == "down"){
  spider_color = custom.colours[["Down"]]
}else{
  stop("up_down variable is neither 'up' or 'down")
}

if(baseline_training == "baseline"){
  path_to_plot <- paste0(path_to_project_file, "/Plots", 
                         "/", baseline_training, 
                         "/", comparison[1], " vs ", comparison[2], 
                         "/bulk-deconvolution comparison")
}else{
  path_to_plot <- paste0(path_to_project_file, "/Plots", 
                         "/", baseline_training, 
                         "/bulk-deconvolution comparison")
}

# spider web chart
pdf(paste0(path_to_plot, "/spider_web_chart_", group_name, ".pdf"), width = 4.5, height = 3.5)
spider_web_chart <- radarchart(data, axistype=1,
           #custom polygon
           pcol=spider_color, pfcol = scales::alpha(spider_color, 0.5), plwd=2,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", calcex=0.7, caxislabels=round(seq(0, length(bulk_genes), length = 5)), cglwd=0.8,
           #custom labels
           vlcex = 0.7)
dev.off()
# bar plot
data_long$group_name <- factor(data_long$group_name, levels = c("ECM_related", "immune_related", "fatty_acid_related", "TCA_ATP_respiratory_related"))
pdf(paste0(path_to_plot, "/bar_plot.pdf"), width = 3.5, height = 3.5)
ggplot(data_long, aes(variable, value)) +
  geom_bar(aes(fill = color), stat = "identity", width = 0.5) +
  facet_wrap(. ~ group_name, ncol = 2, scales = "free_y") +
  ylab("Number of cell-type-specific DEGs") +
  xlab("") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# plot the number of training DEGs for each cell type in obese and lean
deconvoluted_DEG_num <- data.frame()
for (i in 1:length(selection)) {
    for (temp.cell in cells) {
    temp.deconvoluted_DEGs <- as.data.frame(read_excel(paste0(paste0(path_to_project_file, "/Data Analysis/", baseline_training, "/", bulk_deconvolution, "/", wilcoxon_DESeq2, "/", temp.cell), "/DEGs.xlsx"), sheet = selection[i], col_names = TRUE))
    # drop the first column "..1"
    temp.deconvoluted_DEGs[1] <- NULL
    for (temp.dir in c("up", "down")) {
      if (temp.dir == "up") {
        temp.sig_genes <- sum(temp.deconvoluted_DEGs$logFC > 0)
      } else {
        temp.sig_genes <- sum(temp.deconvoluted_DEGs$logFC < 0)
      }
      deconvoluted_DEG_num[i, paste0(temp.cell, "_", temp.dir)] <- temp.sig_genes
    }
    assign(paste0("deconvoluted_DEGs_", selection[i], "_", temp.cell), temp.deconvoluted_DEGs)
  }
  deconvoluted_DEG_num[i, "selection"] <- selection[i]
}
rownames(deconvoluted_DEG_num) <- deconvoluted_DEG_num$selection
DEG_num_long <- as.data.frame(t(deconvoluted_DEG_num[, 1:14]))
DEG_num_long$cell_type <- gsub("_.*", "", rownames(DEG_num_long))
DEG_num_long$dir <- gsub(".*_", "", rownames(DEG_num_long))
DEG_num_long$obese_signed <- ifelse(DEG_num_long$dir == "down", -DEG_num_long$obese, DEG_num_long$obese)
DEG_num_long$lean_signed <- ifelse(DEG_num_long$dir == "down", -DEG_num_long$lean, DEG_num_long$lean)
DEG_num_long_plot <- melt(DEG_num_long[, -c(1,2)], id.vars = c("cell_type", "dir"), measure.vars = c("obese_signed", "lean_signed"), variable.name = "DEG_num")
DEG_num_long_plot$cell_type <- factor(DEG_num_long_plot$cell_type, levels = rev(cell_types.full_anno))
DEG_num_long_plot$DEG_num <- factor(DEG_num_long_plot$DEG_num, levels = c("lean_signed", "obese_signed"))
pdf("Plots/training/deconvolution/wilcoxon/DEG_num_summary.pdf", width = 6, height = 3)
ggplot(DEG_num_long_plot, aes(cell_type, value)) +
  geom_col(aes(fill = dir), width = 0.8) +
  coord_flip() +
  ylab("Number of DEGs") +
  xlab("") +
  scale_fill_manual(values = c("blue", "red")) +
  facet_grid(. ~ DEG_num) +
  theme_classic() +
  theme(text = element_text(size = 14))
dev.off()
