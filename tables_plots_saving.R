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

# Import libraries
library(xlsx)

# import sex specific genes detected at baseline
sex_specific_gene <- as.data.frame(read.xlsx(paste0(path_to_project_file, "/figures/supp_tables/TableS1_baseline_sex_specific_genes.xlsx"), sheetIndex = 1))$gene.name

# Set what and how to analyze
save_to_excel = TRUE
save_plots = TRUE
bulk_deconvolution <- c("bulk") # bulk or deconvolution
wilcoxon_DESeq2 <- c("wilcoxon") # wilcoxon or DESeq2
counts_tpm <- c("tpm")


### Create tables
if(bulk_deconvolution == "bulk"){
  
  # Select cts or tpm dataframe
  if(counts_tpm == "cts"){
    print("--> COUNTS")
    data <- cts
  }else if(counts_tpm == "tpm"){
    print("--> TPM")
    data <- tpm
  }else{
    stop("Input 'counts_tpm' should be either 'cts' or 'tpm")
  }
  
  for(baseline_training in c("baseline", "training")){
    
    if(baseline_training == "baseline"){
      comparisons <- list(c("obese", "lean"),
                          c("T2D obese", "lean"),
                          c("T2D obese", "non-diabetic obese"),
                          c("non-diabetic obese", "lean"),
                          c("male", "female"))
    }else{
      comparisons <- list(c("B", "A"))
    }
    
    # For every possible comparison (e.g. obese vs lean)
    for(i in 1:length(comparisons)){
      comparison <- comparisons[[i]]
      
      # Baseline obese vs lean + other similar configurations (e.g. obese T2D vs lean OR obese T2D vs non-diabetic obese)
      if(comparison[1] %in% c("obese", "non-diabetic obese", "T2D obese", "lean", "MIT lean", "HIIT lean")){
        selections <- list(c("all"), c("m"), c("f"))
      }
      # Baseline male vs female
      if(comparison[1] %in% c("m", "male", "f", "female")){
        selections <- list(c("all"), c("obese"), c("lean"))
      }
      # Training
      if(comparison[1] %in% c("A", "B")){
        selections <- list(c("all"), c("obese"), c("obese", "female"), c("obese", "male"),
                           c("lean"), c("lean", "female"), c("lean", "male"),
                           c("female"), c("male"),
                           c("non-diabetic obese"), c("non-diabetic obese", "female"), c("non-diabetic obese", "male"),
                           c("T2D obese"), c("T2D obese", "female"), c("T2D obese", "male"),
                           c("HIIT lean"), c("HIIT lean", "female"), c("HIIT lean", "male"),
                           c("MIT lean"), c("MIT lean", "female"))
      }
      
      # Save tables to excel
      if(save_to_excel == TRUE){
        wb.DEG <- xlsx::createWorkbook()
        wb.pathway <- xlsx::createWorkbook()
        wb.TF_activity <- xlsx::createWorkbook()
      }
      
      # Select specific subjects
      for(i in 1:length(selections)){
        
        selection <- selections[[i]]
        
        print(paste0(baseline_training, "/", comparison[1], " vs ", comparison[2], 
                     "/", selection))
        
        # CLEAN data and metadata
        res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
        clean_data <- res[[1]]
        clean_coldata <- res[[2]]
        # SELECT specific individuals for the comparison
        res <- select_individuals(selection, clean_data, clean_coldata)
        clean_data <- res[[1]]
        clean_coldata <- res[[2]]
        # WILCOXON RANK SUM TEST => DEGs
        clean_data <- wilcoxon_stats(clean_data, clean_coldata, comparison, baseline_training)
        # significant genes with many expression levels equal to zero as are no longer significant
        clean_data <- drop_rows_with_zeros(clean_data, clean_coldata)
        clean_data$baseline_sex_specific <- ifelse(clean_data$gene.name %in% sex_specific_gene, "Y", "")
        # PATHWAY analysis
        pathways <- data.frame()
        for (temp.gene_set in c("CP", "REACTOME", "KEGG", "GO")) {
          temp.pathways <- pathway_analysis(clean_data, clean_coldata, temp.gene_set, path_to_project_file)
          pathways <- rbind(pathways, temp.pathways)
        }
        ### GENE REGULATORY NETWORKS (GRN) analysis
        tf_activity <- TF_activity_with_Dorothea(clean_data, clean_coldata)
        
        if(save_plots == TRUE){
          if(baseline_training == "baseline"){
            if(length(selection) > 1){
              path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                     "/", comparison[1], " vs ", comparison[2], "/", bulk_deconvolution, 
                                     "/", wilcoxon_DESeq2, "/", selection[1], "_", selection[2])
            }else{
              path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                     "/", comparison[1], " vs ", comparison[2], "/", bulk_deconvolution, 
                                     "/", wilcoxon_DESeq2, "/", selection)
            }
          }else{
            if(length(selection) > 1){
              path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                     "/", bulk_deconvolution, "/", wilcoxon_DESeq2,
                                     "/", selection[1], "_", selection[2])
            }else{
              path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                     "/", bulk_deconvolution, "/", wilcoxon_DESeq2,
                                     "/", selection)
            }
          }
          
          dir.create(file.path(path_to_plot), showWarnings = FALSE)
          # plot DEG Volcano plot
          DEG_volcano_plot(clean_data, nb_labeled_genes=20, pval_threshold=0.05, path_to_plot = path_to_plot, logFC_threshold = 0.5, x_axis_cutoff = 3)
          # plot DEG MA plot
          DEG_MA_plot(clean_data, clean_coldata, label_custom_genes=FALSE, nb_labeled_genes=20, pval_threshold=0.05, path_to_plot)
          # plot pathways Ridge plots
          pathways_ridge_plot(clean_data, clean_coldata, pathways, nb_pathways=15, path_to_plot)
          # plot Dorothea TF activity
          TF_bar_plot(tf_activity, clean_data, nb_TF=20, path_to_plot)
          
        }
        
        
        if(save_to_excel == TRUE){
          if(length(selection) > 1){
            sheetName = paste0(selection[1], " ", selection[2])
          }else{
            sheetName = selection
          }
          
          # DEGs to excel file
          DEG.sheet <- xlsx::createSheet(wb.DEG, sheetName=sheetName)
          xlsx::addDataFrame(clean_data[clean_data$p.value < 0.05,], sheet=DEG.sheet, startColumn=1, row.names=TRUE, col.names=TRUE)
          # Pathways to excel file
          Pathway.sheet <- xlsx::createSheet(wb.pathway, sheetName=sheetName)
          xlsx::addDataFrame(pathways, sheet=Pathway.sheet, startColumn=1, row.names=FALSE, col.names=TRUE)
          # Dorothea TF activity (GRN) to excel file
          TF.sheet <- xlsx::createSheet(wb.TF_activity, sheetName=sheetName)
          xlsx::addDataFrame(tf_activity, sheet=TF.sheet, startColumn=1, row.names=TRUE, col.names=TRUE)
        }
        rm(res, clean_data, clean_coldata, pathways, tf_activity)
      }
      if(save_to_excel == TRUE){
        if(baseline_training == "baseline"){
          path_to_excel_table <- paste0(path_to_project_file, "/Data Analysis/", baseline_training, 
                                        "/", comparison[1], " vs ", comparison[2], "/", bulk_deconvolution, 
                                        "/", wilcoxon_DESeq2)
        }else{
          path_to_excel_table <- paste0(path_to_project_file, "/Data Analysis/", baseline_training, 
                                        "/", bulk_deconvolution, "/", wilcoxon_DESeq2)
        }
        xlsx::saveWorkbook(wb.DEG, paste0(path_to_excel_table, "/DEGs.xlsx"))
        xlsx::saveWorkbook(wb.pathway, paste0(path_to_excel_table, "/Pathways.xlsx"))
        xlsx::saveWorkbook(wb.TF_activity, paste0(path_to_excel_table, "/TF_activity_Dorothea.xlsx"))
        rm(wb.DEG, wb.pathway, wb.TF_activity)
      }
      # END COMPARISON LOOP
    }
    # END BASELINE-TRAINING LOOP
  }
  # END BULK
### TO DO 
}else if(bulk_deconvolution == "deconvolution"){
  cell_types <- c("Adipocyte", "ASC", "Fibroblast", "Myeloidcell", "Lymphoidcell", "EC", "SM")
  cell_types.full_anno <- c("Adipocyte", "Adipose Stem Cell", "Fibroblast", "Myeloid", "Lymphoid", "Endothelial", "Smooth muscle")
  
  cell.id <- 1
  # for each cell type
  for(cell in cell_types.full_anno){
    print(cell)
    path_to_intial_data <- paste0(path_to_project_file, "Data")
    data <- read.csv(paste0(path_to_intial_data, "/CIBERSORTx_Job19_output/CIBERSORTxHiRes_Job19_", cell_types[cell.id], "_Window28.txt"),
                     sep="\t")
    
    # pre-clean raw CIBERSORT data
    data <- data[!is.na(data$A1) & data$A1!=1,]
    rownames(data) <- data$GeneSymbol
    data = data[, !names(data) %in% c("GeneSymbol")]
    
    
    for(baseline_training in c("baseline", "training")){
      
      if(baseline_training == "baseline"){
        comparisons <- list(c("obese", "lean"),
                            c("T2D obese", "lean"),
                            c("T2D obese", "non-diabetic obese"),
                            c("non-diabetic obese", "lean"),
                            c("male", "female"))
      }else{
        comparisons <- list(c("B", "A"))
      }
      
      # For every possible comparison (e.g. obese vs lean)
      for(i in 1:length(comparisons)){
        comparison <- comparisons[[i]]
        
        # Baseline obese vs lean + other similar configurations (e.g. obese T2D vs lean OR obese T2D vs non-diabetic obese)
        if(comparison[1] %in% c("obese", "non-diabetic obese", "T2D obese", "lean", "MIT lean", "HIIT lean")){
          selections <- list(c("all"), c("m"), c("f"))
        }
        # Baseline male vs female
        if(comparison[1] %in% c("m", "male", "f", "female")){
          selections <- list(c("all"), c("obese"), c("lean"))
        }
        # Training
        if(comparison[1] %in% c("A", "B")){
          selections <- list(c("all"), c("obese"), c("obese", "female"), c("obese", "male"),
                             c("lean"), c("lean", "female"), c("lean", "male"),
                             c("female"), c("male"),
                             c("non-diabetic obese"), c("non-diabetic obese", "female"), c("non-diabetic obese", "male"),
                             c("T2D obese"), c("T2D obese", "female"), c("T2D obese", "male"),
                             c("HIIT lean"), c("HIIT lean", "female"), c("HIIT lean", "male"),
                             c("MIT lean"), c("MIT lean", "female"))
        }
        
        # Save tables to excel
        if(save_to_excel == TRUE){
          # wb.DEG <- xlsx::createWorkbook()
          # wb.pathway <- xlsx::createWorkbook()
          wb.TF_activity <- xlsx::createWorkbook()
        }
        
        # Select specific subjects
        for(i in 1:length(selections)){
          
          selection <- selections[[i]]
          
          print(paste0(baseline_training, "/", comparison[1], " vs ", comparison[2], 
                       "/", selection))
          
          # CLEAN data and metadata
          res <- clean_data_for_analysis(data, comparison, coldata, baseline_training)
          clean_data <- res[[1]]
          clean_coldata <- res[[2]]
          # SELECT specific individuals for the comparison
          res <- select_individuals(selection, clean_data, clean_coldata)
          clean_data <- res[[1]]
          clean_coldata <- res[[2]]
          # WILCOXON RANK SUM TEST => DEGs
          clean_data <- wilcoxon_stats(clean_data, clean_coldata, comparison, baseline_training)
          # significant genes with many expression levels equal to zero as are no longer significant
          clean_data <- drop_rows_with_zeros(clean_data, clean_coldata)
          # PATHWAY analysis
          # gene_set <- c("CP") # CP, REACTOME, KEGG, GO
          # pathways <- pathway_analysis(clean_data, clean_coldata, gene_set, path_to_project_file)
          ### GENE REGULATORY NETWORKS (GRN) analysis
          tf_activity <- TF_activity_with_Dorothea(clean_data, clean_coldata)
          
          if(save_plots == TRUE){
            if(baseline_training == "baseline"){
              if(length(selection) > 1){
                path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                       "/", comparison[1], " vs ", comparison[2], "/", bulk_deconvolution, 
                                       "/", wilcoxon_DESeq2, "/", cell)
                dir.create(file.path(path_to_plot), showWarnings = FALSE)
                path_to_plot <- paste0(path_to_plot, "/", selection[1], "_", selection[2])
              }else{
                path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                       "/", comparison[1], " vs ", comparison[2], "/", bulk_deconvolution, 
                                       "/", wilcoxon_DESeq2, "/", cell)
                dir.create(file.path(path_to_plot), showWarnings = FALSE)
                path_to_plot <- paste0(path_to_plot, "/", selection)
              }
            }else{
              if(length(selection) > 1){
                path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                       "/", bulk_deconvolution, "/", wilcoxon_DESeq2,
                                       "/", cell)
                dir.create(file.path(path_to_plot), showWarnings = FALSE)
                path_to_plot <- paste0(path_to_plot, "/", selection[1], "_", selection[2])
              }else{
                path_to_plot <- paste0(path_to_project_file, "/Plots/", baseline_training, 
                                       "/", bulk_deconvolution, "/", wilcoxon_DESeq2,
                                       "/", cell)
                dir.create(file.path(path_to_plot), showWarnings = FALSE)
                path_to_plot <- paste0(path_to_plot, "/", selection)
              }
            }
            
            dir.create(file.path(path_to_plot), showWarnings = FALSE)
            # plot DEG Volcano plot
            DEG_volcano_plot(clean_data, label_custom_genes=clean_data[clean_data$gene.name %in% c(tf1,tf2,tf3,pathway1,pathway2,pathway3,global_ol_oppo) & clean_data$p.value < 0.05, "gene.name"], nb_labeled_genes=0, pval_threshold=0.05, path_to_plot, logFC_threshold = 0, x_axis_cutoff = 0.5)
            # plot DEG MA plot
            # DEG_MA_plot(clean_data, clean_coldata, label_custom_genes=FALSE, nb_labeled_genes=20, pval_threshold=0.05, path_to_plot)
            # plot pathways Ridge plots
            # pathways_ridge_plot(clean_data, clean_coldata, pathways, nb_pathways=15, path_to_plot)
            # plot Dorothea TF activity
            TF_bar_plot(tf_activity, clean_data, nb_TF=20, path_to_plot)
          }
          
          if(save_to_excel == TRUE){
            if(length(selection) > 1){
              sheetName = paste0(selection[1], " ", selection[2])
            }else{
              sheetName = selection
            }
            
            # DEGs to excel file
            # DEG.sheet <- xlsx::createSheet(wb.DEG, sheetName=sheetName)
            # xlsx::addDataFrame(clean_data[clean_data$p.value < 0.05,], sheet=DEG.sheet, startColumn=1, row.names=TRUE, col.names=TRUE)
            # Pathways to excel file
            # Pathway.sheet <- xlsx::createSheet(wb.pathway, sheetName=sheetName)
            # xlsx::addDataFrame(pathways, sheet=Pathway.sheet, startColumn=1, row.names=FALSE, col.names=TRUE)
            # Dorothea TF activity (GRN) to excel file
            TF.sheet <- xlsx::createSheet(wb.TF_activity, sheetName=sheetName)
            xlsx::addDataFrame(tf_activity, sheet=TF.sheet, startColumn=1, row.names=TRUE, col.names=TRUE)
          }
          rm(res, clean_data, clean_coldata, pathways, tf_activity)
        }
        
        if(save_to_excel == TRUE){
          if(baseline_training == "baseline"){
            path_to_excel_table <- paste0(path_to_project_file, "/Data Analysis/", baseline_training, 
                                          "/", comparison[1], " vs ", comparison[2], "/", bulk_deconvolution, 
                                          "/", wilcoxon_DESeq2,
                                          "/", cell)
            dir.create(file.path(path_to_excel_table), showWarnings = FALSE)
          }else{
            path_to_excel_table <- paste0(path_to_project_file, "/Data Analysis/", baseline_training, 
                                          "/", bulk_deconvolution, "/", wilcoxon_DESeq2,
                                          "/", cell)
            dir.create(file.path(path_to_excel_table), showWarnings = FALSE)
          }
          # xlsx::saveWorkbook(wb.DEG, paste0(path_to_excel_table, "/DEGs.xlsx"))
          # xlsx::saveWorkbook(wb.pathway, paste0(path_to_excel_table, "/Pathways.xlsx"))
          xlsx::saveWorkbook(wb.TF_activity, paste0(path_to_excel_table, "/TF_activity_Dorothea.xlsx"))
          rm(wb.DEG, wb.pathway, wb.TF_activity)
        }
        # END COMPARISON LOOP
      }
      # END BASELINE-TRAINING LOOP
    }
    # END CELL DECONVOLUTION
    cell.id <- cell.id + 1
  }
}

