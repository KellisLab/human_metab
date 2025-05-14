# Import libraries
library(org.Hs.eg.db)
library(fgsea)
library(dorothea)


# CLEAN DATA => COMPARISONS
clean_data_for_analysis <- function(data, comparison, coldata, baseline_training){
  rna_data <- data
  
  clean_coldata <- coldata
  clean_coldata[ , "comparison"] <- NA
  
  if(baseline_training == "baseline"){
    for(i in 1:length(comparison)){
      if(comparison[i] %in% c("obese", "lean")){
        clean_coldata[clean_coldata$large_group == comparison[i], ]$comparison <- comparison[i]
      }else if(comparison[i] %in% c("T2D obese", "non-diabetic obese", "MIT lean", "HIIT lean")){
        clean_coldata[clean_coldata$group_name == comparison[i], ]$comparison <- comparison[i]
      }else if(comparison[i] %in% c("male", "female", "m", "f")){
        if(comparison[i] == "male" | comparison[i] == "m"){
          clean_coldata[clean_coldata$Sex == "m", ]$comparison <- comparison[i]
        }else if(comparison[i] == "female" | comparison[i] == "f"){
          clean_coldata[clean_coldata$Sex == "f", ]$comparison <- comparison[i]
        }
      }else{
        print(paste0("comparison => ", comparison[i], " !!!"))
        stop("ERROR: unknown comparison or not found in coldata column names")
      }
    }
    clean_coldata <- clean_coldata[!is.na(clean_coldata$comparison) & clean_coldata$time == "A",]
    clean_data <- rna_data[, rownames(clean_coldata)]
  }else if(baseline_training == "training"){
    clean_data <- rna_data
    clean_coldata[ , "comparison"] <- clean_coldata$time
  }else{
    stop("ERROR: kati den ekana sosta")
  }
  
  return(list(clean_data, clean_coldata))
}


# SELECT INDIVIDUALS
select_individuals <- function(selection, clean_data, clean_coldata){
  for(i in 1:length(selection)){
    print(selection[i])
    if(selection[i] == "all"){
      return(list(clean_data, clean_coldata)) 
    }
    
    if(selection[i] %in% clean_coldata$Sex | selection[i] %in% c("male", "female")){
      if(selection[i] %in% clean_coldata$Sex){
        clean_coldata <- clean_coldata[clean_coldata$Sex == selection[i], ]
      }else{
        if(selection[i] == "male"){
          clean_coldata <- clean_coldata[clean_coldata$Sex == "m", ]
        }else{
          clean_coldata <- clean_coldata[clean_coldata$Sex == "f", ]
        }
      }
      print(paste0(i, " ", rownames(clean_coldata)))
    }else if(selection[i] %in% clean_coldata$large_group){
      clean_coldata <- clean_coldata[clean_coldata$large_group == selection[i], ]
    }else if(selection[i] %in% clean_coldata$group_name){
      clean_coldata <- clean_coldata[clean_coldata$group_name == selection[i], ]
      print(paste0(i, " ", rownames(clean_coldata)))
    }else if(selection[i] %in% clean_coldata$sex_large_group){
      clean_coldata <- clean_coldata[clean_coldata$sex_large_group == selection[i], ]
    }else{
      stop("ERROR: selection is not known")
    }
    
  }
  clean_data <- clean_data[, rownames(clean_coldata)]
  return(list(clean_data, clean_coldata)) 
}


### WILCOXON RANK SUM TEST for DEG analysis
wilcoxon_stats <- function(data, coldata, comparison, baseline_training){
  
  sampleSize <- ncol(data)
  data$gene.name <- row.names(data)
  data$p.value <- NA
  
  if(baseline_training == "baseline"){
    wilcox_paired=FALSE
  }else if(baseline_training == "training"){
    wilcox_paired=TRUE
  }else{
    stop("ERROR: kati den ekana sosta")
  }
  
  for(i.row in 1:nrow(data)){
    x <- as.numeric(data[i.row, rownames(coldata[coldata$comparison == comparison[1], ])])
    y <- as.numeric(data[i.row, rownames(coldata[coldata$comparison == comparison[2], ])])
    data[i.row, "p.value"] <- wilcox.test(x, y, paired=wilcox_paired)$p.value
  }
  data[, paste0("mean.", comparison[1])] <- rowMeans(data[, rownames(coldata[coldata$comparison == comparison[1], ])])
  data[, paste0("mean.", comparison[2])] <- rowMeans(data[, rownames(coldata[coldata$comparison == comparison[2], ])])
  data$FC <- data[, paste0("mean.", comparison[1])]/data[, paste0("mean.", comparison[2])]
  
  data$logFC <- log2(data$FC)
  data <- data[order(data$p.value),]
  data$min.log.pval <- -log10(data$p.value)
  data$p.adj <- p.adjust(data$p.value, method="fdr")
  
  data <- data[!is.infinite(data$logFC), ]
  data <- data[!is.na(data$logFC), ]
  data <- data[!is.na(data$p.value), ]
  
  return(data)
}


# DROP ZEROS AFTER DIFFERENTIAL EXPRESSION ANALYSIS
drop_rows_with_zeros <- function(clean_data, clean_coldata){
  clean_data$zeros <- apply(clean_data[, rownames(clean_coldata)], 1, function(x) sum(as.numeric(x) == 0))
  clean_data$zeros <- clean_data$zeros/nrow(clean_coldata)
  if(nrow(clean_data[clean_data$zeros > 0.3, ]) == 0){
    clean_data <- clean_data[, colnames(clean_data) != "zeros"]
    return(clean_data)
  }
  clean_data[clean_data$zeros > 0.3, ]$p.value = 1
  clean_data[clean_data$zeros > 0.3, ]$min.log.pval = 0
  clean_data <- clean_data[, colnames(clean_data) != "zeros"]
  clean_data <- clean_data[order(clean_data$p.value),]
  
  return(clean_data)
}


### PATHWAY ANALYSIS
pathway_analysis <- function(data, coldata, gene_set, path_to_project_file){
  path_to_intial_data <- paste0(path_to_project_file, "/Data")
  
  data$EntrezID <- mapIds(org.Hs.eg.db, toupper(data$gene.name), 'ENTREZID', 'SYMBOL')
  # Drop NA EntrezIDs
  data <- data[!is.na(data$EntrezID), ]
  
  ranks <- -log10(data$p.value) * sign(data$logFC)
  names(ranks) <- data$EntrezID
  
  if(gene_set == "CP"){
    pathways <- gmtPathways(paste0(path_to_intial_data, "/MSigDB gene sets/c2.cp.v7.5.entrez.gmt.txt"))
  }else if(gene_set == "REACTOME"){
    pathways <- gmtPathways(paste0(path_to_intial_data, "/MSigDB gene sets/c2.cp.reactome.v7.5.entrez.gmt.txt"))
  }else if(gene_set == "KEGG"){
    pathways <- gmtPathways(paste0(path_to_intial_data, "/MSigDB gene sets/c2.cp.kegg.v7.5.entrez.gmt.txt"))
  }else if(gene_set == "GO"){
    pathways <- gmtPathways(paste0(path_to_intial_data, "/MSigDB gene sets/c5.go.v7.5.entrez.gmt.txt"))
  }else{
    stop("ERROR: gene_set is neither 'CP', 'REACTOME', 'KEGG', nor 'GO'")
  }
  #fgseaRes <- fgsea(pathways, ranks, nperm=10000, minSize=15, maxSize=500)
  #fgseaRes <- fgsea(pathways, ranks, nPermSimple=10000, minSize=15)
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    eps      = 0.0,
                    minSize  = 15,
                    maxSize  = 500)
  
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][pval < 0.05], pathways, ranks)
  pathwayEnrichment <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), ]
  # pathwayEnrichment <- as.data.frame(fgseaRes[order(pval)][pval < 0.05])
  if(nrow(pathwayEnrichment) > 0){
    for(i.row in 1:nrow(pathwayEnrichment)){
      pathwayEnrichment[i.row, "leadingEdge"] <- paste(mapIds(org.Hs.eg.db, pathwayEnrichment[[i.row, "leadingEdge"]], 'SYMBOL', 'ENTREZID'), collapse=", ")
    }
  }
  
  return(pathwayEnrichment)
}


### GENE REGULATORY NETWORKS (GRN) ANALYSIS
TF_activity_with_Dorothea <- function(data, coldata){
  
  ## We load Dorothea Regulons
  data(dorothea_hs, package = "dorothea")
  regulons <- dorothea_hs %>%
    dplyr::filter(confidence %in% c("A", "B"))
  
  # We modify the format of the data table to make it suitable for running DoRothEA
  signature <- data[, names(data) %in% c("min.log.pval", "logFC")]
  signature$sign_logFC.min_log_pval <- signature$min.log.pval * sign(signature$logFC)
  signature <- signature[, names(signature) %in% c("sign_logFC.min_log_pval"), drop=FALSE]
  
  tf_activity <- dorothea::run_viper(signature, regulons,
                                     options =  list(minsize = 5, eset.filter = FALSE, 
                                                     cores = 1, verbose = FALSE, nes = TRUE))
  
  #comp1 <- as.matrix(data[, rownames(clean_coldata[clean_coldata$comparison == comparison[1],])])
  #comp2 <- as.matrix(data[, rownames(clean_coldata[clean_coldata$comparison == comparison[2],])])
  #dnull <- viper::ttestNull(comp1, comp2, per=100)
  #tf_activity_test <- msviper(ges=signature, regulon=regulon, nullmodel=dnull, verbose = FALSE)
  
  colnames(tf_activity) <- c("NES")
  tf_activity <- as.data.frame(tf_activity)
  tf_activity$TF <- rownames(tf_activity)
  tf_activity <- tf_activity[order(abs(tf_activity[,"NES"]), decreasing = TRUE),]
  for (i.row in 1:nrow(tf_activity)){
    target_list <- regulons$target[regulons$tf == tf_activity[i.row,"TF"]]
    target_list_up <- target_list[target_list %in% rownames(data[data$p.value < 0.1 & data$logFC > 0,])]
    target_list_down <- target_list[target_list %in% rownames(data[data$p.value < 0.1 & data$logFC < 0,])]
    tf_activity[i.row, "targets_up"] <- paste(as.vector(target_list_up),collapse= ", ")
    tf_activity[i.row, "targets_down"] <- paste(as.vector(target_list_down),collapse= ", ")
  }
  
  return(tf_activity)
}


# CUSTOM COLORS & SHAPES
custom.colours <- c("#E69EB3", "#F16A77", "#F9FF69", "#3CA28F", "#B7F385", "#AA99D6", "#FFCB89",
                    "#FF4747", "#5C7AFF", "#FF4747", "#5C7AFF", 
                    "#6FEAEC", "#DE3F8C", "#FF7733", "#73D260",
                    "#7681B3", "#7681B3",
                    "#F8766D", "#C77CFF", "#7CAE00", "#00BFC4", "#E6AB02", "#A6761D")
names(custom.colours) <- c("Adipocyte", "Adipose Stem Cell", "Fibroblast", "Myeloid", "Lymphoid", "Endothelial", "Smooth muscle",
                           "Up", "Down", "upregulated", "downregulated",
                           "Male", "Female", "obese", "lean", 
                           "Male + Female", "Lean + Obese",
                           "T2D obese", "non-diabetic obese", "MIT lean", "HIIT lean", "A", "B")
custom.shapes <- c(16, 17)
names(custom.shapes) <- c("f", "m")
