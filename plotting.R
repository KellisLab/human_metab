### PLOTS

library(ggplot2)
library(ggridges)
library(ggrepel)
library(stringr)
library(gridExtra)
library(ggnewscale)
library(patchwork)


# Empty plot 
empty_plot <- function(){
  
  empty_plot <- ggplot() + 
    theme_void() +
    geom_text(aes(0,0,label='N/A')) +
    xlab(NULL)
  
  return(empty_plot)
}

# Volcano plot
DEG_volcano_plot <- function(data, label_custom_genes=NULL, nb_labeled_genes=20, pval_threshold=0.05, path_to_plot, logFC_threshold=0.5, x_axis_cutoff=Inf){
  if (is.infinite(x_axis_cutoff)) {
    plotXlim <- max(abs(data$logFC))
  } else {
    plotXlim <- x_axis_cutoff
  }

  # select top upregulated and top downregulated genes
  if(is.null(label_custom_genes)){
    select_genes <- data[order(data$p.value, data$logFC), ]
    select_genes <- select_genes[select_genes$p.value<pval_threshold & select_genes$logFC>logFC_threshold, "gene.name"][1:(nb_labeled_genes/2)]
    select_downreg_genes <- data[order(data$p.value, -data$logFC), ]
    select_downreg_genes <- select_downreg_genes[select_downreg_genes$p.value<pval_threshold & select_downreg_genes$logFC<(-logFC_threshold), "gene.name"][1:(nb_labeled_genes/2)]
    select_genes <- append(select_genes, select_downreg_genes)
  }else{
    select_genes <- label_custom_genes
  }
  
  DEG_volcano_plot <- ggplot(data, aes(x=logFC , y=min.log.pval)) +
    labs(y = bquote(~-log[10]~ "p-value"), x = bquote(~log[2]~ "FC"),
         caption = as.expression(bquote("n"[sigDOWN]==.(nrow(data[data$p.value<pval_threshold & data$logFC< -logFC_threshold,]))~
                                          "  n"[sigUP]==.(nrow(data[data$p.value<pval_threshold & data$logFC> logFC_threshold,]))))) +
    xlim(c(-plotXlim, plotXlim)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(data$min.log.pval) + 0.2)) +
    geom_point(data = data[data$logFC > logFC_threshold & data$p.value<pval_threshold & data$logFC < x_axis_cutoff,],
               aes(x=logFC, y=min.log.pval), fill=custom.colours[["Up"]], color=custom.colours[["Up"]], alpha = 0.6, size=1) +
    geom_point(data = data[data$logFC > x_axis_cutoff & data$p.value<pval_threshold,],
               aes(x=x_axis_cutoff, y=min.log.pval), fill=custom.colours[["Up"]], color=custom.colours[["Up"]], alpha = 0.6, size=1, shape=17) +
    geom_point(data = data[data$logFC < (-logFC_threshold) & data$p.value<pval_threshold & data$logFC > (-x_axis_cutoff),],
               aes(x=logFC, y=min.log.pval), fill=custom.colours[["Down"]], color=custom.colours[["Down"]], alpha = 0.6, size=1) +
    geom_point(data = data[data$p.value<pval_threshold & data$logFC < (-x_axis_cutoff),],
               aes(x=(-x_axis_cutoff), y=min.log.pval), fill=custom.colours[["Down"]], color=custom.colours[["Down"]], alpha = 0.6, size=1, shape=17) +
    geom_point(data = data[abs(data$logFC)<x_axis_cutoff & data$p.value>pval_threshold,],
               aes(x=logFC, y=min.log.pval), shape=21, fill="grey40", color="grey45", alpha = 0.5, size=0.5) +
    geom_point(data = data[data$logFC>x_axis_cutoff & data$p.value>pval_threshold,],
               aes(x=x_axis_cutoff, y=min.log.pval), fill="grey40", color="grey45", alpha = 0.5, size=0.5, shape=17) +
    geom_point(data = data[data$logFC<(-x_axis_cutoff) & data$p.value>pval_threshold,],
               aes(x=(-x_axis_cutoff), y=min.log.pval), fill="grey40", color="grey45", alpha = 0.5, size=0.5, shape=17) +
    geom_text_repel(data = data[data$gene.name %in% select_genes, ],
                    aes(x=logFC, y=min.log.pval, label=gene.name), size=3, color="black", max.overlaps = getOption("ggrepel.max.overlaps", default = 30), min.segment.length = 0, seed = 42, box.padding = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold), size = 0.4, linetype = "longdash") +
    geom_vline(xintercept = logFC_threshold, size = 0.4, linetype = "longdash") +
    geom_vline(xintercept = -logFC_threshold, size = 0.4, linetype = "longdash") +
    theme_classic() +
    theme(text = element_text(size = 14))
  
  ggsave(paste0(path_to_plot, "/DEG_volcano_plot.pdf"), DEG_volcano_plot, width=4.5, height=4)
  
}

# MA plot

DEG_MA_plot <- function(data, coldata, label_custom_genes=FALSE, nb_labeled_genes=20, pval_threshold=0.05, path_to_plot, logFC_threshold=0.5){
  # add avgExp (average expression column)
  data$avg_exp <- rowMeans(data[, names(data) %in% rownames(coldata)])
  data$log.avg_exp <- log2(data$avg_exp)
  
  # select top upregulated and top downregulated genes
  if(label_custom_genes == FALSE){
    select_genes <- data[order(data$p.value, data$logFC), ]
    select_genes <- select_genes[select_genes$p.value<pval_threshold & abs(select_genes$logFC)>logFC_threshold, "gene.name"][1:(nb_labeled_genes/2)]
    select_downreg_genes <- data[order(data$p.value, -data$logFC), ]
    select_downreg_genes <- select_downreg_genes[select_downreg_genes$p.value<pval_threshold & abs(select_downreg_genes$logFC)>logFC_threshold, "gene.name"][1:(nb_labeled_genes/2)]
    select_genes <- append(select_genes, select_downreg_genes)
  }else{
    select_genes <- label_custom_genes
  }
  
  MA_plot <- ggplot(data, aes(x=log.avg_exp, y=logFC)) +
    labs(x = bquote(~log[2]~ "mean expression"), y = bquote(~log[2]~ "FC"),
         caption = as.expression(bquote("n"[sig_downreg_genes]==.(nrow(data[data$p.value<pval_threshold & data$logFC< -logFC_threshold,]))~
                                     "  n"[sig_upreg_genes]==.(nrow(data[data$p.value<pval_threshold & data$logFC> logFC_threshold,]))))) +
    geom_point(data = data[!(abs(data$logFC)>logFC_threshold & data$p.value<pval_threshold),],
               aes(x=log.avg_exp, y=logFC), shape=21, fill="grey40", color="grey45", alpha = 0.5, size=0.5) +
    geom_point(data = data[data$logFC > logFC_threshold & data$p.value<pval_threshold,],
               aes(x=log.avg_exp, y=logFC), fill=custom.colours[["Up"]], color=custom.colours[["Up"]], alpha = 0.6, size=1) +
    geom_point(data = data[data$logFC < -logFC_threshold & data$p.value<pval_threshold,],
               aes(x=log.avg_exp, y=logFC), fill=custom.colours[["Down"]], color=custom.colours[["Down"]], alpha = 0.6, size=1) +
    geom_text_repel(data = data[data$gene.name %in% select_genes, ],
                    aes(x=log.avg_exp, y=logFC, label=gene.name), size=3, color="black", max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
    geom_hline(yintercept = logFC_threshold, size = 0.3, linetype = "dashed") +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_hline(yintercept = -logFC_threshold, size = 0.3, linetype = "dashed") +
    theme_classic()
  
  ggsave(paste0(path_to_plot, "/DEG_MA_plot.pdf"), MA_plot, width=4.5, height=4)
}


### PATHWAYS RIDGE PLOTS
pathways_ridge_plot <- function(data, coldata, pathways, nb_pathways=15, path_to_plot, pval_threshold = 0.05){
  
  ranks <- -log10(data$p.value) * sign(data$logFC)
  names(ranks) <- data$gene.name
  
  pathways$pathway <- str_replace_all(pathways$pathway, "_", " ")
  top_upreg_pathways <- pathways[(pathways$pval < pval_threshold) & sign(pathways$NES)==1,]
  top_downreg_pathways <- pathways[(pathways$pval < pval_threshold) & sign(pathways$NES)==-1,]
  
  
  for(up_down in c("up", "down")){
    if(up_down == "up"){
      # if there are more than "nb_pathways" upreg pathways in the dataframe
      if(nrow(top_upreg_pathways) >= nb_pathways){
        top_upreg_pathways <- top_upreg_pathways[1:nb_pathways,]
      # if there are no upreg pathways
      }else if(nrow(top_upreg_pathways) == 0){
        ridge.plot <- empty_plot()
        ggsave(paste0(path_to_plot, "/top", nb_pathways, "_", up_down, "reg_pathways_ridge_plot.pdf"), ridge.plot, width=5, height=4)
        ridge.df.up <- NULL
        next
      }
      top_pathways <- top_upreg_pathways
    }else{
      # if there are more than "nb_pathways" downreg pathways in the dataframe
      if(nrow(top_downreg_pathways) >= nb_pathways){
        top_downreg_pathways <- top_downreg_pathways[1:nb_pathways,]
      # if there are no downreg pathways
      }else if(nrow(top_downreg_pathways) == 0){
        ridge.plot <- empty_plot()
        ggsave(paste0(path_to_plot, "/top", nb_pathways, "_", up_down, "reg_pathways_ridge_plot.pdf"), ridge.plot, width=5, height=4)
        ridge.df.down <- NULL
        next
      }
      top_pathways <- top_downreg_pathways
    }
    start <- TRUE
    for(i.row in 1:nrow(top_pathways)){
      leading_edges <- unlist(strsplit(top_pathways[[i.row, "leadingEdge"]], ", "))
      for(gene in leading_edges){
        if(start == TRUE){
          ridge.df <- data.frame(pathway=top_pathways[[i.row, "pathway"]],
                                 pval = top_pathways[[i.row, "pval"]],
                                 padj = top_pathways[[i.row, "padj"]],
                                 entrezID = gene,
                                 value = ranks[[gene]])
          start <- FALSE
        }else{
          ridge.df <- rbind(ridge.df, noquote(c(top_pathways[[i.row, "pathway"]], top_pathways[[i.row, "pval"]], 
                                                top_pathways[[i.row, "padj"]], gene, ranks[[gene]])))
        }
      }
    }
    ridge.df$pval <- as.numeric(ridge.df$pval)
    ridge.df$padj <- as.numeric(ridge.df$padj)
    ridge.df$value <- as.numeric(ridge.df$value)
    ridge.df$minus_log10pval <- -log10(ridge.df$pval)
    
    if(up_down == "up"){
      ridge.df.up <- ridge.df
    }else if(up_down == "down"){
      ridge.df.down <- ridge.df
    }else{
      stop("updown variable is neither 'up' or 'down'")
    }
  }
  
  # find the min and max pathway pvalue and the min and max gene ranking metric value
  if(is.data.frame(ridge.df.up) & is.data.frame(ridge.df.down)){
    min_pathway_minus_log10pval = min(min(ridge.df.up$minus_log10pval), min(ridge.df.down$minus_log10pval))
    max_pathway_minus_log10pval = max(max(ridge.df.up$minus_log10pval), max(ridge.df.down$minus_log10pval))
    min_gene_value = min(min(ridge.df.up$value), min(ridge.df.down$value))
    max_gene_value = max(max(ridge.df.up$value), max(ridge.df.down$value))
  # if there aren't any significant downregulated pathways
  }else if(is.data.frame(ridge.df.up)){
    min_pathway_minus_log10pval = min(ridge.df.up$minus_log10pval)
    max_pathway_minus_log10pval = max(ridge.df.up$minus_log10pval)
    min_gene_value = min(ridge.df.up$value)
    max_gene_value = max(ridge.df.up$value)
  # if there aren't any significant upregulated pathways
  }else if(is.data.frame(ridge.df.down)){
    min_pathway_minus_log10pval = min(ridge.df.down$minus_log10pval)
    max_pathway_minus_log10pval = max(ridge.df.down$minus_log10pval)
    min_gene_value = min(ridge.df.down$value)
    max_gene_value = max(ridge.df.down$value)
  # if there aren't any significant up- or down- regulated pathways
  }else{
    return(NULL)
  }
  
  abs_max_gene_value <- max(abs(min_gene_value), abs(max_gene_value))
  
  for(up_down in c("up", "down")){
    if(up_down == "up"){
      ridge.df <- ridge.df.up
      top_color = "#FF4747"
    }else if(up_down == "down"){
      ridge.df <- ridge.df.down
      top_color = "#5c7aff"
    }
    if(!is.data.frame(ridge.df)){
      next
    }
    ridge.plot <- ggplot(ridge.df, aes(x = value, y = reorder(pathway, minus_log10pval), fill=minus_log10pval)) + 
      geom_density_ridges(scale = 0.8, size=0.4, rel_min_height=0.001,
                          jittered_points = TRUE,
                          position = position_points_jitter(width = 0.05, height = 0),
                          point_shape = '|', point_size = 2, point_alpha = 0.4) +
      theme_classic()+
      theme(legend.position = "right",
            axis.line.x = element_line(size = 0.35),
            axis.ticks.y = element_blank(),
            axis.text.y=element_text(size=6, colour = "black"),
            axis.line.y = element_blank()) +
      scale_fill_gradient(low="white", high=top_color, limits = c(min_pathway_minus_log10pval-0.25, max_pathway_minus_log10pval+0.25))+
      scale_x_continuous(name="Genes -log10(p.value) * sign(logFC)") + 
      ylab("Pathways") +
      labs(fill = "Pathway log10(p-value)") +
      xlim(-abs_max_gene_value-0.1, abs_max_gene_value+0.1)
    
    if(up_down == "up"){
      gg_up <- ridge.plot + scale_y_discrete(labels = function(x) str_wrap(x, width = 50), position = "right")
      ggsave(paste0(path_to_plot, "/top", nb_pathways, "_", up_down, "reg_pathways_ridge_plot.pdf"), gg_up, width=7, height=4)
    }else if(up_down == "down"){
      gg_down <- ridge.plot + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
      ggsave(paste0(path_to_plot, "/top", nb_pathways, "_", up_down, "reg_pathways_ridge_plot.pdf"), gg_down, width=7, height=4)
    }
  }
  # if there are only downregulated pathways
  if(nrow(top_upreg_pathways) == 0){
    ggsave(paste0(path_to_plot, "/pathways_ridge_plot.pdf"), gg_down, width=7, height=4)
  # if there are only upregulated pathways
  }else if(nrow(top_downreg_pathways) == 0){
    ggsave(paste0(path_to_plot, "/pathways_ridge_plot.pdf"), gg_up, width=7, height=4)
  }else{
    plot <- gg_down + gg_up + plot_layout(ncol = 2)
    ggsave(paste0(path_to_plot, "/pathways_ridge_plot.pdf"), plot, width=15, height=7)
  }
}



### TF BAR PLOTS

TF_bar_plot <- function(tf_activity, clean_data, nb_TF=20, path_to_plot){
  
  top_upreg_TF <- tf_activity[sign(tf_activity$NES) == 1,]
  top_downreg_TF <- tf_activity[sign(tf_activity$NES) == -1,]
  
  # if there are more than "nb_TF" TF in the dataframe
  if(nrow(top_upreg_TF) >= nb_TF){
    top_upreg_TF <- top_upreg_TF[1:nb_TF,]
    max_NES <- max(top_upreg_TF$NES)
  }else if(nrow(top_upreg_TF) == 0){
    TF.plot <- empty_plot()
    ggsave(paste0(path_to_plot, "/Dorothea_upreg_TF_activity_barplot.pdf"), TF.plot, width=5, height=4)
    skip_up = 1
  }else{
    max_NES <- max(top_upreg_TF$NES)
  }
  if(nrow(top_downreg_TF) >= nb_TF){
    top_downreg_TF <- top_downreg_TF[1:nb_TF,]
    if(exists("max_NES")){
      max_NES <- max(max_NES, abs(min(top_downreg_TF$NES)))
    }else{
      max_NES <- abs(min(top_downreg_TF$NES))
    }
  }else if(nrow(top_downreg_TF) == 0){
    TF.plot <- empty_plot()
    ggsave(paste0(path_to_plot, "/Dorothea_downreg_TF_activity_barplot.pdf"), TF.plot, width=5, height=4)
    skip_down = 1
  }else{
    if(exists("max_NES")){
      max_NES <- max(max_NES, abs(min(top_downreg_TF$NES)))
    }else{
      max_NES <- abs(min(top_downreg_TF$NES))
    }
  }
  
  for(up_down in c("up", "down")){
    if(up_down == "up"){
      top_TF <- top_upreg_TF
      if(exists("skip_up")){
        next
      }
    }else{
      top_TF <- top_downreg_TF
      if(exists("skip_down")){
        next
      }
    }
  
    # top_TF <- merge(top_TF, clean_data[rownames(clean_data) %in% top_TF$TF, c("logFC", "min.log.pval")],
    #                 by.x = "TF", by.y = "row.names")
    if(nrow(top_TF) == 0){
      TF.plot <- empty_plot()
      ggsave(paste0(path_to_plot, "/Dorothea_", up_down, "reg_TF_activity_barplot.pdf"), TF.plot, width=5, height=4)
      next
    }
    # top_TF$sign_logFC.min_log_pval <- sign(top_TF$logFC) * top_TF$min.log.pval
    
    if(up_down == "up"){
      plot <- ggplot(top_TF,aes(x = NES, y = reorder(TF, NES))) +
        scale_y_discrete(position = "right")
    }else{
      plot <- ggplot(top_TF,aes(x = NES, y = reorder(TF, -NES)))
    }
    
    plot <- plot +
      geom_bar(aes(fill = NES), stat = "identity") +
      #geom_bar(aes(fill = sign_logFC.min_log_pval), stat = "identity") +
      scale_fill_gradient2(low = "#5c7aff", high = "#FF4747", 
                           mid = "whitesmoke", midpoint = 0,
                           limits = c(-max_NES-0.5, max_NES+0.5)) + 
      xlab("NES") +
      ylab("Transcription Factors") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold", size = 12),
            axis.text.x = 
              element_text( hjust = 1, size =10),
            axis.text.y = element_text(size =10),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      scale_x_continuous(limits = c(-max_NES-0.5, max_NES+0.5))
    
    ggsave(paste0(path_to_plot, "/Dorothea_", up_down, "reg_TF_activity_barplot.pdf"), plot, width=4, height=3)
  }
}


