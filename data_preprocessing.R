# import libraries
library(readxl)
library(dplyr)

# define path
# path_to_project_file <- "/Users/Danae/Documents/MIT Internship/Jan-Willem Human Study"
path_to_project_file <- "~/Dropbox (MIT)/metabolism_Joslin/Jan_Willem/Jan-Willem Human Study_RNAseq/"
path_to_intial_data <- paste0(path_to_project_file, "/Data")
# set working directory
setwd(path_to_intial_data)

# read bulk RNA-seq data
bulk_data <- as.data.frame(read_excel("Bulk_gene_count_matrix.xlsx", sheet = 1, col_names = TRUE))
bulk_data <- bulk_data %>% arrange(`HGNC symbol`, desc(MPG))
bulk_data <- bulk_data[!duplicated(bulk_data$`HGNC symbol`) & !is.na(bulk_data$`HGNC symbol`), ]

# tpm dataframe
tpm <- bulk_data[bulk_data$`Gene type` == "protein_coding", grepl("l2tpm", colnames(bulk_data))]
rownames(tpm) <- bulk_data[bulk_data$`Gene type` == "protein_coding", "HGNC symbol"]
colnames(tpm) <- sapply(strsplit(colnames(tpm), "\\."), `[`, 1)

# counts dataframe
cts <- bulk_data[bulk_data$`Gene type` == "protein_coding", grepl("intCt", colnames(bulk_data))]
rownames(cts) <- bulk_data[bulk_data$`Gene type` == "protein_coding", "HGNC symbol"]
colnames(cts) <- sapply(strsplit(colnames(cts), "\\."), `[`, 1)

# metadata about each sample
coldata <- as.data.frame(colnames(tpm))
colnames(coldata) <- "sample"
bulk_data_grouping <- as.data.frame(read_excel("Human adipose tissue samples Joslin - for QA.xlsx", sheet = 4, col_names = TRUE, range = "A1:C61"))
colnames(bulk_data_grouping)[3] <- "subject"
coldata <- cbind(coldata, bulk_data_grouping[, c("Group", "subject")])
rownames(coldata) <- coldata$sample
coldata <- coldata[, -1]
coldata$Group <- as.factor(coldata$Group)
coldata$time <- substr(rownames(coldata), 1, 1)
coldata <- coldata[order(coldata$Group),]
coldata[coldata$subject %in% c(5,6,15,17,29,30), "Subject_recode"] <- seq(1, 6)
coldata[coldata$subject %in% c(2,3,4,7,9,10,18,21,24), "Subject_recode"] <- seq(1, 9)
coldata[coldata$subject %in% c(11,13,19,20,22,23,25,26), "Subject_recode"] <- seq(1, 8)
coldata[coldata$subject %in% c(1,8,12,14,16,27,28), "Subject_recode"] <- seq(1, 7)
covariate <- as.data.frame(read_excel("Clinical Parameters 6.8.2023_UPDATED.xlsx", sheet = 1, col_names = TRUE, range = "B2:DE32"))
coldata <- merge(coldata, covariate[, c(1, 3:42, 44:ncol(covariate))], by.x = "subject", by.y = "Subject Number")
colnames(coldata) <- gsub("Δ ", "delta_", colnames(coldata))
colnames(coldata) <- gsub("Δ", "delta_", colnames(coldata))
colnames(coldata) <- gsub(" ", "_", colnames(coldata))
colnames(coldata) <- gsub("%", "percent", colnames(coldata))
colnames(coldata) <- gsub("\\)", "", gsub("\\(", "", colnames(coldata)))
colnames(coldata) <- gsub("\\/", "", colnames(coldata))
colnames(coldata) <- gsub("-_", "", colnames(coldata))
rownames(coldata) <- paste0(coldata$time, coldata$subject)
coldata[coldata=="N/A"] <- NA
coldata[, 8:ncol(coldata)] <- lapply(coldata[, 8:ncol(coldata)], as.numeric)
coldata$subject <- as.character(coldata$subject)

# add days between last exercise session
ex_biopsy_visit <- as.data.frame(read_excel("Exercise biopsy visit.xlsx", sheet = 1, col_names = TRUE))
colnames(ex_biopsy_visit)[2] <- "Days_between_last_ex_session"
coldata <- merge(coldata, ex_biopsy_visit, by.x = 0, by.y = "Subject", all = TRUE)
# coldata$Days_between_last_ex_session <- as.factor(coldata$Days_between_last_ex_session)
coldata <- coldata[match(colnames(tpm), coldata$Row.names),]
rownames(coldata) <- coldata$Row.names
coldata <- coldata[, -1]

tpm <- tpm[, match(rownames(coldata), colnames(tpm))]
cts <- cts[, match(rownames(coldata), colnames(cts))]
coldata$large_group <- ifelse(coldata$Group == 1 | coldata$Group == 4, "lean", "obese")
# coldata$sex_large_group <- paste0(coldata$Sex, ' ', coldata$large_group)
coldata$sex_large_group <- ifelse(coldata$Sex == "m", paste0(coldata$large_group, ' male'), 
                                  paste0(coldata$large_group, ' female'))
# exclude subject 3
tpm <- tpm[, !colnames(tpm) %in% c("A3", "B3")]
cts <- cts[, !colnames(cts) %in% c("A3", "B3")]
coldata <- coldata[!rownames(coldata) %in% c("A3", "B3"),]

coldata$group_name <- NA
coldata[coldata$Group == "1", ]$group_name <- "MIT lean"
coldata[coldata$Group == "2", ]$group_name <- "non-diabetic obese"
coldata[coldata$Group == "3", ]$group_name <- "T2D obese"
coldata[coldata$Group == "4", ]$group_name <- "HIIT lean"

rm(bulk_data, bulk_data_grouping, covariate, ex_biopsy_visit)
