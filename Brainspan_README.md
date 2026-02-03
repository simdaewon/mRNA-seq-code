# mRNA-seq-code

library(dplyr)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

Cortex <- read.delim("data/Cortex_RPKM.txt")
nrow(Cortex) #52376
ncol(Cortex) #525
colnames(Cortex)

Cortex <- Cortex %>% dplyr::select(!starts_with("AMY"))
ncol(Cortex) #492
colnames(Cortex)

duplicated_genes <- Cortex$gene_symbol[duplicated(Cortex$gene_symbol)]

Cortex_final <- Cortex %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

Cortex_final <- as.data.frame(Cortex_final)

rownames(Cortex_final) <- Cortex_final$gene_symbol
Cortex_final["Y_RNA",]

#write.csv(Cortex_final, "data/Cortex_filter.csv")

Our <- read.delim("data/human_rawcounts.txt")

Our$gene_length <- abs(Our$end - Our$start)

count_data <- Our[, 5:(ncol(Our)-1)]

gene_info <- data.frame(
  Geneid = Our[[1]],
  Length = abs(Our$end - Our$start))


dge <- DGEList(counts = count_data, genes = gene_info)

Our_rpkm <- rpkm(dge, gene.length = dge$genes$Length)
Our_rpkm <- cbind(gene_info[,1, drop=FALSE], Our_rpkm)
rownames(Our_rpkm) <- Our_rpkm$Geneid
Our_rpkm <- Our_rpkm[,-1]

common_genes <- intersect(rownames(Cortex_final), rownames(Our_rpkm))
common_genes <- sort(common_genes)

Cortex_common <- Cortex_final[common_genes, ]
Our_common <- Our_rpkm[common_genes, ]

merged_rpkm <- cbind(Cortex_common, Our_common)

rownames(merged_rpkm) <- common_genes

#write.csv(merged_rpkm,"data/merged_rpkm.csv")

merged_rpkm <- merged_rpkm[,-1]
colnames(merged_rpkm)

sampleinfo <- read.csv("data/sampleinfo.csv", row.names = 1)

sampleinfo <- sampleinfo %>% filter(Region_Class != "AMY")

nrow(sampleinfo)
ncol(merged_rpkm)

Dorsal_Telencephalon <- c("TBR1", "SATB2", "NEUROD6", "MEF2C", "SLC17A7","EMX1","KCNV1","FEZF2","NEUROD2",
                          "PRDM8","CYP26A1","TSHZ3","SLA","NR4A2","KLHL1","ADRA2A","SLN")

Ventral_Telencephalon <- c("PDYN", "SIX3", "RXRG", "DRD2","TAC1","DLX6","GPR6","PBX3","ISL1","FAM40B",
                           "GPR88","DLX5","NKX2-1","DRD1","POU3F4","RARB","SERTAD4","ZNF503","NTN1")

Thalamus <- c("TCF7L2","RGS16","GBX2","LHX9","NTNG1","OTX2","SHOX2","RGS8","CPNE9","ZIC4","CPNE7",
              "KITLG","SOX14")

HIP <- c("NTS","FOXJ1","TNFAIP8L3","MSTN","NRP1")

Cerebellum <- c("CBLN1","ZIC1","PCP2","GABRA6","BARHL1","EN2","ZIC2","MAB21L1","GRM4","TLX3",
                "FAT2","CDH15","CA8","LCAT","TIMP4","BARHL2","LHX1","CALB1")

existing_genes <- intersect(c(Dorsal_Telencephalon, Ventral_Telencephalon, Thalamus, HIP, Cerebellum), rownames(merged_rpkm))

group_info <- sampleinfo$Region_Class[match(colnames(merged_rpkm), sampleinfo$Sample_ID)]

group_mean_df <- aggregate(t(merged_rpkm[existing_genes, ]), 
                           by = list(Group = group_info), 
                           FUN = mean)


rownames(group_mean_df) <- group_mean_df$Group
group_mean_mat <- as.matrix(group_mean_df[, -1])

target_order <- c("CBC", "HIP", "MD", "STR","NCX", "TP_NPCs")

existing_order <- intersect(target_order, rownames(group_mean_mat))
group_mean_mat <- group_mean_mat[existing_order, ]

existing_dorsal <- intersect(Dorsal_Telencephalon, existing_genes)
existing_ventral <- intersect(Ventral_Telencephalon, existing_genes)
existing_Thalamus <- intersect(Thalamus, existing_genes)
existing_HIP <- intersect(HIP, existing_genes)
existing_Cerebellum <- intersect(Cerebellum, existing_genes)


all_marker_genes <- c(existing_dorsal,existing_ventral, existing_Thalamus, existing_HIP, existing_Cerebellum)

group_mean_mat_extended <- group_mean_mat[, all_marker_genes]
colnames(group_mean_mat_extended) <- make.unique(all_marker_genes)

annotation_col_gene <- data.frame(
  Marker = factor(c(rep("Dorsal", length(existing_dorsal)), 
                    rep("Ventral", length(existing_ventral)),
                    rep("Thalamus", length(existing_Thalamus)),
                    rep("HIP", length(existing_HIP)),
                    rep("Cerebellum", length(existing_Cerebellum))),
                  levels = c("Dorsal", "Ventral", "Thalamus", "HIP", "Cerebellum")))

rownames(annotation_col_gene) <- colnames(group_mean_mat_extended)
plot_data <- group_mean_mat_extended

annotation_row_group <- data.frame(
  Group = factor(rownames(plot_data), levels = target_order))

rownames(annotation_row_group) <- rownames(plot_data)

ann_colors = list(
  Marker = c(Dorsal = "#1B9E77", Ventral = "#D95F02", Thalamus = "#7570B3", HIP = "#ba68c8", Cerebellum = "#ff80ab"))


png(filename = "figure/Non_quantile_norm.png",  width = 1200, height = 300, res = 80)
pheatmap(plot_data, 
         annotation_row = annotation_row_group,
         annotation_col = annotation_col_gene,
         annotation_colors = ann_colors,
         scale = "column",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_row = if("TP_NPCs" %in% rownames(plot_data)) which(rownames(plot_data) == "TP_NPCs") - 1 else NULL,
         gaps_col = cumsum(c(length(existing_dorsal), 
                             length(existing_ventral), 
                             length(existing_Thalamus), 
                             length(existing_HIP))),
         angle_col = "45",
         labels_col = all_marker_genes,
         display_numbers = FALSE,
         main = "Average Expression",
         color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
         fontsize_row = 11,
         fontsize_col = 8,
         border_color = "white")
dev.off()

####################################################################################################################
library(preprocessCore)

norm_mat <- normalize.quantiles(as.matrix(merged_rpkm))

rownames(norm_mat) <- rownames(merged_rpkm)
colnames(norm_mat) <- colnames(merged_rpkm)

merged_rpkm_norm <- as.data.frame(norm_mat)

group_mean_df <- aggregate(t(merged_rpkm_norm[existing_genes, ]), 
                           by = list(Group = group_info), 
                           FUN = mean)

rownames(group_mean_df) <- group_mean_df$Group
group_mean_mat <- as.matrix(group_mean_df[, -1])
group_mean_mat <- group_mean_mat[existing_order, ]

plot_data_Quantile <- group_mean_mat[, all_marker_genes]
colnames(plot_data_Quantile) <- rownames(annotation_col_gene)

png(filename = "figure/Quantile_norm_final.png",  width = 2600, height = 700, res = 200)
pheatmap(plot_data_Quantile,
         annotation_col = annotation_col_gene,
         annotation_colors = ann_colors,
         scale = "column",
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         gaps_row = if("TP_NPCs" %in% rownames(plot_data_Quantile)) which(rownames(plot_data_Quantile) == "TP_NPCs") - 1 else NULL,
         gaps_col = cumsum(c(length(existing_dorsal), 
                             length(existing_ventral), 
                             length(existing_Thalamus), 
                             length(existing_HIP))),
         labels_col = all_marker_genes,
         angle_col = "45",
         main = "Comparison of marker gene expression",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 11, fontsize_col = 8,
         border_color = "white",
         annotation_names_col = FALSE,
         cellwidth = 10, cellheight = 20,)
dev.off()
#####################################################################################################################################
#####################################################################################################################################
target_ages <- c("8 pcw", "9 pcw", "12 pcw", "13 pcw", "16 pcw", "17 pcw", 
                 "19 pcw", "21 pcw", "24 pcw", "25 pcw", "26 pcw", "35 pcw", "37 pcw")

target_order <- c("CBC", "HIP", "MD", "STR","NCX", "TP_NPCs")

comparison_info <- sampleinfo %>%
  filter(age %in% target_ages | Region_Class == "TP_NPCs")

comparison_info_final <- comparison_info %>%
  filter(Region_Class %in% target_order)

target_samples <- comparison_info_final$Sample_ID
valid_target_samples <- intersect(target_samples, colnames(merged_rpkm_norm))
subset_rpkm <- merged_rpkm_norm[existing_genes, valid_target_samples]

region_group_info <- comparison_info_final$Region_Class[match(colnames(subset_rpkm), comparison_info_final$Sample_ID)]

region_mean_df <- aggregate(t(subset_rpkm), 
                            by = list(Group = region_group_info), 
                            FUN = mean)

rownames(region_mean_df) <- region_mean_df$Group
region_mean_mat <- as.matrix(region_mean_df[, -1])


existing_order <- intersect(target_order, rownames(region_mean_mat))
region_mean_mat <- region_mean_mat[existing_order, ]

plot_data_final <- region_mean_mat[, all_marker_genes]

pheatmap(plot_data_final,
         annotation_col = annotation_col_gene,
         annotation_colors = ann_colors,
         scale = "column",
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         gaps_row = if("TP_NPCs" %in% rownames(plot_data_final)) which(rownames(plot_data_final) == "TP_NPCs") - 1 else NULL,
         gaps_col = cumsum(c(length(existing_dorsal), 
                             length(existing_ventral), 
                             length(existing_Thalamus), 
                             length(existing_HIP))),
         labels_col = all_marker_genes,
         angle_col = "45",
         main = "Marker Gene Expression by Target Regions (Filtered Ages)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 11, 
         fontsize_col = 8,
         border_color = "white",
         cellwidth = 12, 
         cellheight = 20)

##################################################################################################################################

