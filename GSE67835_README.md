#A survey of human brain transcriptome diversity at the single cell level

library(GEOquery)
library(data.table)
library(dplyr)
library(purrr)

getGEOSuppFiles("GSE67835")

dir.create("GSE67835_RAW_FILES", showWarnings = FALSE)

untar("GSE67835/GSE67835_RAW.tar", exdir = "GSE67835_RAW_FILES")

file_list <- list.files("GSE67835_RAW_FILES", pattern = "\\.csv.gz$", full.names = TRUE)
length(file_list)


read_sample <- function(f) {
  dt <- fread(f)
  sample_id <- regmatches(basename(f), regexpr("GSM[0-9]+", basename(f)))
  colnames(dt) <- c("Gene", sample_id)
  return(unique(dt, by = "Gene"))
}

message("파일 읽는 중... 총 ", length(file_list), "개")
all_data_list <- lapply(file_list, read_sample)

final_counts_dt <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), all_data_list)
final_counts_dt[is.na(final_counts_dt)] <- 0
final_counts <- as.data.frame(final_counts_dt)
rownames(final_counts) <- final_counts$Gene
final_counts$Gene <- NULL

dim(final_counts)
head(final_counts[, 1:5])


gse <- getGEO("GSE67835", GSEMatrix = TRUE, getGPL = FALSE)


meta1 <- pData(gse[[1]])
meta2 <- pData(gse[[2]])


final_meta_all <- bind_rows(meta1, meta2)


common_samples <- intersect(colnames(final_counts), rownames(final_meta_all))
message("전체 샘플 수: ", length(common_samples), "개")


final_counts <- final_counts[, common_samples]
final_meta <- final_meta_all[common_samples, ]


if(ncol(final_counts) == 466) {
  message("모든 샘플(466개)이 통합.")
} else {
  message("샘플 수가 부족. 현재 수: ", ncol(final_counts))}


all(colnames(final_counts) == rownames(final_meta)) #TRUE

final_counts
final_meta
write.csv(final_counts, file = "GSE67835/GSE67835_merged_counts.csv", row.names = TRUE)
write.csv(final_meta, file = "GSE67835/GSE67835_merged_metadata.csv", row.names = TRUE)

############################################################################################################################
#Single-cell alternative splicing analysis with Expedition reveals splicing dynamics during neuron differentiation

#GSE85908
getGEOSuppFiles("GSE85908")

library(data.table)

counts_85908 <- fread("GSE85908/GSE85908_expression.csv.gz")
counts_85908_t <- t(counts_85908[, -1])
print(counts_85908_t[1:10, 1:10])

sample_names <- counts_85908$V1
colnames(counts_85908_t) <- sample_names

print(counts_85908_t[1:5, 1:5])   #Transcripts per million (TPM) 값으로 돼있음

meta_85908 <- fread("GSE85908/GSE85908_metadata.csv.gz")
colnames(meta_85908)

final_meta_85908 <- meta_85908[, .(sample_id = V1, phenotype)]

all(colnames(counts_85908_t) == final_meta_85908$sample_id) 

final_meta_85908 <- final_meta_85908[match(colnames(counts_85908_t), sample_id)]
table(final_meta_85908$phenotype)

