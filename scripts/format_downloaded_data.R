# wget -O ~/Documents/ICBCuration/data_source/Van_Allen/TPM_RSEM_VAScience2015.txt https://github.com/vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM/raw/master/TPM_RSEM_VAScience2015.txt

library(data.table)
library(readxl) 
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

# CLIN.txt
clin <- read_excel(
  file.path(work_dir, 'tables2.clinical_and_genome_characteristics_each_patient.xlsx'), 
  sheet='exome analysis (n=110)'
)

clin_rna <- read_excel(
  file.path(work_dir, 'tables2.clinical_and_genome_characteristics_each_patient.xlsx'), 
  sheet='transcriptome analysis (n=42)'
)
clin_rna$patient_number <- NULL

common_cols <- colnames(clin)[colnames(clin) %in% colnames(clin_rna)]
clin <- clin[, c(common_cols, colnames(clin)[!colnames(clin) %in% common_cols])]
clin_rna <- clin_rna[, c(common_cols, colnames(clin_rna)[!colnames(clin_rna) %in% common_cols])]
clin_cols <- colnames(clin)[!colnames(clin) %in% common_cols]
clin_rna_cols <- colnames(clin_rna)[!colnames(clin_rna) %in% common_cols]

merged_cols <- unique(c(colnames(clin), colnames(clin_rna)))
patients <- unique(c(clin$patient, clin_rna$patient))
clin_merged <- data.frame(matrix(ncol = length(merged_cols), nrow = length(patients)))
colnames(clin_merged) <- merged_cols
rownames(clin_merged) <- patients

for(patient in rownames(clin_merged)){
  if(patient %in% clin$patient){
    clin_merged[rownames(clin_merged) == patient, common_cols] <- unlist(lapply(common_cols, function(col){
      return(clin[clin$patient == patient, col])
    }))
    clin_merged[rownames(clin_merged) == patient, clin_cols] <- unlist(lapply(clin_cols, function(col){
      return(clin[clin$patient == patient, col])
    }))
  }else{
    clin_merged[rownames(clin_merged) == patient, common_cols] <- unlist(lapply(common_cols, function(col){
      return(clin_rna[clin_rna$patient == patient, col])
    }))
  }
  if(patient %in% clin_rna$patient){
    clin_merged[rownames(clin_merged) == patient, clin_rna_cols] <- unlist(lapply(clin_rna_cols, function(col){
      return(clin_rna[clin_rna$patient == patient, col])
    }))
  }
}

write.table(clin_merged, file=file.path(work_dir, 'CLIN.txt'), quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE)


# EXPR.txt.gz
expr <- read.csv(file.path(work_dir, 'TPM_RSEM_VAScience2015.txt'), sep="\t")
colnames(expr) <- str_replace(colnames(expr), 'MEL.IPI_', '')
gz <- gzfile(file.path(work_dir, 'EXPR.txt.gz'), "w")
write.table( expr , file=gz , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
close(gz)

# SNV.txt.gz
snv <- read_excel(
  file.path(work_dir, 'tables1.mutation_list_all_patients.xlsx'), 
  sheet='S1_MEL-All-Muts.csv'
)
gz <- gzfile(file.path(work_dir, 'SNV.txt.gz'), "w")
write.table( snv , file=gz , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
close(gz)

