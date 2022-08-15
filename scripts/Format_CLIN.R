library(tibble)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
annot_dir <- args[3]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )
selected_cols <- c( "patient","age_start","RECIST","overall_survival","progression_free","primary","histology","stage","gender","dead","progression"  )
clin = cbind( clin_original[ , selected_cols ] , "CTLA4" , NA , NA , NA , NA )
colnames(clin) = c( "patient" , "age" , "recist" , "t.os" ,"t.pfs" , "primary" ,"histo" , "stage" , "sex" , "os" , "pfs"  , "drug_type" , "dna" , "rna" , "response.other.info" , "response" )

clin$t.os = clin$t.os/30.5
clin$t.pfs = clin$t.pfs/30.5
clin$sex = ifelse(clin$sex %in% "female" , "F" , "M")
clin$stage = ifelse(clin$stage %in% "Stage 4" , "IV" , "III")
clin$recist[ clin$recist %in% "X" ] = NA
clin$response = Get_Response( data=clin )


case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )
clin$rna[ clin$patient %in% case[ case$expr %in% 1 , ]$patient ] = "tpm"
clin$dna[ clin$patient %in% case[ case$cna %in% 1 , ]$patient ] = "wes"

clin$primary = "Melanoma"

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]

clin <- format_clin_data(clin_original, 'patient', selected_cols, clin)
clin <- clin[, colnames(clin)[colnames(clin) != 'X']]

# Tissue and drug annotation
annotation_tissue <- read.csv(file=file.path(annot_dir, 'curation_tissue.csv'))
clin <- annotate_tissue(clin=clin, study='Van_Allen', annotation_tissue=annotation_tissue, check_histo=FALSE)

clin <- add_column(clin, unique_drugid='', .after='unique_tissueid')

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
