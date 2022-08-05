library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

clin = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t" )

patient = sort( clin$patient )
case = as.data.frame( cbind( patient , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) , rep( 0 , length(patient) ) ) )
colnames(case) = c( "patient" , "snv" , "cna" , "expr" )
rownames(case) = patient
case$snv = as.numeric( as.character( case$snv ) )
case$cna = as.numeric( as.character( case$cna ) )
case$expr = as.numeric( as.character( case$expr ) )

expr <- as.matrix( fread( file.path(input_dir, "EXPR.txt.gz") , stringsAsFactors=FALSE  , sep="\t" , dec=',') )
expr_patients <- colnames(expr)[-1]
snv = as.data.frame( fread( file.path(input_dir, "SNV.txt.gz") , stringsAsFactors=FALSE , sep="\t" ))
snv_patients <- unique(snv$patient)
# cna = as.data.frame( fread( file.path(input_dir, "gistic/all_thresholded.by_genes.txt.gz") , stringsAsFactors=FALSE , sep="\t" , dec="," ))
# cna_patients <- colnames(cna)[-c(1:3)]

for( i in 1:nrow(case)){
	if( rownames(case)[i] %in% expr_patients ){
	  case$expr[i] = 1
	}
	if( rownames(case)[i] %in% snv_patients ){
		case$snv[i] = 1
	}
#   if( rownames(case)[i] %in% cna_patients ){
#     case$cna[i] = 1
#   }
}

write.table( case , file=file.path(output_dir, "cased_sequenced.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
