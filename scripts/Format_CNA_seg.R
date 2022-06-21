library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

case = read.csv( file.path(output_dir, "cased_sequenced.csv"), stringsAsFactors=FALSE , sep=";" )

cna = as.data.frame( fread( file.path(input_dir, "CNA_seg.txt.gz") , stringsAsFactors=FALSE , sep="\t" , dec="," ))
cna$chrom = sapply( cna$chrom , function(x){ paste( "chr" , x , sep="" ) } )

cna = cna[ cna$sample %in% case[ case$cna %in% 1 , ]$patient , ]

write.table( cna , file= file.path(output_dir, "CNA_seg.txt") , quote=FALSE , sep="\t" , col.names=TRUE , row.names=FALSE )
