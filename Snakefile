from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Van_Allen-data/main/"

rule get_MultiAssayExp:
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        # S3.remote(prefix + "processed/CNA_gene.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/SNV.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData")
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb=4000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v19.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Van_Allen", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v19.annotation.RData")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v19.annotation.RData?raw=true -O {prefix}annotation/Gencode.v19.annotation.RData 
        """

rule format_snv:
    input:
        S3.remote(prefix + "download/SNV.txt.gz"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    output:
        S3.remote(prefix + "processed/SNV.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_SNV.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_expr:
    input:
        S3.remote(prefix + "download/EXPR.txt.gz"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    output:
        S3.remote(prefix + "processed/EXPR.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_EXPR.R \
        {prefix}download \
        {prefix}processed \
        """

# rule format_cna_seg:
#     output:
#         S3.remote(prefix + "processed/CNA_seg.txt")
#     input:
#         S3.remote(prefix + "processed/cased_sequenced.csv"),
#         S3.remote(prefix + "download/CNA_seg.txt.gz")
#     resources:
#         mem_mb=2000
#     shell:
#         """
#         Rscript scripts/Format_CNA_seg.R \
#         {prefix}download \
#         {prefix}processed \
#         """

# rule format_cna_gene:
#     input:
#         S3.remote(prefix + "processed/cased_sequenced.csv"),
#         S3.remote(prefix + "download/gistic/all_thresholded.by_genes.txt.gz")
#     output:
#         S3.remote(prefix + "processed/CNA_gene.csv")
#     resources:
#         mem_mb=2000
#     shell:
#         """
#         Rscript scripts/Format_CNA_gene.R \
#         {prefix}download \
#         {prefix}processed \
#         """

rule format_clin:
    input:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "download/CLIN.txt")
    output:
        S3.remote(prefix + "processed/CLIN.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_CLIN.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_cased_sequenced:
    input:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt.gz"),
        S3.remote(prefix + "download/SNV.txt.gz"),
        # S3.remote(prefix + "download/gistic/all_thresholded.by_genes.txt.gz")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_cased_sequenced.R \
        {prefix}download \
        {prefix}processed \
        """

rule format_downloaded_data:
    input:
        S3.remote(prefix + "download/tables1.mutation_list_all_patients.xlsx"),
        S3.remote(prefix + "download/tables2.clinical_and_genome_characteristics_each_patient.xlsx"),
        S3.remote(prefix + "download/TPM_RSEM_VAScience2015.txt")
    output:
        S3.remote(prefix + "download/CLIN.txt"),
        S3.remote(prefix + "download/EXPR.txt.gz"),
        S3.remote(prefix + "download/SNV.txt.gz")
    shell:
        """
        Rscript scripts/format_downloaded_data.R {prefix}download
        """

rule download_data:
    output:
        S3.remote(prefix + "download/tables1.mutation_list_all_patients.xlsx"),
        S3.remote(prefix + "download/tables2.clinical_and_genome_characteristics_each_patient.xlsx"),
        S3.remote(prefix + "download/TPM_RSEM_VAScience2015.txt"),
        # S3.remote(prefix + "download/gistic/all_thresholded.by_genes.txt.gz"),
    resources:
        mem_mb=2000
    shell:
        """
        wget {data_source}tables1.mutation_list_all_patients.xlsx -O {prefix}download/tables1.mutation_list_all_patients.xlsx
        wget {data_source}tables2.clinical_and_genome_characteristics_each_patient.xlsx -O {prefix}download/tables2.clinical_and_genome_characteristics_each_patient.xlsx
        # wget {data_source}gistic/all_thresholded.by_genes.txt.gz -O {prefix}download/gistic/all_thresholded.by_genes.txt.gz
        wget -O {prefix}download/TPM_RSEM_VAScience2015.txt https://github.com/vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM/raw/master/TPM_RSEM_VAScience2015.txt 
        """ 