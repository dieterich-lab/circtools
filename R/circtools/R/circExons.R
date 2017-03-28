if ( FALSE ) {
library(dplyr)
library(biomaRt)
library(Biostrings)


get_exons <- function(t_chr, t_start, t_end, mart){
    s <- getSequence(chr=t_chr, start=t_start, end=t_end, mart=human,
        type=c(
            "ensembl_gene_id",
            "ensembl_exon_id",
            "ensembl_transcript_id", 
            "chromosome_name",
            "exon_chrom_start", 
            "exon_chrom_end",
            "strand", 
            "external_gene_name"), 
        seqType="gene_exon") 
    left_exon<- s %>% filter(exon_chrom_start==t_start) 
    right_exon<- s %>% filter(exon_chrom_end==t_end) 
    list(left_exon=left_exon,right_exon=right_exon)
}

}