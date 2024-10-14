library(data.table)

all_genes <- fread('./raw/hg19_biomart_all_genes.txt')
colnames(all_genes) <- c('ens_ID', 'start', 'end', 'chr', 'hgnc_symbol', 'biotype')
all_trans <- fread('./raw/hg19_biomart_all_transcripts.txt')
colnames(all_trans) <- c('ens_ID', 'transcript_ID', 'transcript_biotype', 'refseq_ID')
all_exons <- fread('./raw/hg19_biomart_all_exons.txt')
colnames(all_exons) <- c('transcript_ID', 'exon_start', 'exon_end', 'exon_ID')

system('mkdir processed')

# ensemble to hgnc
fwrite(unique(all_genes[hgnc_symbol != '', .(ens_ID, hgnc_symbol)]),
              'gene_symbols.txt', sep ='\t')
# protein coding only
fwrite(all_genes[biotype == 'protein_coding', ][, .(ens_ID, chr, start, end)],
       './processed/protein_coding_genes.txt', sep = '\t')
# exons in RefSeq transcripts
fwrite(all_exons[transcript_ID %chin% all_trans[refseq_ID != '', transcript_ID], ],
       'exons_in_refseq_transcripts.txt', sep = '\t')
