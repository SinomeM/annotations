library(data.table)

all_genes <- fread('./raw/hg38_biomart_all_genes.txt')
colnames(all_genes) <- c('ens_ID', 'start', 'end', 'chr', 'biotype', 'hgnc_symbol')
all_trans <- fread('./raw/hg38_biomart_all_transcripts.txt')
colnames(all_trans) <- c('ens_ID', 'transcript_ID', 'transcript_biotype',
                         'refseq_ID')
all_exons <- fread('./raw/hg38_biomart_all_exons.txt.gz')
colnames(all_exons) <- c('transcript_ID', 'exon_ID', 'exon_start', 'exon_end')

gnomad <- fread('./raw/gnomad.v4.1.constraint_metrics.tsv')

system('mkdir processed')

# ensemble to hgnc
fwrite(unique(all_genes[hgnc_symbol != '', .(ens_ID, hgnc_symbol)]),
              './processed/gene_symbols.txt', sep ='\t')
# protein coding only
fwrite(all_genes[biotype == 'protein_coding', ][, .(ens_ID, chr, start, end)],
       './processed/protein_coding_genes.txt', sep = '\t')
# exons in RefSeq transcripts
dt <- merge(all_exons, all_trans[, .(ens_ID, transcript_ID)],
            by = 'transcript_ID')
dt <- merge(dt, all_genes[, .(ens_ID, chr)], by = 'ens_ID')
fwrite(dt[transcript_ID %chin% all_trans[refseq_ID != '', transcript_ID], ][
            , .(exon_ID, ens_ID, chr, exon_start, exon_end)],
       './processed/exons_in_refseq_transcripts.txt', sep = '\t')

# gene constrain
dt <- unique(gnomad[, .(gene, gene_id, transcript, lof.oe_ci.upper)])
setnames(dt, 'lof.oe_ci.upper', 'LOEUF')
fwrite(dt, './processed/gene_contrain_gnomad.txt', sep = '\t')




# New version

gnomad <- fread('./raw/gnomad.v4.1.constraint_metrics.tsv')[, .(gene_id, transcript, lof.oe_ci.upper)]
colnames(gnomad) <- c('ens_ID', 'trans_ID', 'LOEUF')
all_genes <- fread('./raw/hg38_biomart_all_genes.txt')
colnames(all_genes) <- c('ens_ID', 'start', 'end', 'chr', 'biotype', 'hgnc_symbol')
fread('./raw/genes_coords.txt')
gt_biotypes <- fread('./raw/gene_trans_biotypes.txt')
colnames(gt_biotypes) <- c('ens_ID', 'trans_ID', 'trans_type', 'biotype')
exons_coord <- fread('./raw/exons_coord.txt')
colnames(exons_coord) <- c('trans_ID', 'exon_start', 'exon_end', 'exon_ID')
gte_ids <- fread('./raw/gene_trans_exon_IDs.txt')


dir.create('./processed_new')
fwrite(gnomad[!is.na(LOEUF), ], './processed_new/genes_trans_loeuf.txt', sep = '\t')
fwrite(all_genes[, .(ens_ID, hgnc_symbol, chr, start, end)],
       './processed_new/all_genes.txt', sep = '\t')

dt <- gt_biotypes[trans_type == biotype &
                  biotype %in% c('lncRNA', 'miRNS', 'snRNA', 'misc_RNA', 'snoRNA',
                                 'protein_coding'), .(ens_ID, trans_ID, biotype)]
fwrite(dt, './processed_new/genes_trans_biotype_select.txt', sep = '\t')

exons <- merge(exons_coord, dt[, .(ens_ID, trans_ID)])
fwrite(exons[, .(trans_ID, exon_ID, exon_start, exon_end)], 'processed_new/exons_select.txt', sep = '\t')
