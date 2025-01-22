export PATH="/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin:$PATH"
source /ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/etc/profile.d/conda.sh

gff=
genome=

# extract fasta from gff
gffread $gff -g $genome -w transcripts.fa

# transform gff to gff3
gffread $gff -T -o gff2gtf.gtf
gtf_to_alignment_gff3.pl gff2gtf.gtf >transcripts.gff3

# extract the long open reading frames
TransDecoder.LongOrfs -t transcripts.fa -m 200

# predict the likely coding regions
TransDecoder.Predict -t transcripts.fa --single_best_only

# generate a genome-based coding region annotation file
cdna_alignment_orf_to_genome_orf.pl transcript.fa.transdecoder.gff3 transcripts.gff3 transcripts.fa > TransDecoder.genome.gff3
