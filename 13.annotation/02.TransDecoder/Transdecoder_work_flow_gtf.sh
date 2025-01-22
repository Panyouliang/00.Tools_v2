export PATH="/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin:$PATH"
source /ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/etc/profile.d/conda.sh

gtf=
genome=

# extract fasta from gtf
gtf_genome_to_cdna_fasta.pl $gtf $genome > transcripts.fa

# transform gtf to gff3
gtf_to_alignment_gff3.pl $gtf >transcripts.gff3

# extract the long open reading frames
TransDecoder.LongOrfs -t transcripts.fa -m 200

# predict the likely coding regions
TransDecoder.Predict -t transcripts.fa --single_best_only

# generate a genome-based coding region annotation file
cdna_alignment_orf_to_genome_orf.pl transcript.fa.transdecoder.gff3 transcripts.gff3 transcripts.fa > TransDecoder.genome.gff3
