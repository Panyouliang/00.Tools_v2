#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip



seq_d = {}
chrid,seq = '',[]
with gzip.open(sys.argv[1],'rt') as f:
    for line in f:
        if line.strip():
            if line[0] == '>':
                seq_d[chrid] = ''.join(x for x in seq)
                chrid = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())

    seq_d[chrid] = ''.join(x for x in seq)


for k,v in seq_d.items():
    if k == 'chrMT' or k == 'chrM':
        print(k,'ReName','gene',1,len(v),'.','+','.','gene_id "MitoG01"; gene_name "MitoG01"; gene_biotype "protein-coding";',sep='\t')
        print(k,'ReName','transcript',1,len(v),'.','+','.','gene_id "MitoG01"; transcript_id "MitoT01"; gene_name "MitoG01"; gene_biotype "protein-coding";',sep='\t')
        print(k,'ReName','exon',1,len(v),'.','+','.','gene_id "MitoG01"; transcript_id "MitoT01"; gene_name "MitoG01"; exon_id "MitoT01.exon"; transcript_name "MitoT01"; gene_biotype "protein-coding";',sep='\t')
        print(k,'ReName','CDS',1,len(v),'.','+','.','gene_id "MitoG01"; transcript_id "MitoT01"; gene_name "MitoG01"; exon_id "MitoT01.exon"; transcript_name "MitoT01"; gene_biotype "protein-coding";',sep='\t')


        print(k,'ReName','gene',1,len(v),'.','-','.','gene_id "MitoG02"; gene_name "MitoG02"; gene_biotype "protein-coding";',sep='\t')
        print(k,'ReName','transcript',1,len(v),'.','-','.','gene_id "MitoG02"; transcript_id "MitoT02"; gene_name "MitoG02"; gene_biotype "protein-coding";',sep='\t')
        print(k,'ReName','exon',1,len(v),'.','-','.','gene_id "MitoG02"; transcript_id "MitoT02"; gene_name "MitoG02"; exon_id "MitoT02.exon"; transcript_name "MitoT02"; gene_biotype "protein-coding";',sep='\t')
        print(k,'ReName','CDS',1,len(v),'.','-','.','gene_id "MitoG02"; transcript_id "MitoT02"; gene_name "MitoG02"; exon_id "MitoT02.exon"; transcript_name "MitoT02"; gene_biotype "protein-coding";',sep='\t')
