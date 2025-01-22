#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip


def readgtf(F):
    gened = {}
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')
    for line in f:
        if line[0] == '#':
            pass
        else:
            line = line.strip().split('\t')
            if line[2] == 'gene':
                pass
            else:
                geneID = line[8][line[8].find('gene_id'):].split('"')[1]
                rnaID  = line[8][line[8].find('transcript_id'):].split('"')[1]
                if line[2] == 'transcript':
                    mline = line[0:8].copy()
                    mline.append(rnaID)
                    mline[2] = 'mRNA'
                    gened[geneID] = [mline]
                if line[2] == 'CDS' or line[2] == 'exon':
                    cds = line[0:8].copy()
                    cds.append(rnaID)
                    gened[geneID].append(cds)
    return gened


gened = readgtf(sys.argv[1])

for k,v in gened.items():
    for line in v:
        if line[2] == 'mRNA':
            print('\t'.join(x for x in line[0:8]),'ID='+line[8]+';Parent='+k+';',sep='\t')
        else:
            print('\t'.join(x for x in line[0:8]),'Parent='+line[8]+';',sep='\t')


