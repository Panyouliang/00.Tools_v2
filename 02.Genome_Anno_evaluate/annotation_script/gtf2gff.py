#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This is a script for transform gtf to gff file.
Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-06-13, yyyy-mm-dd
======================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gtf', metavar='annotation file', type=str, required=True, help='Please input the gtf file')
#parser.add_argument('-sp', metavar='Specie name', type=str, required=True, help='Please input the Specie name')
args = parser.parse_args()
#================================================================================


def read(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')


    UTR_list = ['five_prime_UTR','five_prime_utr','three_prime_utr','three_prime_UTR','UTR3','UTR5']
    
    gene_dict,all_dict = {},{}

    for line in f:
        if line[0] != '#':
            line = line.strip().split('\t')
            if line[2] == 'gene':

                gidx = line[8].find('gene_id')
                gid = line[8][gidx:].split(';')[0].split('"')[1]

                if 'gene_biotype' in line[8]:
                    tyidx = line[8].find('gene_biotype')
                else:
                    tyidx = line[8].find('gene_type')

                typ = line[8][tyidx:].split(';')[0].split('"')[1]
                

                if  typ == 'protein_coding':
                    d = line[0:8]
                    d.append('ID='+gid+';')
                    gene_dict[gid] = d
                    all_dict[gid] = {}

            elif line[2] == 'transcript':
                gidx = line[8].find('gene_id')
                gid = line[8][gidx:].split(';')[0].split('"')[1]

                cds,exon,utr = 0,0,0
                line[2] = 'mRNA'
                tidx,ggidx = line[8].find('transcript_id'),line[8].find('gene_id')
                tid,ggid = line[8][tidx:].split(';')[0].split('"')[1],line[8][ggidx:].split(';')[0].split('"')[1]

                
                trans = line[0:8]
                trans.append('ID='+tid+';'+'Parent='+gid+';')
                if ggid in all_dict.keys():
                    all_dict[gid][tid] = [trans]

            elif line[2] == 'exon':
                exon += 1
                idx,ggidx = line[8].find('transcript_id'),line[8].find('gene_id')
                pid,ggid = line[8][idx:].split(';')[0].split('"')[1],line[8][ggidx:].split(';')[0].split('"')[1]
                
                exons = line[0:8]
                exons.append('Parent='+tid+';')
                if ggid in all_dict.keys():
                    all_dict[ggid][pid].append(exons)

            elif line[2] == 'CDS':

                cds += 1
                cidx,ggidx = line[8].find('transcript_id'),line[8].find('gene_id')
                cid,ggid = line[8][cidx:].split(';')[0].split('"')[1],line[8][ggidx:].split(';')[0].split('"')[1]
                
                cdss = line[0:8]
                cdss.append('Parent='+tid+';')
                if ggid in all_dict.keys():
                    all_dict[ggid][cid].append(cdss)


            elif line[2] in UTR_list:

                utr += 1
                uidx,ggidx = line[8].find('transcript_id'),line[8].find('gene_id')
                uid,ggid = line[8][uidx:].split(';')[0].split('"')[1],line[8][ggidx:].split(';')[0].split('"')[1]
                
                utrs = line[0:8]
                utrs.append('Parent='+tid+';')
                if ggid in all_dict.keys():
                    all_dict[ggid][uid].append(utrs)

    return gene_dict,all_dict



def main():
    gtf = args.gtf
    gene_d,all_d = read(gtf)

    for k,v in all_d.items():
        print('\t'.join(x for x in gene_d[k]))
        for mRNA,value in v.items():
            for line in value:
                print('\t'.join(x for x in line))



if __name__ == '__main__':
    main()
