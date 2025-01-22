#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse
import numpy as np

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
============================================================
This script is used to rename the gene ID of the genome annotation file.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-07-07, yyyy-mm-dd
============================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gff', metavar='annotation file', type=str, required=False, help='Please input the annotation file(None isoform)')
parser.add_argument('-name', metavar='Specie name', type=str, required=True, help='Please input the abbreviation of the species that you want to Rename')
args = parser.parse_args()
#=================================================================================

def read_gff(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    chr_dict,sort_gene = {},{}
    utr5,utr3 = ['five_prime_UTR','five_prime_utr','UTR5'],['three_prime_UTR','three_prime_utr','UTR3']
    old_id = {}
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                line[1] = 'ReName'
                if line[2] == 'gene':
                    chrid = line[0]
                    geneID = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    gid = geneID
                    if chrid in sort_gene.keys():
                        sort_gene[chrid].append([geneID,int(line[3]),int(line[4]),line[6]])
                    else:
                        sort_gene[chrid] = [[geneID,int(line[3]),int(line[4]),line[6]]]

                if line[2] == 'mRNA':
                    idx = line[8].find('Parent') if 'Parent' in line[8] else line[8].find('gene_id')
                    chrid,tid = line[0],line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    geneID = line[8][idx:].split('=')[1].split(';')[0]
                    old_id[geneID] = tid
                    if chrid not in sort_gene:
                        sort_gene[chrid] = [[geneID,int(line[3]),int(line[4]),line[6]]]
                    else:
                        sort_gene[chrid].append([geneID,int(line[3]),int(line[4]),line[6]])

                    if chrid in chr_dict.keys():
                        if geneID in chr_dict[chrid].keys():
                            if tid in chr_dict[chrid][geneID].keys():
                                chr_dict[chrid][geneID][tid].append(line)
                            else:
                                chr_dict[chrid][geneID][tid] = [line]
                        else:
                            chr_dict[chrid][geneID] = {tid:[line]}
                    else:
                        chr_dict[chrid] = {geneID:{tid:[line]}}


                if line[2] == 'CDS':
                    pchrid = line[0]
                    ptid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    if pchrid in chr_dict.keys() and ptid in chr_dict[pchrid][geneID].keys():
                        chr_dict[pchrid][geneID][ptid].append(line)
                if line[2] == 'exon':
                    pchrid = line[0]
                    ptid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    if pchrid in chr_dict.keys() and ptid in chr_dict[pchrid][geneID].keys():
                        chr_dict[pchrid][geneID][ptid].append(line)

                if line[2] in utr5:
                    pchrid = line[0]
                    ptid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    if pchrid in chr_dict.keys() and ptid in chr_dict[pchrid][geneID].keys():
                        chr_dict[pchrid][geneID][ptid].append(line)

                if line[2] in utr3:
                    pchrid = line[0]
                    ptid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    if pchrid in chr_dict.keys() and ptid in chr_dict[pchrid][geneID].keys():
                        chr_dict[pchrid][geneID][ptid].append(line)

    sort_gene_dict = {}
    for k,v in sort_gene.items():
        sort = sorted(v,key=lambda x: int(x[1]))
        sort_gene_dict[k] = sort
        #print(sort)


    for chrs,gvalue in chr_dict.items():
        for gid,tvalue in gvalue.items():
            for tid,line in tvalue.items():
                tagm = line[0]
                if tagm[6] == '-':
                    line = sorted(line[1:], key=lambda x:int(x[3]),reverse=True)
                else:
                    line = sorted(line[1:], key=lambda x:int(x[3]),reverse=False)

                line.insert(0,tagm)


    return chr_dict,sort_gene_dict,old_id

def out_p(chr_dict,gene_sort,old_id,name):
    print('##gff-version3')
    gnum = 0
    out = open('result_rename.lis','w')
    out.write('#old-gene-id\tnew-gene-id\n')
    for chr,gid_lis in gene_sort.items():
        #print(chr,gid_lis)
        for gid in gid_lis:
            gnum += 1
            iso = 0
            geneID = str(gnum).zfill(6)
            out.write(gid[0]+'\t'+name+'G'+geneID+'\n')
            print(chr,'ReName','gene',gid[1],gid[2],'.',gid[3],'.','ID='+name+'G'+geneID+';',sep='\t')
            for tid,tvalue in chr_dict[chr][gid[0]].items():
                #print(tid,tvalue)
                iso += 1
                utr5,utr3,cds,exon = 0,0,0,0
                for line in tvalue:
                    #print(line)
                    if line[2] == 'mRNA':
                        print('\t'.join(x for x in line[0:8]),'ID='+name+'G'+geneID+'.'+str(iso)+';Parent='+name+'G'+geneID+';',sep='\t')
                    if line[2] == 'UTR5' or line[2] == 'five_prime_UTR':
                        utr5 += 1
                        print('\t'.join(x for x in line[0:8]),'Parent='+name+'G'+geneID+'.'+str(iso)+';',sep='\t')
                    if line[2] == 'exon':
                        exon += 1
                        print('\t'.join(x for x in line[0:8]),'Parent='+name+'G'+geneID+'.'+str(iso)+';',sep='\t')
                    if line[2] == 'CDS':
                        cds += 1
                        print('\t'.join(x for x in line[0:8]),'Parent='+name+'G'+geneID+'.'+str(iso)+';',sep='\t')
                    if line[2] == 'UTR3' or line[2] == 'three_prime_UTR':
                        utr3 += 1
                        print('\t'.join(x for x in line[0:8]),'Parent='+name+'G'+geneID+'.'+str(iso)+';',sep='\t')

def main():
    gff,name = args.gff,args.name

    chr_dict,gene_sort,old_id=read_gff(gff)

    out_p(chr_dict,gene_sort,old_id,name)

if __name__ == '__main__':
    main()
