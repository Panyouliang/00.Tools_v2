#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,gzip,argparse
import numpy as np

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
============================================================
This is a script for evaluating the quality of genome annotation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-06-13, yyyy-mm-dd
============================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-gff', metavar='annotation file', type=str, required=False, help='Please input the annotation file(None isoform)')
parser.add_argument('-ref',  metavar='genome file', type=str, default="None", required=False, help='Please input the genome file (.fasta)')
parser.add_argument('-sp', metavar='Specie name', type=str, required=True, help='Please input the Specie name')
args = parser.parse_args()
#=================================================================================



def get_gff(gff):

    cds_dict, orf_lis, codingLengthlis, tran_dict, cds_pos, gene_num = {},[],[],{},{},0

    feature_d = {}

    if gff[-3:] == '.gz':
        f = gzip.open(gff,'rt')
    else:
        f = open(gff,'r')
    codingL = 0
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    pass
                elif line[2] == 'mRNA' or line[2] == 'transcript':
                    gene_num += 1
                    flag = 0
                    g_id = line[8].split(';')[0].split('=')[1]
                    tran_id = g_id +' '+ line[6]
                    cds_dict[g_id] = []
                    orf_lis.append(int(line[4])-int(line[3])+1)
                    tran_list = []
                    codingLengthlis.append(codingL)
                    codingL = 0
                    feature_d[g_id] = []
                    cds_pos[g_id] = []

                else:
                    #if line[2] == 'exon' or line[2] == 'CDS':
                    feature_d[g_id].append([int(line[3]),int(line[4])])

                    if line[2] == 'CDS':
                        cds_dict[g_id].append(abs(int(line[4])-int(line[3]))+1)
                        cds_pos[g_id] += [int(line[3]),int(line[4])]
                        tran_pos = [line[0],line[3],line[4]]
                        tran_list.append(tran_pos)
                        tran_dict[tran_id] = tran_list
                        codingL += int(line[4])-int(line[3])+1

    codingLengthlis.append(codingL)
    gff_dict = {}
    for k,v in tran_dict.items():
        str_list = sorted(v, key = (lambda x:int(x[1])))
        gff_dict[k] = str_list

    feature_adj = {}
    for k,v in feature_d.items():
        vsort = sorted(v, key=lambda x:int(x[1]))
        feature_adj[k] = merge_adjacent_coordinates(vsort)

    return cds_dict, orf_lis, gff_dict, gene_num, codingLengthlis, feature_adj, cds_pos


def merge_adjacent_coordinates(coordinates):
    merged_coordinates = []
    if not coordinates:
        return merged_coordinates

    current_start, current_end = coordinates[0]
    for coord in coordinates[1:]:
        if coord[0] <= current_end + 1:
            current_end = max(current_end, coord[1])
        else:
            merged_coordinates.append([current_start, current_end])
            current_start, current_end = coord

    merged_coordinates.append([current_start, current_end])

    return merged_coordinates



def get_genome(F1):
    if F1[-3:] == '.gz':
        f = gzip.open(F1,'rt')
    else:
        f = open(F1,'r')
    genome_dict,seqlist,chr = {},[],''

    for line in f:
        if line.strip():
            line = line.strip()
            if line[0] == '>':
                genome_dict[chr] = ''.join(x for x in seqlist)
                chr = line.split()[0][1:]
                seqlist = []
            else:
                seqlist.append(line)

    genome_dict[chr] = ''.join(x for x in seqlist)
    return genome_dict


def get_mRNA(gff_dict,genome_dict):
    sequence = {}
    for k,v in gff_dict.items():
        if k:
            gid = k
            sequence[gid] = ''
            fragment = ''
            for line in v:
                fragment = genome_dict[line[0]][int(line[1])-1:int(line[2])]
                sequence[gid] += fragment
    return sequence

def rev(seq):
    base = {
            'A':'T','T':'A','G':'C','C':'G','N':'N','n':'n','a':'t','t':'a','c':'g','g':'c'
            }
    seq_list = list(reversed(seq))
    seq_rev = [base[k] for k in seq_list]
    seq_list = ''.join(seq_rev)
    return seq_list

def extract_mRNA(fragment):
    mRNA = {}
    for k,v in fragment.items():
        if k.endswith('-'):
            string = rev(v)
            mRNA[k] = string.upper()
        elif k.endswith('+'):
            mRNA[k] = v.upper()

    return mRNA


code = {
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'TGC':'C','TGT':'C',
    'GAC':'D','GAT':'D',
    'GAA':'E','GAG':'E',
    'TTC':'F','TTT':'F',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'CAC':'H','CAT':'H',
    'ATA':'I','ATC':'I','ATT':'I',
    'AAA':'K','AAG':'K',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L','TTA':'L','TTG':'L',
    'ATG':'M',
    'AAC':'N','AAT':'N',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R','AGA':'R','AGG':'R',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S','AGC':'S','AGT':'S',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'TGG':'W',
    'TAC':'Y','TAT':'Y',
    'TAA':'*','TAG':'*','TGA':'*',
    'NNN':'*'
    }

stop_codon = ['TAA','TAG','TGA']


def main():

    gff,ref,name = args.gff,args.ref,args.sp

    if gff:

        cds_dict, orf_lis, gff_dict, gene_num, codingLengthlis, feature_adj, cds_pos = get_gff(gff)

        # total genes
        total_gene = gene_num

        # single exon genes
        #single = open('single-exons-gene.lis','w')
        single_exon_gene = 0
        for k,v in feature_adj.items():
            if len(v) == 1:
                single_exon_gene += 1
                #single.write(k+'\n')
        #single.close()

        # exons per gene
        # exon length (mean)
        exon_length,exons = 0,0
        for v in feature_adj.values():
            for cut in v:
                exons += 1
                exon_length += int(cut[1])-int(cut[0])+1

        exon_length_mean = exon_length / exons
        exon_per_gene = exons / gene_num


        # coding length (mean)
        cds_list = sum(cds_dict.values(), [])
        mRNA_length_mean = sum(cds_list) / gene_num
        mRNA_length_median = np.median(np.array(codingLengthlis))

        # locus length (mean)
        gene_length_mean = sum(orf_lis) / len(orf_lis)
        locus_length_median = np.median(np.array(orf_lis))

        orfspanlis = []
        for k,v in cds_pos.items():
            maxv,minv = max(v),min(v)
            orfspanlis.append(maxv-minv+1)
        orfspan_length_mean = np.mean(orfspanlis)
        orfspan_length_median = np.median(orfspanlis)

        # intron length (mean)
        intron_list = []
        for v in feature_adj.values():
            if len(v) > 1:
                for i in range(len(v)-1):
                    intron_list.append(abs(int(v[i+1][0])-int(v[i][1])))
        intron_length = sum(intron_list) / len(intron_list)


    intact_orf,pseudogene = 0,0

    if ref != 'None':
        fragment = get_mRNA(gff_dict, get_genome(ref))
        all_mRNA = extract_mRNA(fragment)

        for k,v in all_mRNA.items():
            protein = ''
            for i in range(0, len(v), 3):
                codon = v[i:i+3]
                if codon in code:
                    protein += code[codon]
                else:
                    protein += 'X'

            judge_list = list(filter(None,protein.split('*')))
            if len(judge_list) > 1:
                pseudogene += 1

        intactORF_list = []
        for k,v in all_mRNA.items():
            if v[-3:] in stop_codon:
                intactORF_list.append(k)


        for k,v in gff_dict.items():
            if k not in intactORF_list:
                if k.endswith('-'):
                    if int(v[0][1]) > 3:
                        v[0][1] = int(v[0][1]) - 3
                else:
                    v[-1][2] = int(v[-1][2]) + 3


        orf_RNA = extract_mRNA(get_mRNA(gff_dict, get_genome(ref)))

        #intactorf = open('intact-ORF.gene.lis','w')
        for k,v in orf_RNA.items():
            if v[0:3] == 'ATG':
                if v[-3:] in stop_codon:
                    intact_orf += 1
                    #intactorf.write(k.split()[0]+'\n')

    #print('species','total genes','intact-ORF genes', 'intact-ORF proportion', 'pre-terminated genes', 'pre-terminated proportion', 'single exon genes', 'single exon proportion', 'gene length avg (bp)', 'gene length median (bp)', 'orf span length avg (bp)', 'orf span length median (bp)','coding region length avg (bp)', 'coding region length median (bp)','exon length avg (bp)','exons per gene','intron length avg (bp)',sep='\t')

    #print(name,total_gene, intact_orf, format(intact_orf/total_gene,'.2f'), pseudogene, format(pseudogene/total_gene,'.4f'), single_exon_gene, format(single_exon_gene/total_gene,'.2f'), format(gene_length_mean,'.2f'), format(locus_length_median,'.2f'),format(orfspan_length_mean,'.2f'), format(orfspan_length_median,'.2f'), format(mRNA_length_mean,'.2f'), format(mRNA_length_median,'.2f'), format(exon_length_mean,'.2f'), format(exon_per_gene,'.2f'), format(intron_length,'.2f'), sep='\t')

    print('species\t'+str(name))
    print('total genes\t'+str(total_gene))
    print('intact-ORF number\t'+str(intact_orf))
    print('intact-ORF rate\t'+str(format(intact_orf/total_gene,'.2f')))
    print('pre-terminated genes\t'+str(pseudogene))
    print('pre-terminated rate\t'+str(format(pseudogene/total_gene,'.4f')))
    print('mono-exonic genes\t'+str(single_exon_gene))
    print('mono-exonic gene rate\t'+str(format(single_exon_gene/total_gene,'.2f')))
    print('gene length avg (bp)\t'+str(format(gene_length_mean,'.2f')))
    print('gene length median (bp)\t'+str(format(locus_length_median,'.2f')))
    print('orf span length avg (bp)\t'+str(format(orfspan_length_mean,'.2f')))
    print('orf span length median (bp)\t'+str(format(orfspan_length_median,'.2f')))
    print('coding region length avg (bp)\t'+str(format(mRNA_length_mean,'.2f')))
    print('coding region length median (bp)\t'+str(format(mRNA_length_median,'.2f')))
    print('exon length avg (bp)\t'+str(format(exon_length_mean,'.2f')))
    print('exons per gene\t'+str(format(exon_per_gene,'.2f')))
    print('intron length avg (bp)\t'+str(format(intron_length,'.2f')))




if __name__ == '__main__':
    main()

