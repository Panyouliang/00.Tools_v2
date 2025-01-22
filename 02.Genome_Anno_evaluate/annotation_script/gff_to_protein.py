#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import gzip,sys,re,argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================
This is a script for extracted pep and mRNA from genome.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2022-04-21, yyyy-mm-dd
==========================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('-gff',type=str,required=False,default='PAIRED',help='Please input the GFF annotation files')
parser.add_argument('-ref',type=str,required=True,help='Please input the reference genome')
parser.add_argument('-sp',type=str,required=True,help='Please input the species name')
parser.add_argument('-IDtype',type=str,required=True,default='T',help='Please input the G|T, "G" represent geneID, "T" represent transcriptID')

args = parser.parse_args()
#=========================parameters end=============================================================================


def get_genome(files):
    g_dict,seq,name = {},[],''
    if files[-3:] == '.gz':
        F = gzip.open(files,'rt')
    else:
        F = open(files)
    fileiter = (x.strip() for x in F)
    for line in fileiter:
        if line.startswith('>'):
            if name != '':
                g_dict[name] = ''.join(seq)
            name = line.split(' ')[0]
            seq = []
        else:
            seq.append(line)
    g_dict[name] = ''.join(seq)
    return g_dict


def gff_get(file2):
    if file2[-3:] == '.gz':
        F = gzip.open(file2,'rt')
    else:
        F = open(file2)
    gff_d,g2t,d = {},{},[]
    for line in F:
        if line[0] != '#':
            line = line.strip().split()
            if len(line) != 0:
                if line[2] == 'mRNA':
                    if 'Parent' in line[8]:
                        idx = line[8].find('Parent')
                    else:
                        idx = line[8].find('gene_id')

                    P = line[8][idx:].split(';')[0].split('=')[1]
                    tid = line[8].split(';')[0].split('=')[1]
                    gid = tid +line[6]
                    g2t[tid] = P
                    d = []
                elif line[2] == 'CDS':
                    pos = ['>'+line[0],line[3],line[4]]
                    d.append(pos)
                    gff_d[gid] = d

    gff_d2 = {}
    for k,v in gff_d.items():
        gff_d2[k] = sorted(v,key = (lambda x:int(x[1])))
    return gff_d2,g2t


def get_fragment(gff_d,g_seq):
    seq = {}
    for k,v in gff_d.items():
        if k:
            gid = k
            seq[gid] = ''
            fragment = ''
            for line in v:
                fragment = g_seq[line[0]][int(line[1])-1:int(line[2])]
                seq[gid] += fragment
    return seq


def rev(seq):
    base = {
            'A':'T','T':'A','G':'C','C':'G','N':'N','n':'n','a':'t','t':'a','c':'g','g':'c'
            }
    seq_list = list(reversed(seq))
    seq_rev = [base[k] for k in seq_list]
    seq_list = ''.join(seq_rev)
    return seq_list

def break_sequence(sequence, length):
    subsequences = []
    for i in range(0, len(sequence), length):
        subsequence = sequence[i:i+length]
        subsequences.append(subsequence)
    return subsequences



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
        'TAA':'U','TAG':'U','TGA':'U'
        }

stop_code = {'TAA':'','TAG':'','TGA':''}


def main():
    gff,ref,sp,IDtype = args.gff,args.ref,args.sp,args.IDtype

    g_data = get_genome(ref)
    gff_data,g2t = gff_get(gff)
    fragment = get_fragment(gff_data,g_data)

    g_mRNA = {}
    for k,v in fragment.items():
        if k.endswith('-'):
            frag = rev(v)
            g_mRNA[k] = frag
        elif k.endswith('+'):
            g_mRNA[k] = v

    fcds = open(sp+'.cds.fa','w')
    for k,v in g_mRNA.items():
        if IDtype == 'G':
            fcds.write('>'+g2t[k[0:-1]]+'\n')
        if IDtype == 'T':
            fcds.write('>'+k[0:-1]+'\n')
        fcds.write(v+'\n')

    fcds.close()

    fpep = open(sp+".protein.fa",'w')
    for k,v in g_mRNA.items():
        if IDtype == 'G':
            fpep.write('>'+g2t[k[0:-1]]+'\n')
        if IDtype == 'T':
            fpep.write('>'+k[0:-1]+'\n')

        protein = ''
        subseq = break_sequence(v.upper(), 3)

        for codon in subseq[0:-1]:
            if codon in code:
                protein += code[codon]
            else:
                protein += 'X'

        if subseq[-1] in stop_code:
            protein += stop_code[subseq[-1]]
        else:
            if subseq[-1] in code:
                protein += code[subseq[-1]]
            else:
                protein += ''
        fpep.write(protein+'\n')

    fpep.close()

if __name__ == '__main__':
    main()

