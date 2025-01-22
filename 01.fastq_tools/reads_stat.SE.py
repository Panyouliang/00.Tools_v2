#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,re,gzip,argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================
This is a script for counting the genome information in fasta format.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2022-04-21, yyyy-mm-dd
==========================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('--fq', metavar='fastq',type=str,required=True,help='Please input the fastq(gzip) format')
parser.add_argument('--sample',metavar='sample_ID',type=str,help='Set the sample ID')
args = parser.parse_args()
#=========================parameters end=============================================================================

fq,sample=[args.fq,args.sample]


f = gzip.open(fq,'rt')
#f2 = gzip.open(fq2,'rt')

def get_main(inp):
    l_n = 0
    reads = 0
    total_base = 0
    g_count,c_count,n_total = 0,0,0
    length = 0
    read_length_lists = []
    for line in inp:
        l_n += 1
        if l_n %4 == 2:
            reads += 1
            line = line.strip('\n')
            line = line.upper()
            length += len(line)
            read_length_lists.append(len(line))
            g_count += line.count('G')
            c_count += line.count('C')
            n_total += line.count('N')
            total_base += len(line)

    Reads = reads
    gc_rate = (g_count + c_count) / (total_base - n_total)
    max_read_len = max(read_length_lists)
    min_read_len = min(read_length_lists)
    return length,Reads,total_base,gc_rate,max_read_len,min_read_len

def get_than1(r1,r2):
    if r1 > r2:
        value1 = r1
    else:
        value1 = r2
    return value1

def get_than2(t1,t2):
    if t1 < t2:
        value2 = t1
    else:
        value2 = t2
    return value2



r_len,r_reads,r_total,gc_rate,max_read_len,min_read_len = get_main(f)


print(sample, max_read_len, min_read_len, format(r_len/r_reads,'.2f'), r_reads, r_total, format(gc_rate,'.3f'), sep='\t')

