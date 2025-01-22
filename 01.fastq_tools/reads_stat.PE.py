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
parser.add_argument('--fq1', metavar='fastq',type=str,required=True,help='Please input the fastq(gzip) format')
parser.add_argument('--fq2', metavar='fastq',type=str,required=True,help='Please input the fastq(gzip) format')
parser.add_argument('--sample',metavar='sample_ID',type=str,help='Set the sample ID')
parser.add_argument('--land',metavar='land_ID',type=str,help='Set the land ID')
args = parser.parse_args()
#=========================parameters end=============================================================================

fq1,fq2,sample,land=[args.fq1,args.fq2,args.sample,args.land]


f1 = gzip.open(fq1,'rt')
f2 = gzip.open(fq2,'rt')

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



r1_len,r1_reads,r1_total,gc1_rate,max1_read_len,min1_read_len = get_main(f1)
r2_len,r2_reads,r2_total,gc2_rate,max2_read_len,min2_read_len = get_main(f2)

max_len = get_than1(max1_read_len,max2_read_len)
min_len = get_than2(min1_read_len,min2_read_len)

#print('Sample','read_len','read1_num','read2_num','total_base','GC_rate',sep='\t')
print(sample, land, max_len, min_len, r1_len/r1_reads, r1_reads+r2_reads, r1_total + r2_total, format((gc1_rate + gc2_rate)/2,'.3f'), sep='\t')

