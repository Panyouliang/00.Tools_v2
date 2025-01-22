#! /usr/bin/env python3
# -*- coding:utf-8 -*-
#
import sys,gzip,re,argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================================================
This is a script for filtering the CNV information in VCF format from lumpyexpress solfware.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2022-04-22, yyyy-mm-dd
==========================================================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('--vcf', metavar=': <vcf.gz>',type=str,required=True,help='Please input the vcf(gzip) format')
parser.add_argument('--AN',metavar=': <integer>',type=int,help='mininum of alleles in called genotypes',default=1)
parser.add_argument('--AC',metavar=': <integer>',type=int,help='mininum of the varians alleles count of genotypes',default=1)
parser.add_argument('--qual',metavar=': <float>',type=float,help='low quality threshold (default: [5])',default=5)
parser.add_argument('--svlen',metavar=': <integer>',type=int,help='low sequence length(bp) of CNV (default: [100] bp)',default=100)
parser.add_argument('--out',metavar=': <name>',type=str,required=True,help='output file name')
args = parser.parse_args()
#=========================parameters end=============================================================================
vcf_path,ancut,accut,svlength,output,qua = [args.vcf,args.AN,args.AC,args.svlen,args.out,args.qual]
f_out = gzip.open(output+'.vcf.gz','wt')


F = gzip.open(vcf_path,'rt')

vcf_dict = {}
for line in F:
    line = line.strip()
    if line[1] == '#':
        f_out.write(line+'\n')
    elif line[1] == 'C':
        f_out.write('##INFO=<ID=SF,Number=.,Type=String,Description="Source File (index to sourceFiles, f when filtered)">'+'\n'+'##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes">'+'\n'+'##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'+'\n')
        f_out.write('\t'.join(x for x in line.split())+'\n')
    elif line.split()[4][0] == 'N' or line.split()[4][-1] == 'N':
        pass
    else:
        line = line.split()
        lstr = re.split(r';',line[7])
        typer = lstr[0].split('=')[1]
        cnv = line[0]+'-'+line[1]+'-'+typer
        line[7] = lstr
        vcf_dict[cnv] = line

for key,value in vcf_dict.items():
    AN,AC,SF = 0,0,[]
    for line in value[9:]:
        sv_type = re.split(r':',line)
        if sv_type[0] == '0/0':
            AN += 2
        elif sv_type[0] == '0/1':
            AN += 2
            AC += 1
            SF.append(value[9:].index(line))
        elif sv_type[0] == '1/1':
            AN += 2
            AC += 2
            SF.append(value[9:].index(line))
    value[7].insert(5,'SF='+','.join(str(x) for x in SF))
    value[7].insert(0,'AN='+str(AN))
    value[7].insert(0,'AC='+str(AC))
    if AC >= accut and AN >= ancut:
        if float(value[5]) >= qua:
            start = int(value[1])
            pend = int(value[7][4].split('=')[1])
            if (pend - start + 1) >= svlength:
                f_out.write('\t'.join(x for x in value[0:7])+'\t'+';'.join(x for x in value[7])+'\t'+'\t'.join(x for x in value[8:])+'\n')



