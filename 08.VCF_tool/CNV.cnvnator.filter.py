#! /usr/bin/env python3
# -*- coding:utf-8 -*-
#
import sys,gzip,re,argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==================================================================
This is a script for filtering the CNV information in VCF format.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2022-04-22, yyyy-mm-dd
==================================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('--vcf', metavar='<vcf.gz>',type=str,required=True,help='Please input the vcf(gzip) format')
parser.add_argument('--AN',metavar='<integer>',type=int,help='Set the minimum number of alleles in called genotypes',default=1)
parser.add_argument('--AC',metavar='<integer>',type=int,help='Set the mininum number of varians alleles count of genotypes',default=1)
parser.add_argument('--out',metavar='<output>',type=str,required=True,help='Set output name')
args = parser.parse_args()
#=========================parameters end=============================================================================
vcf_path,ANcut,ACcut,output=[args.vcf,args.AN,args.AC,args.out]

f_out = gzip.open(output+'.vcf.gz','wt')



F = gzip.open(vcf_path,'rt')
vcf_dict = {}
for line in F:
    line = line.strip()
    if line[1] == '#':
        f_out.write(line+'\n')
    elif line[1] == 'C':
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
    AN,AC,SF,n = 0,0,[],0
    for i in range(len(value[9:])):
        if value[9:][i] == '.':
            AN += 2
            value[int(i)+9] = '0/0:.'

        else:
            sv_type = re.split(r':',value[int(i)+9])
            if sv_type[0] == './1':
                AN += 2
                AC += 1
                SF.append(i)
                value[int(i)+9] = '0/1:'+sv_type[1]
            elif sv_type[0] == '0/1':
                AN += 2
                AC += 1
                SF.append(i)
            elif sv_type[0] == '1/1':
                AN += 2
                AC += 2
                SF.append(i)
            elif sv_type[0] == './2':
                AN += 2
                AC += 1
                SF.append(i)
                value[int(i)+9] = '0/2:'+sv_type[1]
            elif sv_type[0] == '0/2':
                AN += 2
                AC += 1
                SF.append(i)
            elif sv_type[0] == '2/2':
                AN += 2
                AC += 1
                SF.append(i)


    value[7][4] = 'SF='+','.join(str(x) for x in SF)
    value[7][0] = 'AC='+str(AC)
    value[7][1] = 'AN='+str(AN)
    if AC >= ACcut and AN >= ANcut:
        f_out.write('\t'.join(str(x) for x in value[0:7])+'\t'+';'.join(str(x) for x in value[7])+'\t'+'\t'.join(str(x) for x in value[8:])+'\n')



