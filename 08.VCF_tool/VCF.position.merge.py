#! /usr/bin/env python3
# -*- coding:utf-8 -*-
#

import sys,re,gzip,argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================
This is a script for merge the CNV position of reference in VCF format.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2022-04-21, yyyy-mm-dd
==========================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('--vcf', metavar='vcf.gz',type=str,required=True,help='Please input the vcf(gzip) format')
parser.add_argument('--c',metavar='(default: [1])',type=int,help='Set the minimum number of alleles in called genotypes',default=1)
parser.add_argument('--out',metavar='name',type=str,required=True,help='Set output name')
args = parser.parse_args()
#=========================parameters end=============================================================================

vcf_path,cutoff,outpath=[args.vcf,args.c,args.out]
f_out = gzip.open(outpath+'.vcf.gz','wt')



F = gzip.open(vcf_path,'rt')
cut = int(cutoff)
pos_dict = {}
for line in F:
    if line[1] == '#':
        f_out.writelines(line.strip()+'\n')
    elif line[1] == 'C':
        f_out.write('\t'.join(str(x) for x in line.strip().split())+'\n')
    else:
        SF = []
        line = line.strip().split()
        if line[4] == '<DUP>,<DEL>' or line[4] == '<DEL>,<DUP>':
            pass
        else:
            lstr = re.split(r';',line[7])
            AC = int(lstr[0].split('=')[1])
            AN = int(lstr[1].split('=')[1])
            pend = lstr[2].split('=')[1]
            sf = re.split(r',',lstr[4].split('=')[1])
            svlen = int(pend)-int(line[1])+1
            for i in sf:
                SF.append(int(i))


            p_info = [line[1],pend,svlen,AC,AN,SF]
            vcf_line = []
            for i in line:
               vcf_line.append(i)


            p_info.append(vcf_line)
            if line[4] == '<DEL>':
                later = '-'
                if line[0]+'DEL' not in pos_dict.keys():
                    chrm = [p_info]
                    pos_dict[line[0]+'DEL'] = chrm
                else:
                    chrm.append(p_info)
                    pos_dict[line[0]+'DEL'] = chrm

            elif line[4] == '<DUP>':
                later = ''
                if line[0]+'DUP' not in pos_dict.keys():
                    chrm = [p_info]
                    pos_dict[line[0]+'DEL'] = chrm
                else:
                    chrm.append(p_info)
                    pos_dict[line[0]+'DUP'] = chrm



for key,value in pos_dict.items():
    index = []
    for i in range(len(value)):
        index.append(int(i))
    
    if len(index) < 2:
        #pass
        if len(value[0][5]) >= cut: 
            #print(key,'\t'.join(str(x) for x in value[0][:3]),','.join(str(x) for x in value[0][3]),sep='\t')
            f_out.write('\t'.join(str(x) for x in value[0][6][0:])+'\n')
        
    else:
        t = int(len(index)-1)
        index.pop(t)
        
        out_index = []
        dele_index = []
        for ind in index:
            ind = int(ind)

            covlen = int(value[0+ind][1])-int(value[1+ind][0])
            h_svlen = (int(value[0+ind][2])+int(value[1+ind][2]))/4

            if covlen >= h_svlen:
                dele_index.append(ind)
                value[1+ind][1] = value[0+ind][1]
                value[1+ind][2] = int(value[1+ind][1])-int(value[1+ind][0])+1

                for sample in value[0+ind][5]:
                    value[1+ind][5].append(sample)
                    value[1+ind][6][int(sample)+9] = value[0+ind][6][int(sample)+9]
            
            dedup = set(value[1+ind][5])
            value[1+ind][5] = sorted(dedup)

            vcf_info = re.split(r';',value[1+ind][6][7])
            vcf_info[0] = 'AC='+str(value[1+ind][3]+value[0+ind][3])
            vcf_info[1] = 'AN='+str(value[1+ind][4]+value[0+ind][3])
            vcf_info[2] = 'END='+str(value[1+ind][1])
            vcf_info[4] = 'SF='+','.join(str(x) for x in value[1+ind][5])
            vcf_info[5] = 'SVLEN='+later+str(value[1+ind][2])
            value[1+ind][6][7] = ';'.join(x for x in vcf_info)


        for vcf_del in reversed(dele_index):
            value.pop(vcf_del)

        for i in range(len(value)):
            out_index.append(int(i))

        for o in out_index:
            if len(value[o][5]) >= cut:
                #print(key,'\t'.join(str(x) for x in value[o][:3]),','.join(str(x) for x in value[o][3]),sep='\t')
                f_out.write('\t'.join(x for x in value[o][6][0:])+'\n')













