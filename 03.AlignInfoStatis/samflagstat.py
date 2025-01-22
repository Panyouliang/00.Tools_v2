#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,subprocess,argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================
    Author: Youliang Pan, panyouliang@genomics.cn
    Version: v1.0
    Date: 2023-06-21, yyyy-mm-dd
==========================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('-run',type=str,required=True,help='Please input run name types')
parser.add_argument('-flag',type=str,required=True,help='Please input the samtools flagstat result')
parser.add_argument('-bam',type=str,required=True,help='Please input the bamfiles')

args = parser.parse_args()

#=========================parameters end=============================================================================

run,flag,bamfile = args.run,args.flag,args.bam



with open(flag,'r') as f:
    lines = f.readlines()


    total = int(lines[1].split()[0])

    maread = int(lines[6].split()[0])
    pmread = int(lines[7].split()[0])
    ppread = int(lines[11].split()[0])


uniqs = subprocess.run(['samtools', 'view', '-c', '-q', '30', '-f 2', bamfile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

uniqs = int(uniqs.stdout.decode())

#print('run', 'primary mapped rate', 'properly paired rate', 'unique mapped rate', sep='\t')
print(run, str(pmread/total), str(ppread/pmread), str(uniqs/pmread), sep='\t')

