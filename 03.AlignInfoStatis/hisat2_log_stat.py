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
parser.add_argument('-log',type=str,required=True,help='Please input the samtools flagstat result')
parser.add_argument('-seqtype',type=str,required=True,help='Please input the PAIRED or SINGLE')

args = parser.parse_args()

#=========================parameters end=============================================================================


def PE(F):
    with open(F,'r') as f:
        lines = f.readlines()
        total = int(lines[1].split()[0])
        primary_rate = float(lines[14].split()[0].strip('%'))/100
        exactrd = int(lines[3].split()[0])
        seconrd = int(lines[4].split()[0])
        dis_uniq = int(lines[7].split()[0])
        dis_single_uniq = int(lines[12].split()[0])
        properly_rate = (exactrd+seconrd)/(primary_rate*total) ### properly rate: properly reads / total reads
        uniq_primary_rate = (exactrd*2+dis_uniq*2+dis_single_uniq)/(primary_rate*total*2) ### unique properly rate: unique reads / primary reads
        uniq_properly_rate= exactrd/(exactrd+seconrd) ### unique mapped rate: unique reads / properly reads
    return primary_rate, properly_rate, uniq_primary_rate, uniq_properly_rate


def SE(F):
    with open(F,'r') as f:
        lines = f.readlines()
        total = int(lines[1].split()[0])
        primary_rate = float(lines[6].split()[0].strip('%'))/100
        exactrd = int(lines[4].split()[0])
        seconrd = int(lines[5].split()[0])
        properly_rate = (exactrd+seconrd)/(exactrd+seconrd)
        uniq_primary_rate = exactrd*2/(exactrd*2+seconrd*2)
        uniq_properly_rate= exactrd/(exactrd+seconrd)

    return primary_rate, properly_rate, uniq_primary_rate, uniq_properly_rate





def main():
    run,flag,seqtype = args.run,args.log,args.seqtype
    if seqtype == 'PAIRED':
        primary_rate, properly_rate, uniq_primary_rate, uniq_properly_rate = PE(flag)
    else:
        primary_rate, properly_rate, uniq_primary_rate, uniq_properly_rate = SE(flag)


    #print('run', 'primary/total rate', 'properly/primary rate', 'unique/primary rate', 'unique/properly rate', sep='\t')
    print(run, str(format(primary_rate,'.4f')), str(format(properly_rate,'.4f')), str(format(uniq_primary_rate,'.4f')), str(format(uniq_properly_rate,'.4f')), sep='\t')



if __name__ == '__main__':
    main()

