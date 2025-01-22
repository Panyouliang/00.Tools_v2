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
#parser.add_argument('--type',type=str,required=False,default='PAIRED',help='Please input PAIRED|SINGLE types')
parser.add_argument('--filelist',type=str,required=True,help='Please input the SOAPnuke result statistic files lists')
args = parser.parse_args()
#=========================parameters end=============================================================================


def get_sp(F,PE):
    fqQC = {}
    with open(F,'r') as f:
        for line in f:
            if line.strip():
                line = line.strip().split('\t')
                if PE == 'PAIRED':
                    fqQC[line[0].strip()] = [line[1].strip(),line[2].strip(),line[3].strip(),line[4].strip()]
                else:
                    fqQC[line[0].strip()] = [line[1].strip(),line[2].strip()]
    return fqQC


print('Run_ids','raw_reads_length(fq1)','raw_reads_length(fq2)','clean_reads_length(fq1)','clean_reads_length(fq2)','raw_reads_number','clean_reads_number','raw_total_base','clean_total_base','filter_base(%)','Adapter_related(%)','raw_data_Q20(%)','clean_data_Q20(%)','raw_data_Q30(%)','clean_data_Q30(%)',sep='\t')

#=====================
def get_PE_stat(D):
    lis = D
    raw_read1L,raw_read2L,cl_read1L,cl_read2L = lis['Read length'][0],lis['Read length'][2],lis['Read length'][1],lis['Read length'][3]

    raw_readN,cl_readN = int(lis['Total number of reads'][0].split()[0]) + int(lis['Total number of reads'][2].split()[0]),int(lis['Total number of reads'][1].split()[0]) + int(lis['Total number of reads'][3].split()[0])

    raw_total_B,cl_total_B = int(lis['Total number of bases'][0].split()[0]) + int(lis['Total number of bases'][2].split()[0]),int(lis['Total number of bases'][1].split()[0]) + int(lis['Total number of bases'][3].split()[0])

    filter_b,ad_content = (int(lis['Number of filtered bases (%)'][0].split()[0]) + int(lis['Number of filtered bases (%)'][2].split()[0]))/raw_total_B,(int(lis['Reads related to Adapter and Trimmed (%)'][0].split()[0]) + int(lis['Reads related to Adapter and Trimmed (%)'][2].split()[0]))/raw_total_B

    raw_q20,raw_q30 = (int(lis['Number of base calls with quality value of 20 or higher (Q20+) (%)'][0].split()[0]) + int(lis['Number of base calls with quality value of 20 or higher (Q20+) (%)'][2].split()[0]))/raw_total_B,(int(lis['Number of base calls with quality value of 20 or higher (Q20+) (%)'][1].split()[0]) + int(lis['Number of base calls with quality value of 20 or higher (Q20+) (%)'][3].split()[0]))/cl_total_B

    cl_q20,cl_q30 = (int(lis['Number of base calls with quality value of 30 or higher (Q30+) (%)'][0].split()[0]) + int(lis['Number of base calls with quality value of 30 or higher (Q30+) (%)'][2].split()[0]))/raw_total_B,(int(lis['Number of base calls with quality value of 30 or higher (Q30+) (%)'][1].split()[0]) + int(lis['Number of base calls with quality value of 30 or higher (Q30+) (%)'][3].split()[0]))/cl_total_B

    outlis = [str(raw_read1L),str(raw_read2L),str(cl_read1L),str(cl_read2L),str(raw_readN),str(cl_readN),str(raw_total_B),str(cl_total_B),str(format(filter_b,'.2f')),str(format(ad_content,'.2f')),str(format(raw_q20,'.2f')),str(format(raw_q30,'.2f')),str(format(cl_q20,'.2f')),str(format(cl_q30,'.2f'))]

    return outlis

def get_SE_stat(D):
    lis = D
    raw_read1L,cl_read1L = lis['Read length'][0],lis['Read length'][1]

    raw_readN,cl_readN = int(lis['Total number of reads'][0].split()[0]),int(lis['Total number of reads'][1].split()[0])

    raw_total_B,cl_total_B = int(lis['Total number of bases'][0].split()[0]),int(lis['Total number of bases'][1].split()[0])

    filter_b,ad_content = int(lis['Number of filtered bases (%)'][0].split()[0])/raw_total_B,int(lis['Reads related to Adapter and Trimmed (%)'][0].split()[0])/raw_total_B

    raw_q20,raw_q30 = int(lis['Number of base calls with quality value of 20 or higher (Q20+) (%)'][0].split()[0])/raw_total_B,int(lis['Number of base calls with quality value of 20 or higher (Q20+) (%)'][1].split()[0])/cl_total_B

    cl_q20,cl_q30 = int(lis['Number of base calls with quality value of 30 or higher (Q30+) (%)'][0].split()[0])/raw_total_B,int(lis['Number of base calls with quality value of 30 or higher (Q30+) (%)'][1].split()[0])/cl_total_B

    outlis = [str(raw_read1L),0,str(cl_read1L),0,str(raw_readN),str(cl_readN),str(raw_total_B),str(cl_total_B),str(format(filter_b,'.2f')),str(format(ad_content,'.4f')),str(format(raw_q20,'.2f')),str(format(raw_q30,'.2f')),str(format(cl_q20,'.2f')),str(format(cl_q30,'.2f'))]

    return outlis



def main():
    filelis = args.filelist
    for i in open(filelis,'r'):
        i = i.strip().split()
        ids = i[0].split('/')[-2]
        if i[1] == 'PAIRED':
            outlis = get_PE_stat(get_sp(i[0],i[1]))
            print(ids,'\t'.join(str(x) for x in outlis),sep='\t')
        else:
            outlis = get_SE_stat(get_sp(i[0],i[1]))
            print(ids,'\t'.join(str(x) for x in outlis),sep='\t')


if __name__ == '__main__':
    main()

