#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt

import argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Cyclone raw data statistic''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('--fq', metavar='fastq',type=str,required=True,help='Please input the fastq(gzip) format')
parser.add_argument('--sample',metavar='sample_ID',type=str,help='Set the sample ID')
parser.add_argument('--species', metavar='species', required=True, default='name', type=str,help='species')
args = parser.parse_args()
#=========================parameters end=============================================================================




def readfq(F,sample,species):

    reads, length, rlist, score = 0, 0, [], []
    G_count, C_count, N_count = 0, 0, 0
    line_num = 0

    k10,k20,k30,k40,k50,k80 = 0,0,0,0,0,0

    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    for line in f:
        line_num += 1
        line = line.strip().upper()
        if line_num %4 == 1:
            vscore = float(line.split('_')[-1])
            score.append(vscore)
        if line_num %4 == 2:
            reads += 1
            length += len(line)
            rlist.append(len(line))

            G_count += line.count('G')
            C_count += line.count('C')
            N_count += line.count('N')
            if len(line) > 10000:
                k10 += len(line)
            if len(line) > 20000:
                k20 += len(line)
            if len(line) > 30000:
                k30 += len(line)
            if len(line) > 40000:
                k40 += len(line)
            if len(line) > 50000:
                k50 += len(line)
            if len(line) > 80000:
                k80 += len(line)


    #=====#statistic value #=====#
    GC_rate = (G_count + C_count)/(length - N_count)
    max_read = max(rlist)
    min_read = min(rlist)

    rlist.sort(reverse=True)

    nprlis = np.array(rlist)
    rmean = np.mean(nprlis)
    rmedian = np.median(nprlis)
    k10rate = format(k10/length*100,'.2f')
    k20rate = format(k20/length*100,'.2f')
    k30rate = format(k30/length*100,'.2f')
    k40rate = format(k40/length*100,'.2f')
    k50rate = format(k50/length*100,'.2f')
    k80rate = format(k80/length*100,'.2f')

    npscore = np.array(score)

    score_mean,score_median = np.mean(npscore), np.median(npscore)
    max_score,min_score = max(score),min(score)
    #print(score_mean, score_median, max_score, min_score)



    N50_length = length/2
    Nnum = 0
    for i in rlist:
        Nnum += i
        if Nnum >= N50_length:
            N50_reads = i
            break


    print(sample,reads,max_read,min_read,k10rate,k20rate,k30rate,k40rate,k50rate,N50_reads,format(rmean,'.2f'),format(rmedian,'.2f'),length,format(GC_rate,'.2f'),score_mean,score_median,max_score,min_score,sep='\t')


    #plt.hist(nprlis, bins=1000, edgecolor='black')
    #plt.xlim((-1,100000))
    #plt.xlabel('read length')
    #plt.ylabel('Frequency')
    #plt.title('Histogram of Cyclone')

    #plt.savefig(sample+'.pdf', bbox_inches='tight')


def main():
    readfq(args.fq, args.sample, args.species)


if __name__ == '__main__':
    main()






