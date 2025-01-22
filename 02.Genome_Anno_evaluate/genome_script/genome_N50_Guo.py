#!/ldfssz1/ST_DIVERSITY/PUB/USER/guoqunfei/local/Python-3.8.1/bin/python
# -*- coding:utf-8 -*-


import gzip,argparse


def getting_N_dt(len_l,cutoff):

    len_l.sort(reverse=True)
    Total_size,Longest = str(sum(len_l)),str(len_l[0])
    len_cutoff = str(len([x for x in len_l if x >= cutoff]))
    len2000 = str(len([x for x in len_l if x >= 2000]))

    N_dt = {}
    for i in range(10,100,10):
        num = 0
        for j in range(len(len_l)):
            num += len_l[j]
            if num >= int(Total_size)*i/100:
                N_dt[i] = [str(len_l[j]),str(j+1)]
                break
    return Longest,Total_size,len_cutoff,len2000,N_dt


def deal_scf(scf,scfGC,scfN,scflen_l,ctgGC,ctglen_l,cutoff):

    scf_seq = ''.join(scf).upper()
    scfGC += scf_seq.count('G') + scf_seq.count('C')
    scfN += scf_seq.count('N')
    scflen_l.append(len(scf_seq))

    ctg_l = [x for x in scf_seq.split('N') if len(x) >= cutoff and len(x) > 0]
    ctgGC += sum([x.count('G') for x in ctg_l]) + sum([x.count('C') for x in ctg_l])
    ctglen_l.extend([len(x) for x in ctg_l])

    return scfGC,scfN,scflen_l,ctgGC,ctglen_l


def getting_scf_ctg(fa_f,cutoff):

    if fa_f[-3:] == ".gz":
        f_in = gzip.open(fa_f,'rt')
    else:
        f_in = open(fa_f,'r')

    scf,scfGC,scfN,scflen_l,ctgGC,ctglen_l = [],0,0,[],0,[]
    for line in f_in:
        l = line.split()

        if len(l) == 0:
            continue

        if l[0][0] == '>':
            if scf:
                scfGC,scfN,scflen_l,ctgGC,ctglen_l = deal_scf(scf,scfGC,scfN,scflen_l,ctgGC,ctglen_l,cutoff)
            scf = []
        else:
            scf.append(l[0])
    scfGC,scfN,scflen_l,ctgGC,ctglen_l = deal_scf(scf,scfGC,scfN,scflen_l,ctgGC,ctglen_l,cutoff)

    f_in.close()

    return scflen_l,ctglen_l,scfGC,scfN,ctgGC


def main():

#----------------------------------parameters start-------------------------------------------------------------------------
    version = "v1.0"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
=====================================================================
This is a script for counting the genome information in fasta format.

Author: Qunfei Guo, guoqunfei@genomics.cn
Version: v1.0
Date: 2020-08-04, yyyy-mm-dd
=====================================================================''')
    parser.add_argument('-v', '--version',action='version', version=version)
    parser.add_argument('--f',metavar='fasta',type=str,required=True,help='Please input the genome in fasta(gzip) format')
    parser.add_argument('--c',metavar='cutoff',type=int,help='Set the minimum length of the contig.default=100',default=100)
    args = parser.parse_args()
#----------------------------------parameters end----------------------------------------------------------------------------

    fa_f,cutoff = [args.f,args.c]

    scflen_l,ctglen_l,scfGC,scfN,ctgGC = getting_scf_ctg(fa_f,cutoff)
    scfLongest,scfTotal_size,scflen_cutoff,scflen2000,scfN_dt = getting_N_dt(scflen_l,cutoff)
    ctgLongest,ctgTotal_size,ctglen_cutoff,ctglen2000,ctgN_dt = getting_N_dt(ctglen_l,cutoff)

    print('\t' + 'scaffold' + '\t\t' + 'contig')
    print('\t'.join(['\tlength(bp)','number','length(bp)','number']))
    for i in range(90,0,-10):
        print('\t'.join(['N' + str(i)] + scfN_dt[i] + ctgN_dt[i]))
    print('\t'.join(['Longest',scfLongest,'\t' + ctgLongest]))
    print('\t'.join(['Total_size',scfTotal_size,'\t' + ctgTotal_size]))
    print('\t'.join(['number>='+str(cutoff)+'bp',scflen_cutoff,'\t' + ctglen_cutoff]))
    print('\t'.join(['number>=2kbp',scflen2000,'\t' + ctglen2000]))
    print('\t'.join(['GC_rate','%.3f' % float(scfGC/(int(scfTotal_size)-scfN)),'%.3f' % float(ctgGC/int(ctgTotal_size))]))

if __name__ == "__main__":
    main()
