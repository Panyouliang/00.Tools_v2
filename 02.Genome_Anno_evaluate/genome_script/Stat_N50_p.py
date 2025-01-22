#! /usr/bin/env python3

import gzip,argparse


def scf_ctg_message_get(files,num):

    scf, scf_seq, scf_total, l, N_total = [],'', 0, 0, 0
    ctg, ctg_seq, ctg_total,G, C = [], '', 0, 0, 0


    if files[-3] == '.gz':
        f_in = gzip.open(files,'rt')
    else:
        f_in = open(files,'r')

    for line in f_in:
        if line[0] == '>':
            scf.append(l)
            l = 0
            for i in ctg_seq.split('N'):
                if len(i) >= int(num):
                    ctg_total += len(i)
                    ctg.append(len(i))
            ctg_seq = ''
        else:
            line = line.upper()
            l += len(line.strip('\n'))
            G += line.count('G')
            C += line.count('C')
            N_total += line.count('N')
            ctg_seq += line.strip('\n')
            scf_total += len(line.strip('\n'))

    for i in ctg_seq.split('N'):
        if len(i) > int(num):
            ctg_total += len(i)
            ctg.append(len(i))
    scf.append(l)
    GC_num = G + C
    scf = sorted(scf,key=lambda x:x,reverse=True)
    ctg = sorted(ctg,key=lambda x:x,reverse=True)

    return scf, ctg, GC_num, scf_total,ctg_total, N_total

def stats(files, num):
    up, up2000 = 0, 0
    for i in files:
        if i > int(num):
            up += 1
            if i >= 2000:
                up2000 += 1
    return up, up2000

def N50_get(files, total):
    seq = {}
    for n in range(10,100,10):
        length, num = 0, 0
        seq["N"+str(n)] = ''
        for i in files:
            length += i
            num += 1
            if length >= total*n/100:
                seq["N"+str(n)] = str(i)+"\t"+str(num)
                break
    return seq


def main():

#-------------------------parameters start---------------------------------------------------------------------------
    version = "v1.0"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================
This is a script for counting the genome information in fasta format.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2021-09-03, yyyy-mm-dd
==========================================================''')
    parser.add_argument('-v', '--version',action='version',version=version)
    parser.add_argument('--f',metavar='fasta',type=str,required=True,help='Please input the genome in fasta(gzip) format')
    parser.add_argument('--c',metavar='cutoff',type=int,help='Set the minimum length of the contig.default=100',default=100)
    args = parser.parse_args()
#------------------------parameters end-----------------------------------------------------------------------------

    path,cutoff = [args.f,args.c]


    scf,ctg,GC_num,scf_total,ctg_total,N_total = scf_ctg_message_get(path,cutoff)

    scf_up,scf_up2000 = stats(scf,cutoff)
    ctg_up,ctg_up2000 = stats(ctg,cutoff)

    scf_seq = N50_get(scf,scf_total)
    ctg_seq = N50_get(ctg,ctg_total)

    print("\t"+"scaffold"+"\t\t"+"contig")
    print("\t"+"length"+"\t"+"number"+"\t"+"length"+"\t"+"number")
    for n in range(90,0,-10):
        print("N"+str(n)+"\t"+scf_seq["N"+str(n)]+"\t"+ctg_seq["N"+str(n)])
    print("Longest"+"\t"+str(scf[0])+"\t"+str(ctg[0]))
    print("Total_size"+"\t"+str(scf_total)+"\t"+str(ctg_total))
    print("number>="+str(cutoff)+"bp"+"\t"+str(scf_up)+"\t"+str(ctg_up))
    print("number>=2kb"+"\t"+str(scf_up2000)+"\t"+str(ctg_up2000))
    print("GC_rate"+"\t"+str(format(GC_num/(scf_total - N_total),'.3f'))+'\t'+str(format(GC_num/ctg_total,'.3f')))

if __name__ == "__main__":
    main()
