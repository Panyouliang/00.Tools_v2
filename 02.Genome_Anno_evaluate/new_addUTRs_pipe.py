#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,gzip,argparse
import subprocess

#================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
================================================================
This is a script for evaluating the quality of genome annotation.

Author: Panyouliang, panyouliang@genomics.cn
Version: v1.0
Date: 2023-07-04, yyyy-mm-dd
================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-core', metavar='annotation file', type=str, required=True, help='Please input the annotation file(Select Top length isoforms before input)')
parser.add_argument('-isogff',  metavar='genome file', type=str, required=False, help='Please input the stringtie assemble file of ISO-seq(gff format)')
parser.add_argument('-ngsgff',  metavar='genome file', type=str, required=False, help='Please input the stringtie assemble file of RNA-seq(gff format)')
parser.add_argument('-sp', metavar='Specie name', type=str, required=False, default="husky", help='Please input the Specie name')
args = parser.parse_args()
#=================================================================================


def merge_adjacent_coordinates(coordinates):
    merged_coordinates = []
    if not coordinates:
        return merged_coordinates

    current_start, current_end = coordinates[0]
    for coord in coordinates[1:]:
        if coord[0] == current_end + 1:
            current_end = max(current_end, coord[1])
        else:
            merged_coordinates.append([current_start, current_end])
            current_start, current_end = coord

    merged_coordinates.append([current_start, current_end])

    return merged_coordinates



def stat_overlap(gff1,gff2):
    result = subprocess.run(['perl', '/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/FindOverlapAtCDSlevel.exon.pl', gff1, gff2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output = result.stdout.decode()
    #f = open('MPHA_GENOME.anno.gff.loc.gff.lis','r')
    overlapd = {}
    #for line in f:
    for line in output.splitlines():
        if line[0] != '#':
            line = line.strip().split('\t')
            pid1,pid2 = line[0],line[1]

            loc1,loc2 = float(line[5]),float(line[6])
            cod1,cod2 = float(line[10]),float(line[11])

            datalis = [loc1,loc2,cod1,cod2]
            if pid1 in overlapd.keys():
                if cod1 > 0.5 and cod2 > 0.20:
                    overlapd[pid1][pid2] = sum(datalis)
            else:
                if cod1 > 0.5 and cod2 > 0.20:
                    overlapd[pid1] = {pid2:sum(datalis)}

    core2rna_d = {}
    for k,v in overlapd.items():
        maxk = max(v, key=v.get)
        core2rna_d[k] = maxk

    return core2rna_d


def geometric_mean(lis):
    total = 1
    for i in lis:
        total *= i
    return pow(total,1/len(lis))


def core_loc_get(F):
    gff,core_cds_loc,core_exon_loc = {},{},{}

    chrs_lis = []
    locus_minus = {}
    locus_plus  = {}
    core2core = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split()
                if line[2] == 'gene':
                    pass
                elif line[2] == 'mRNA':
                    rnaID = line[8].split('=')[1].split(';')[0]
                    idx = line[8].find('Parent') if 'Parent' in line[8] else line[8].find('gene_id')
                    geneID = line[8][idx:].split('=')[1].split(';')[0]
                    core_cds_loc[rnaID+' '+line[6]] = []
                    core_exon_loc[rnaID] = []
                    core2core[rnaID] = rnaID
                    if geneID in gff.keys():
                        gff[geneID][rnaID].append(line)
                    else:
                        gff[geneID] = {rnaID:[line]}

                    '''
                    chrs = line[0]
                    chrs_lis.append(chrs)
                    if line[6] == '+':
                        if chrs in locus_plus.keys():
                            locus_plus[chrs][rnaID] = [int(line[3]),int(line[4])]
                        else:
                            locus_plus[chrs] = {rnaID:[int(line[3]),int(line[4])]}

                    if line[6] == '-':
                        if chrs in locus_minus.keys():
                            locus_minus[chrs][rnaID] = [int(line[3]),int(line[4])]
                        else:
                            locus_minus[chrs] = {rnaID:[int(line[3]),int(line[4])]}
                    '''


                else:
                    core_exon_loc[rnaID].extend([int(line[3]),int(line[4])])
                    if line[2] == 'CDS':
                        gff[geneID][rnaID].append(line)
                        core_cds_loc[rnaID+' '+line[6]].extend([int(line[3]),int(line[4])])

                        chrs = line[0]
                        chrs_lis.append(chrs)
                        if line[6] == '+':
                            if chrs in locus_plus.keys():
                                if rnaID in locus_plus[chrs].keys():
                                    locus_plus[chrs][rnaID].extend([int(line[3]),int(line[4])])
                                else:
                                    locus_plus[chrs][rnaID] = [int(line[3]),int(line[4])]

                            else:
                                locus_plus[chrs] = {rnaID:[int(line[3]),int(line[4])]}

                            locus_plus[chrs][rnaID].sort()


                        if line[6] == '-':
                            if chrs in locus_minus.keys():
                                if rnaID in locus_minus[chrs].keys():
                                    locus_minus[chrs][rnaID].extend([int(line[3]),int(line[4])])
                                else:
                                    locus_minus[chrs][rnaID] = [int(line[3]),int(line[4])]
                            else:
                                locus_minus[chrs] = {rnaID:[int(line[3]),int(line[4])]}

                            locus_minus[chrs][rnaID].sort()





    plus_loc = {}
    for chrs,v in locus_plus.items():
        plus_loc[chrs] = [x[0] for x in sorted(v.items(), key=lambda x:x[1][0],reverse=False)]


    minus_loc = {}
    for chrs,v in locus_minus.items():
        minus_loc[chrs] = [x[0] for x in sorted(v.items(), key=lambda x:x[1][0],reverse=False)]


    chrs_los = list(set(chrs_lis))


    return gff, core_cds_loc, core_exon_loc, plus_loc, minus_loc, chrs_los, core2core



def transcript_loc_get(F):
    exonloc = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] == '#':
                pass
            else:
                line = line.strip().split('\t')
                if line[2] == 'transcript' or line[2] == 'mRNA':
                    rnaID = line[8].split('=')[1].split(';')[0]
                    exonloc[rnaID] = []
                elif line[2] == 'exon':
                    exonloc[rnaID].extend([int(line[3]),int(line[4])])

    return exonloc


def cutUTR(core_cds_loc,trans_exon_loc,core2rna_d):
    UTR3,UTR5 = {},{}
    for k,v in core_cds_loc.items():
        tid,strand = k.split()[0],k.split()[1]
        cds_max,cds_min = max(v),min(v)
        if tid in core2rna_d.keys():
            if strand == '+':
                UTR3[k] = [x for x in trans_exon_loc[core2rna_d[tid]] if x > cds_max]
                UTR5[k] = [x for x in trans_exon_loc[core2rna_d[tid]] if x < cds_min]
                if len(UTR3[k]) %2 == 1:
                    UTR3[k].append(cds_max+10)
                if len(UTR5[k]) %2 == 1:
                    UTR5[k].append(cds_min-10)
            elif strand == '-':
                UTR3[k] = [x for x in trans_exon_loc[core2rna_d[tid]] if x < cds_min]
                UTR5[k] = [x for x in trans_exon_loc[core2rna_d[tid]] if x > cds_max]
                if len(UTR3[k]) %2 == 1:
                    UTR3[k].append(cds_min-10)
                if len(UTR5[k]) %2 == 1:
                    UTR5[k].append(cds_max+10)

    for k,v in UTR3.items():
        v = v.sort()
    for k,v in UTR5.items():
        v = v.sort()

    return UTR3,UTR5


def splitUTR(chrslis,minus_loc,plus_loc,UTR5,UTR3,core_cds_loc):

    splitUTR3,splitUTR5 = {},{}
   

    for chrs in chrslis:
        if chrs in minus_loc.keys():
            if len(minus_loc[chrs]) > 1:
                for mnum in range(len(minus_loc[chrs])-1):
                    lmid, rmid = minus_loc[chrs][mnum], minus_loc[chrs][mnum+1]
                    mabsv = abs(min(core_cds_loc[rmid+' -']) - max(core_cds_loc[lmid+' -']))
                    if mabsv > 10:
                        cut_m_value = int(mabsv*1/3 + max(core_cds_loc[lmid+' -']))
                        if lmid+' -' in UTR5.keys():
                            UTR5[lmid+' -'] = [x for x in UTR5[lmid+' -'] if x < cut_m_value]
                            splitUTR5[lmid] = [x for x in UTR5[lmid+' -'] if x < cut_m_value]
                            if len(splitUTR5[lmid]) %2 == 1:
                                splitUTR5[lmid].append(cut_m_value-15)
                            splitUTR5[lmid].sort()
                            

                        if rmid+' -' in UTR3.keys():
                            UTR3[rmid+' -'] = [x for x in UTR3[rmid+' -'] if x > cut_m_value]
                            splitUTR3[rmid] = [x for x in UTR3[rmid+' -'] if x > cut_m_value]
                            if len(splitUTR3[rmid]) %2 == 1:
                                splitUTR3[rmid].append(cut_m_value+15)
                            splitUTR3[rmid].sort()
                            

            else:
                lrpid = minus_loc[chrs][0]
                splitUTR5[lrpid] = UTR5[lrpid+' -']
                splitUTR3[lrpid] = UTR3[lrpid+' -']
                


        if chrs in plus_loc.keys():
            if len(plus_loc[chrs]) > 1:
                for pnum in range(len(plus_loc[chrs])-1):
                    lpid, rpid = plus_loc[chrs][pnum], plus_loc[chrs][pnum+1]
                    pabsv = abs(min(core_cds_loc[rpid+' +']) - max(core_cds_loc[lpid+' +']))
                    if pabsv > 10:
                        cut_p_value = int(pabsv*2/3 + max(core_cds_loc[lpid+' +']))
                        if lpid+' +' in UTR3.keys():
                            UTR3[lpid+' +'] = [x for x in UTR3[lpid+' +'] if x < cut_p_value]
                            splitUTR3[lpid] = [x for x in UTR3[lpid+' +'] if x < cut_p_value]
                            if len(splitUTR3[lpid]) %2 == 1:
                                splitUTR3[lpid].append(cut_p_value-15)
                            splitUTR3[lpid].sort()
                            

                        if rpid+' +' in UTR5.keys():
                            UTR5[rpid+' +'] = [x for x in UTR3[lpid+' +'] if x > cut_p_value]
                            splitUTR5[rpid] = [x for x in UTR5[rpid+' +'] if x > cut_p_value]
                            if len(splitUTR5[rpid]) %2 == 1:
                                splitUTR5[rpid].append(cut_p_value+15)
                            splitUTR5[rpid].sort()
                            

            else:
                lrpid = plus_loc[chrs][0]
                splitUTR3[lrpid] = UTR3[lrpid+' +']
                splitUTR5[lrpid] = UTR5[lrpid+' +']
                


    return splitUTR3,splitUTR5



def library():
    coregff,core_cds_loc,core_exon_loc,plus_loc,minus_loc,chrs_lis,core2core = core_loc_get(args.core)
    coreUTR3,coreUTR5 = cutUTR(core_cds_loc,core_exon_loc,core2core)
    coreUTR3_sp,coreUTR5_sp = splitUTR(chrs_lis,minus_loc,plus_loc,coreUTR5,coreUTR3,core_cds_loc)

    isoUTR3_sp,isoUTR5_sp,ngsUTR3_sp,ngsUTR5_sp = {},{},{},{}

    if args.isogff is not None and args.ngsgff is not None:
        core2isoID = stat_overlap(args.core,args.isogff)
        iso_exon_loc = transcript_loc_get(args.isogff)
        isoUTR3,isoUTR5 = cutUTR(core_cds_loc,iso_exon_loc,core2isoID)
        isoUTR3_sp,isoUTR5_sp = splitUTR(chrs_lis,minus_loc,plus_loc,isoUTR5,isoUTR3,core_cds_loc)

        core2ngsID = stat_overlap(args.core,args.ngsgff)
        ngs_exon_loc = transcript_loc_get(args.ngsgff)
        ngsUTR3,ngsUTR5 = cutUTR(core_cds_loc,ngs_exon_loc,core2ngsID)
        ngsUTR3_sp,ngsUTR5_sp = splitUTR(chrs_lis,minus_loc,plus_loc,ngsUTR5,ngsUTR3,core_cds_loc)

        return coregff, coreUTR3_sp, coreUTR5_sp, isoUTR3_sp, isoUTR5_sp, ngsUTR3_sp, ngsUTR5_sp

    elif args.isogff is not None and args.ngsgff is None:
        core2isoID = stat_overlap(args.core,args.isogff)
        iso_exon_loc = transcript_loc_get(args.isogff)
        isoUTR3,isoUTR5 = cutUTR(core_cds_loc,iso_exon_loc,core2isoID)
        isoUTR3_sp,isoUTR5_sp = splitUTR(chrs_lis,minus_loc,plus_loc,isoUTR5,isoUTR3,core_cds_loc)

        ngsUTR3_sp,ngsUTR5_sp = {},{}

        return coregff, coreUTR3_sp, coreUTR5_sp, isoUTR3_sp, isoUTR5_sp, ngsUTR3_sp, ngsUTR5_sp

    elif args.ngsgff is not None and args.isogff is None:
        core2ngsID = stat_overlap(args.core,args.ngsgff)
        ngs_exon_loc = transcript_loc_get(args.ngsgff)
        ngsUTR3,ngsUTR5 = cutUTR(core_cds_loc,ngs_exon_loc,core2ngsID)
        ngsUTR3_sp,ngsUTR5_sp  = splitUTR(chrs_lis,minus_loc,plus_loc,ngsUTR5,ngsUTR3,core_cds_loc)

        isoUTR3_sp,isoUTR5_sp = {},{}

        return coregff, coreUTR3_sp, coreUTR5_sp, isoUTR3_sp, isoUTR5_sp, ngsUTR3_sp, ngsUTR5_sp
    elif args.isogff is None and args.ngsgff is None:

        return coregff, coreUTR3_sp, coreUTR5_sp, isoUTR3_sp, isoUTR5_sp, ngsUTR3_sp, ngsUTR5_sp

def outUTR(line,tid,isoUTR,ngsUTR,coreUTR,UTR):
    if tid in isoUTR.keys() and len(isoUTR[tid]) > 0:
        for start,end in zip(isoUTR[tid][:-1:2],isoUTR[tid][1::2]):
            utr = line.copy()
            utr[2] = UTR
            utr[3] = start
            utr[4] = end
            utr[8] = 'Parent='+tid+';'
            print('\t'.join(str(x) for x in utr))
            continue
    else:
        if tid in ngsUTR.keys() and len(ngsUTR[tid]) > 0:
            for start,end in zip(ngsUTR[tid][:-1:2],ngsUTR[tid][1::2]):
                utr = line.copy()
                utr[2] = UTR
                utr[3] = start
                utr[4] = end
                utr[8] = 'Parent='+tid+';'
                print('\t'.join(str(x) for x in utr))
                continue
        else:
            if tid in coreUTR.keys() and len(coreUTR[tid]) > 0:
                for start,end in zip(coreUTR[tid][:-1:2],coreUTR[tid][1::2]):
                    utr = line.copy()
                    utr[2] = UTR
                    utr[3] = start
                    utr[4] = end
                    utr[8] = 'Parent='+tid+';'
                    print('\t'.join(str(x) for x in utr))


def main():
    coregff, coreUTR3, coreUTR5, isoUTR3, isoUTR5, ngsUTR3, ngsUTR5 = library()

    for k,v in coregff.items():
        for tid,fea in v.items():
            for line in fea:
                if line[2] == 'mRNA':
                    print('\t'.join(x for x in line))
                    outUTR(line,tid,isoUTR3,ngsUTR3,coreUTR3,'UTR3')
                    outUTR(line,tid,isoUTR5,ngsUTR5,coreUTR5,'UTR5')
                else:
                    print('\t'.join(x for x in line))


if __name__ == '__main__':
    main()


