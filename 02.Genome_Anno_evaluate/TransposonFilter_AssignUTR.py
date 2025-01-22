#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,os,gzip,argparse
import subprocess
import importlib

#======================================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This pipeline is used for filtering transposons and annotating UTRs.

Version: v1.0
Author: Panyouliang, panyouliang@genomics.cn
Date: 2023-08-09, yyyy-mm-dd
======================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [fa,fasta,fa.gz,fasta.gz]')
parser.add_argument('-core', metavar='gff', type=str, required=True, help='Please input the core files [.gff]')
parser.add_argument('-isogtf', metavar='gtf', type=str, required=False, help='Please input the ISO-seq assemble files [.gtf], (OPTIONAL)')
parser.add_argument('-ngsgtf', metavar='gtf', type=str, required=False, help='Please input the RNA-seq assemble files [.gtf], (OPTIONAL)')
parser.add_argument('-sp', metavar='specie name', type=str, default='hasky', required=False, help='Please input the specie name, default=hasky, (OPTIONAL)')
parser.add_argument('-rename', metavar='rename', type=str, required=False, help='Please input the abbreviation name of the species, (OPTIONAL)')

args = parser.parse_args()
#======================================================================================================

script_path = '/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script'


try:
    gffread =  subprocess.run(["which", 'gffread'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    run_gffread = gffread.stdout.strip()
except subprocess.CalledProcessError:
    print("error: No found the software 'gffread', please install the 'gffread' and add it to .bashrc or .bash_profile")



def filterTransposon():
    subprocess.run(['python', script_path+'/GetTopIsoform.py', '-gff', args.core, '-sp', args.core])
    subprocess.run(['python', script_path+'/gff_to_protein.py', '-gff', args.core+'.top.iso.coding.gff', '-ref', args.ref, '-sp', args.sp, '-IDtype', 'T'])
    subprocess.run(['perl', script_path+'/fastaDeal.pl', '--attr', 'id:lc', args.sp+'.cds.fa'], stdout=open(args.sp+'.lc', "w"))
    subprocess.run(['awk','($2>=0.7)', args.sp+'.lc'], stdout=open(args.sp+'.lc.0.7', 'w'))
    subprocess.run(['python', script_path+'/drop_ID.py', args.sp+'.lc.0.7', args.core+'.top.iso.coding.gff'], stdout=open(args.core+'.top.iso.coding.gff.filter', 'w'))
    print('Transposon filtering completed!')
    subprocess.run(['python', script_path+'/filter_low_length.cdsgene.py', args.core+'.top.iso.coding.gff.filter'], stdout=open(args.core+'.top.iso.coding.gff.filter.rmrepeat', 'w'))
    print('Low length single exon genes filtering completed!')
    subprocess.run(['rm', args.sp+'.lc', args.sp+'.lc.0.7', args.sp+'.protein.fa', args.sp+'.cds.fa'])


def assignUTR():
    if args.isogtf is not None and args.ngsgtf is None:
        subprocess.run([run_gffread, args.isogtf, '-o', args.isogtf+'.gffread.gff'])
        subprocess.run(['python', script_path+'/new_addUTRs_pipe_v2.py', '-core', args.core+'.top.iso.coding.gff.filter.rmrepeat', '-isogff', args.isogtf+'.gffread.gff'], stdout=open(args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR', 'w'))

        subprocess.run(['rm', args.isogtf+'.gffread.gff'])
        print('Annotation of UTRs from ISO-seq data')

    if args.ngsgtf is not None and args.isogtf is None:

        subprocess.run([run_gffread, args.ngsgtf, '-o', args.ngsgtf+'.gffread.gff'])

        subprocess.run(['python', script_path+'/new_addUTRs_pipe_v2.py', '-core', args.core+'.top.iso.coding.gff.filter.rmrepeat', '-ngsgff', args.ngsgtf+'.gffread.gff'], stdout=open(args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR', 'w'))

        subprocess.run(['rm', args.ngsgtf+'.gffread.gff'])
        print('Annotation of UTRs from mRNA-seq data')

    if args.isogtf and args.ngsgtf:
        subprocess.run([run_gffread, args.isogtf, '-o', args.isogtf+'.gffread.gff'])
        subprocess.run([run_gffread, args.ngsgtf, '-o', args.ngsgtf+'.gffread.gff'])
        subprocess.run(['python', script_path+'/new_addUTRs_pipe_v2.py', '-core', args.core+'.top.iso.coding.gff.filter.rmrepeat', '-isogff', args.isogtf+'.gffread.gff', '-ngsgff', args.ngsgtf+'.gffread.gff'], stdout=open(args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR', 'w'))
        subprocess.run(['rm', args.isogtf+'.gffread.gff', args.ngsgtf+'.gffread.gff'])
        print('Annotation of UTRs from ISO-seq and mRNA-seq data')

    if args.isogtf is None and args.ngsgtf is None:
        subprocess.run(['python', script_path+'/new_addUTRs_pipe_v2.py', '-core', args.core+'.top.iso.coding.gff.filter.rmrepeat'],stdout=open(args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR','w'))

    print('UTRs annotation completed!')

def rename_splitUTR():

    if args.rename:
        subprocess.run(['python', script_path+'/Gff_ReName.py', '-gff', args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR', '-name', args.rename], stdout=open('result_'+args.sp+'.gff', 'w'))
    else:
        subprocess.run(['cat', args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR'], stdout=open('result_'+args.sp+'.gff', 'w'))

    #subprocess.run(['python', script_path+'/gff_format.py', args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR.reName'], stdout=open(args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR.reName.format', 'w'))
    #subprocess.run(['perl', '/hwfssz1/ST_EARTH/P18Z10200N0107/P17Z10200N0101_Metazoa_RNA_Editing/liqiye/PC_PA_UN/bin/annotation/personal/UTR/bin/dealUTR.pl', './', './', 'format'])
    #subprocess.run(['sh', args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR.reName.format.sh'])
    #subprocess.run(['python', script_path+'/replace_id.py', args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR.reName', args.core+'.top.iso.coding.gff.filter.rmrepeat.assignUTR.reName.format.pos.olp.mark.gff.pos.olp.mark.gff.mark.gff'], stdout=open('result_'+args.sp+'.gff', 'w'))
    #print('UTRs annotation correction completed!')


    subprocess.run(['python', script_path+'/gff2gtf_debug.py', 'result_'+args.sp+'.gff'], stdout=open('result_'+args.sp+'.gtf', 'w'))
    print('GTF format conversion completed!')



def check_prefix_file_exists(directory, prefix):
    files = os.listdir(directory)
    for file in files:
        if file.startswith(prefix):
            return True
    return False


def main():
    filterTransposon()
    assignUTR()
    rename_splitUTR()

    dirs=os.getcwd()
    prefix = args.core+'.top.iso.coding.gff'
    pattern = f"{prefix}*"
    if check_prefix_file_exists(dirs, prefix):
        rm_command = f"rm -f {pattern}"
        subprocess.run(rm_command,shell=True)
    else:
        print(f"No files with prefix '{prefix}' found in the directory.")

if __name__ == '__main__':
    main()

