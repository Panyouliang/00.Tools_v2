#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,os,gzip,argparse
import subprocess
import glob
from multiprocessing import Pool
import time

#======================================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This pipeline is used for modifing overlaped locus.
Version: v1.0
Author: Panyouliang, panyouliang@genomics.cn
Date: 2023-08-10, yyyy-mm-dd
======================================================================''')
parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [fa,fasta,fa.gz,fasta.gz].')
parser.add_argument('-gff', metavar='gff', type=str, required=True, help='Please input the core files [.gff].')
parser.add_argument('-bam_dir', metavar='bamfile path', type=str, required=False, help='Please input the bamfile path, (OPTIONAL).')
parser.add_argument('-group_id', metavar='groups', type=str,default='P18Z10200N0107',required=False, help='Please assign the groups, (OPTIONAL).')

args = parser.parse_args()
#======================================================================================================

script='/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script'


def run_command(command_list):
    try:
        result = subprocess.run(command_list, stdout=subprocess.PIPE, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def check_job_status(job_id):
    cmd = f"qstat -f | grep {job_id}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    output = result.stdout.strip()
    if job_id not in output:
        return True
    print(f"The job {job_id} is runing...")
    return False

def monitor_job_completed(job_id, interval=10):
    while True:
        if check_job_status(job_id):
            print(f"The job {job_id} has completed!")
            break
        time.sleep(interval)


try:
    gffread = subprocess.run(["which", 'gffread'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    run_gffread = gffread.stdout.strip()
except subprocess.CalledProcessError:
    print("error: No found the software 'gffread', please install the 'gffread' and add it to .bashrc or .bash_profile")
    sys.exit()

try:
    featureCount = subprocess.run(["which", 'featureCounts'], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    run_featureCount = featureCount.stdout.strip()
except subprocess.CalledProcessError:
    print("error: No found the software 'featureCounts', please install the 'featureCounts' and add it to .bashrc or .bash_profile")
    sys.exit()


def filter_overlap_transposon():
    subprocess.run(['perl', script+'/FindOverlapAtCDSlevel.pl', args.gff, args.gff],  stdout=open('overlap.lis', "w"))
    command1 = f"awk '$1!=$2' overlap.lis | awk '$6>0.7 || $7 > 0.7' | awk '{{print $1\"\t\"$2}}' > overlap.lis.id"
    run_command(["bash", '-c', command1])
    command2 = f"awk '{{print $1}}' overlap.lis.id | sort -u > geneID"
    run_command(["bash", '-c', command2])
    subprocess.run(['python', script+'/gff_to_protein.py', '-gff', args.gff, '-ref', args.ref, '-sp', 'output', '-IDtype', 'T'])
    subprocess.run(['/share/app/seqkit/0.14.0-dev/seqkit', 'grep', '-f', 'geneID', 'output.protein.fa'], stdout=open('geneID.fa','w'))
    subprocess.run(['/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin/blastp', '-db', '/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work2/03.align_evolution/00.swissprot/transposon.db', '-query', 'geneID.fa', '-out', 'geneID.fa.m8', '-num_threads', '52', '-qcov_hsp_perc', '60', '-outfmt', '6', '-evalue', '1e-5', '-num_alignments', '1'])
    command3 = f"awk '{{print $1}}' geneID.fa.m8 >overlap_filter.id"
    run_command(['bash', '-c', command3])

def FeatureCount():
    bamdir=glob.glob(args.bam_dir+'/*.bam')
    subprocess.run(['python',script+'/gff2gtf_debug.py',args.gff],stdout=open(args.gff+'.gtf', 'w'))
    featurecount = [run_featureCount, '-T', '16', '-g', 'transcript_id', '-a', args.gff+'.gtf', '-t', 'CDS', '-p', '-o', 'final_featureCounts.txt']
    featurecount += bamdir
    run_feature = open('FeatureCounts.sh','w')
    run_feature.write(' '.join(x for x in featurecount)+'\n')
    run_feature.close()
    os.chmod('FeatureCounts.sh', 0o775)

    qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P',args.group_id,'-l','vf=4g,p=16','-binding','linear:16','FeatureCounts.sh'],stdout=subprocess.PIPE,stderr=subprocess.PIPE, check=True)

    monitor_job_completed(qsubs.stdout.strip().split()[2].decode())

    subprocess.run(['rm','-f',args.gff+'.gtf'])
    subprocess.run(['python', script+'/get_matrix.py', args.gff,'final_featureCounts.txt', 'overlap.lis.id'], stdout=open('overlap_filter.id','a'))

def selectmaxlength():
    subprocess.run(['python', script+'/filter_overlap_gene.py', 'overlap_filter.id', args.gff],stdout=open(args.gff+'.filter.gff','w'))
    subprocess.run(['python', script+'/get_top_length_locus.py', 'overlap.lis.id', args.gff+'.filter.gff'], stdout=open(args.gff+'.filter.gff.selectmax.gff', 'w'))

    subprocess.run(['rm', args.gff+'.filter.gff'])


def trans2gtf():
    subprocess.run(['python', script+'/gff2gtf_debug.py', args.gff+'.filter.gff.selectmax.gff'], stdout=open(args.gff+'.filter.gff.selectmax.gtf','w'))

def check_prefix_file_exists(directory, prefix):
    files = os.listdir(directory)
    for file in files:
        if file.startswith(prefix):
            return True
    return False


def main():

    filter_overlap_transposon()

    if args.bam_dir:
        FeatureCount()

    selectmaxlength()
    trans2gtf()

    dirs=os.getcwd()
    for prefix in ['geneID', 'output', 'overlap']:
        rmfile = f"{prefix}*"
        if check_prefix_file_exists(dirs, prefix):
            rm_command = f"rm -f {rmfile}"
            subprocess.run(rm_command, shell=True)
        else:
            print(f"No files with prefix '{rmfile}' found in the directory.")



if __name__ == '__main__':
    main()
