#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,os
import argparse
import subprocess
import time
#=============================================
version = "v1.0"

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [fa,fasta]')
parser.add_argument('-readlis', metavar='txt', type=str, required=True, help='Please input the reads.path.txt files')
parser.add_argument('-output', metavar='dir', type=str, required=True, default='predict', help='Please set the output dir name')
parser.add_argument('-groupID', metavar='str', type=str, required=False, default='P18Z10200N0107', help='groupID for qsub tasks')

args = parser.parse_args()

#=============================================


def prepare(ref,outdir):
    if os.path.exists(ref+'.sa'):
        print('the reference index existing!')
    else:
        subprocess.run(['/share/app/bwa/0.7.17/bwa', 'index', ref])
    if os.path.exists(outdir):
        print(outdir+' has existing!')
    else:
        os.system('mkdir '+outdir)
    if os.path.exists(outdir+'/align'):
        print('the align has existing')
    else:
        os.system('mkdir '+outdir+'/align')


def get_all_qstat_job():
    qstatlis = []
    cmd = "qstat | tail -n +3 | awk '{print $1}'"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        for line in result.stdout.decode('utf-8').splitlines():
            qstatlis.append(line)
    else:
        print(result.stderr)
    return qstatlis

def moniter_group_completed(qstatlis,qsublis):
    overlap = set(qstatlis) & set(qsublis)
    while True:
        if len(overlap) == 0:
            print('All qsub tasks was ran done!')
            break
        time.sleep(30)
        qstatlis = get_all_qstat_job()
        overlap = set(qstatlis) & set(qsublis)



def call_bwa(ref,outdir):

    qsublis = []
    current_dir = os.getcwd()
    with open(args.readlis,'r') as f:
        for line in f:
            line = line.strip().split()
            name,fq1,fq2 = line[0],line[1],line[2]
            bwa = ['/share/app/bwa/0.7.17/bwa', 'mem', '-t', '8', current_dir+'/'+ref, fq1, fq2, '|', 'samtools', 'view', '-@', '8', '-F', '0x4', '-b', '-', '|', 'samtools', 'fixmate', '-m', '-@', '8', '- -', '|', 'samtools', 'sort', '-@', '8', '-o', name+'.sort.bam', '-']
            run_align = open(outdir+'/align/'+name+'.bwa.sh','w')
            run_align.write(' '.join(str(x) for x in bwa)+'\n')
            run_align.write('echo "$?: bwa align task state!"\n')
            run_align.close()
            os.chmod(outdir+'/align/'+name+'.bwa.sh',0o775)
            os.chdir(outdir+'/align')
            qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=2g,p=2', '-binding', 'linear:2', name+'.bwa.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            qsublis.append(qsubs.stdout.strip().split()[2].decode())
            os.chdir(current_dir)

    print('All bwa align jobs was qsub!')
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)


def mergebam(outdir):

    mergelis = []
    qsublis = []
    path = os.getcwd()+'/'+outdir+'/align'
    current_dir = os.getcwd()
    files = os.listdir(path)
    for f in files:
        if f[-3:] == 'bam':
            file_path = os.path.abspath(os.path.join(path,f))
            if os.path.exists(file_path):
                mergelis.append(file_path)

    samtools_merge = open(outdir+'/bam.merge.sh','w')
    samtools_merge.write('samtools merge -@ 16 -o sgs.sort.bam '+' '.join(x for x in mergelis)+'\n')
    samtools_merge.write('samtools markdup -@ 16 -r sgs.sort.bam sgs.sort.markdup.bam\n')
    samtools_merge.write('rm sgs.sort.bam\n')
    samtools_merge.write('samtools index -@ 16 sgs.sort.markdup.bam\n')
    samtools_merge.close()

    os.chmod(outdir+'/bam.merge.sh',0o775)
    os.chdir(outdir)
    #qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=20g,p=4', '-binding', 'linear:4', 'bam.merge.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=40g,num_proc=2', '-P', args.groupID+'_super', '-binding', 'linear:2', '-q', 'st_supermem.q', 'bam.merge.sh'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    qsublis.append(qsubs.stdout.strip().split()[2].decode())
    os.chdir(current_dir)

    print('The samtools merge jobs was qsub!')
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)


def run_polish(ref,outdir,seed):
    qsublis = []
    current_dir = os.getcwd()
    polish = ['python','/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/guoqunfei/src/02.Assessment_genome/NextPolish/lib/nextpolish1.py', '-g', current_dir+'/'+ref, '-t', seed, '-p', '24', '-s', current_dir+'/'+outdir+'/sgs.sort.markdup.bam', '>', 'polish_genome.fa']

    if os.path.exists(outdir+'/result'):
        pass
    else:
        os.system('mkdir '+outdir+'/result')

    polish_run = open(outdir+'/result/polish_run.sh','w')
    polish_run.write(' '.join(str(x) for x in polish)+'\n')
    polish_run.close()

    os.chmod(outdir+'/result/polish_run.sh', 0o775)
    os.chdir(outdir+'/result')
    qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=100g,num_proc=2', '-P', args.groupID+'_super', '-binding', 'linear:2', '-q', 'st_supermem.q', 'polish_run.sh'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    qsublis.append(qsubs.stdout.strip().split()[2].decode())
    os.chdir(current_dir)
    print('The polish jobs was qsub!')

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)
    print('All works was done!!!')

if __name__ == '__main__':

    current_dir = os.getcwd()

    #round 1
    ref,outdir = args.ref,'round1'
    prepare(ref,outdir)
    call_bwa(ref,outdir)
    mergebam(outdir)
    run_polish(ref,outdir,'1')

    #round 2
    ref2,outdir2 = os.path.join(outdir, 'result', 'polish_genome.fa'),'round2'
    prepare(ref2,outdir2)
    call_bwa(ref2,outdir2)
    mergebam(outdir2)
    run_polish(ref2,outdir2,'2')


