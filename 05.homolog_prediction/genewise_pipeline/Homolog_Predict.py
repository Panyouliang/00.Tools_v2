#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys,os
import argparse
import subprocess
import time


#======================================================================================================
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This pipeline is used for search homologous from genome.

Version: v1.0
Author: Panyouliang, panyouliang@genomics.cn
Date: 2023-09-26, yyyy-mm-dd
======================================================================''')

parser.add_argument('-v', '--version', action='version', version=version)
parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [fa,fasta,fa.gz,fasta.gz]')
parser.add_argument('-pep', metavar='fasta', type=str, required=True, help='Please input the query files [.fasta,.fa,.fasta.gz,.fa.gz]')
parser.add_argument('-output', metavar='str', type=str, required=True, default='predict', help='Please set the output dir name')
parser.add_argument('-query_split', metavar='int', type=str, required=True, default=2, help='split numbers for query sequence [2]')
parser.add_argument('-alignRate', metavar='float', type=float, required=False, default=0, help='blast alignRatio threshold [0-1]')
parser.add_argument('-identity', metavar='int', type=int, required=False, default=0, help='blast alignRatio threshold [0-100]')
parser.add_argument('-gwis_split', metavar='int', type=str, required=False, default='5', help='split work to run, default=5')
parser.add_argument('-groupID', metavar='str', type=str, required=False, default='P18Z10200N0107', help='groupID for qsub tasks')


args = parser.parse_args()
#=======================================================================================================

scripts='/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/05.homolog_prediction/genewise_pipeline/script'
makeblastdb='/share/app/blast/2.11.0/bin/makeblastdb'
tblastn='/share/app/blast/2.11.0/bin/tblastn'



def get_all_qstat_job():
    qstatlis = []
    cmd = "qstat | tail -n +3 | awk '{print $1}'"
    result = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        for line in result.stdout.splitlines():
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
        time.sleep(10)
        qstatlis = get_all_qstat_job()
        overlap = set(qstatlis) & set(qsublis)


def callblast():
    tblastn = []
    qsublis = []
    path=os.getcwd()+'/'+args.output+'/split'
    current_dir = os.getcwd()
    files = os.listdir(path)
    for f in files:
        file_path = os.path.abspath(os.path.join(path, f))
        if os.path.exists(file_path):
            blast = ['/share/app/blast/2.11.0/bin/tblastn','-query', file_path,'-db', current_dir+'/'+args.ref, '-evalue', '1e-5', '-out', f+'.m8', '-outfmt', '6', '-num_threads', '4']
            run_blast = open(args.output+'/blast/'+f+'.blast.sh','w')
            run_blast.write(' '.join(str(x) for x in blast)+'\n')
            run_blast.close()
            os.chmod(args.output+'/blast/'+f+'.blast.sh',0o775)
            os.chdir(args.output+'/blast')
            qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=1g,p=1', '-binding', 'linear:1', f+'.blast.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            qsublis.append(qsubs.stdout.strip().split()[2].decode())
            os.chdir(current_dir)

    print('All tblastn jobs was qsub!')
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)


def combinehit():
    subprocess.run(['python', scripts+'/cat_file.py', '-p', args.output+'/blast/', '-s','.m8', '-n', args.output+'/combine'])
    subprocess.run(['perl','/home/guoqunfei/bin/pipeline/gene_finding/homology-predict/protein-map-genome/bin/solar/solar.pl','-a','prot2genome2','-f','m8',args.output+'/combine.m8'],stdout=open(args.output+'/combine.m8.solar','w'))
    subprocess.run(['perl',scripts+'/solar_add_realLen.pl',args.output+'/combine.m8.solar', args.pep],stdout=open(args.output+'/combine.m8.solar.cor','w'))
    subprocess.run(['perl',scripts+'/solar_add_identity.pl','--solar',args.output+'/combine.m8.solar.cor','--m8',args.output+'/combine.m8'],stdout=open(args.output+'/combine.m8.solar.cor.idAdd', 'w'))
    awk_run = f"awk '$5>{args.alignRate} && $13>{args.identity}' "+args.output+'/combine.m8.solar.cor.idAdd'
    subprocess.run(awk_run, shell=True, text=True, stdout=open(args.output+'/combine.m8.solar.cor.idAdd.filter', 'w'))

def callgwis():
    os.system('mkdir '+args.output+'/gwis')
    subprocess.run(['perl',scripts+'/get_pos.pl',args.output+'/combine.m8.solar.cor.idAdd.filter'],stdout=open(args.output+'/gwis/gwis.pos','w'))
    subprocess.run(['perl',scripts+'/extract_sequence.pl','--pos',args.output+'/gwis/gwis.pos','--fasta', args.ref,'--extent', '2000'], stdout=open(args.output+'/gwis/gwis.nuc','w'))
    subprocess.run(['perl',scripts+'/prepare_pep.pl',args.output+'/gwis/gwis.pos',args.pep],stdout=open(args.output+'/gwis/gwis.pep','w'))

    awk_run = "awk '{print $1\"\t\"$3}' "+args.output+'/gwis/gwis.pos'
    result = subprocess.run(awk_run, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(args.output+'/gwis/gwis.strandList', 'w') as f:
        f.write(result.stdout)

    qsub_run = subprocess.run(['perl', scripts+'/callgenewise.pl','--pep', args.output+'/gwis/gwis.pep', '--nuc', args.output+'/gwis/gwis.nuc', '--list', args.output+'/gwis/gwis.strandList', '--key', 'gwis', '--out', args.output+'/gwis/', '--num', args.gwis_split], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    print('All genewise jobs was qsub!')

    qsublis = []
    qsublis.append(qsub_run.stdout.strip().split()[2].decode())
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)



def mergegwis():
    subprocess.run(['python', scripts+'/cat_file.py', '-p', args.output+'/gwis/result', '-s', '.gw', '-n', args.output+'/gwis/combine'])
    subprocess.run(['perl',scripts+'/gw_parser.pl', '--gw', args.output+'/gwis/combine.gw', '--pep', args.output+'/gwis/gwis.pep', '--ac', '0', '--id', '0', '--type', '1'], stdout=open(args.output+'/gwis/combine.gw.alg', 'w'))
    subprocess.run(['perl',scripts+'/gw_parser.pl', '--gw', args.output+'/gwis/combine.gw', '--pep', args.output+'/gwis/gwis.pep', '--ac', '0', '--id', '0', '--type', '2'], stdout=open(args.output+'/gwis/combine.gw.mut', 'w'))
    subprocess.run(['perl', scripts+'/merge_overlap.pl', args.output+'/gwis/combine.gw.alg', '0.1'], stdout=open(args.output+'/gwis/combine.gw.alg.nr', 'w'))

    awk_run = "awk '$9>0' "+args.output+'/gwis/combine.gw.alg.nr'
    subprocess.run(awk_run, shell=True, text=True, stdout=open(args.output+'/gwis/combine.gw.alg.nr.filter', 'w'))
    if os.path.exists(args.output+'/output'):
        print('the output dir existing!')
    else:
        os.system('mkdir '+args.output+'/output')

    subprocess.run(['perl', scripts+'/select_gw.pl', args.output+'/gwis/result/', args.output+'/gwis/combine.gw.alg.nr.filter'], stdout=open(args.output+'/output/result.gw', 'w'))
    subprocess.run(['perl', scripts+'/gw_parser.pl', '--gw', args.output+'/output/result.gw', '--pep', args.output+'/gwis/gwis.pep', '--ac', '0', '--id', '0', '--type', '1'], stdout=open(args.output+'/output/result.gw.alg', 'w'))
    subprocess.run(['perl', scripts+'/gw_parser.pl', '--gw', args.output+'/output/result.gw', '--pep', args.output+'/gwis/gwis.pep', '--ac', '0', '--id', '0', '--type', '2'], stdout=open(args.output+'/output/result.gw.mut', 'w'))


def gwis2gff():
    subprocess.run(['perl', scripts+'/gw_to_gff.pl', args.output+'/output/result.gw', args.output+'/gwis/gwis.pep'], stdout=open(args.output+'/output/result.gw.gff', 'w'))
    subprocess.run(['perl', scripts+'/getGene.pl', args.output+'/output/result.gw.gff', args.ref], stdout=open(args.output+'/output/result.gw.cds', 'w'))
    subprocess.run(['perl', scripts+'/cds2aa.pl', args.output+'/output/result.gw.cds'],stdout=open(args.output+'/output/result.gw.pep', 'w'))
    print('All works was done!')


def prepare_data():
    if os.path.exists(args.ref+'.ndb'):
        print('the blast DataBase existing!')
    else:
        subprocess.run([makeblastdb,'-in',args.ref,'-dbtype', 'nucl','-out', args.ref])
        print('blast DataBase was building!')

    if os.path.exists(args.output):
        print(args.output+' has existing!')
    else:
        os.system('mkdir '+args.output)

    if os.path.exists(args.output+'/split'):
        print('the query file has be split!')
    else:
        subprocess.run(['python',scripts+'/fasta_split_non_dir.py', args.pep, args.query_split, 'QueryCut', args.output+'/split'])
        print('the query was split!')

    if os.path.exists(args.output+'/blast'):
        print('the blast dir has existing!')
    else:
        os.system('mkdir '+args.output+'/blast')


def main():
    prepare_data()
    callblast()
    combinehit()
    callgwis()
    mergegwis()
    gwis2gff()



if __name__ == '__main__':
    main()
