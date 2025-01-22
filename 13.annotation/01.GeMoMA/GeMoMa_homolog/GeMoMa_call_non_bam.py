import sys,glob
import os,argparse
import subprocess
import time

GEMOMA='/hwfssz1/ST_EARTH/P18Z10200N0160/P18Z10200N0102_GAGA/xiongzj/software/GeMoMa/GeMoMa-1.9.jar'
java='/share/app/java/jdk1.8.0_261/bin/java'
mmseqs='/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin/mmseqs'



currentdir = os.getcwd()
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
        time.sleep(30)
        qstatlis = get_all_qstat_job()
        overlap = set(qstatlis) & set(qsublis)

def qsubwork(shell,mem,num):
    qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf='+str(mem)+'g,num_proc='+str(num), '-P', args.groupID+'_super', '-binding', 'linear:'+str(num), '-q', 'st_supermem.q', shell], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)


def ERE(bamdir,strand):

    qsublis = []

    ERE_tools = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'ERE', 's='+strand, 'maximumcoverage=500', 'outdir=ERE/']
    bamlis = glob.glob(bamdir+'/*bam')
    for bam in bamlis:
        ERE_tools.append('m='+bam)

    ERErun = open('ERE.sh', 'w')
    ERErun.write(' '.join(x for x in ERE_tools)+'\n')
    ERErun.close()
    os.chmod('ERE.sh',0o755)

    #qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=100g,num_proc=8', '-binding', 'linear:4', 'ERE.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=100g,num_proc=4', '-P', args.groupID+'_super', '-binding', 'linear:4', '-q', 'st_supermem.q', 'ERE.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)

    while True:
        if os.path.exists('ERE/introns.gff'):
            print('ERE was ran done!')
            break
        else:
            time.sleep(30)



def DenoiseIntrons(maxintrons):

    qsublis = []

    denoised = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'DenoiseIntrons', 'i=ERE/introns.gff', 'coverage_unstranded=ERE/coverage.bedgraph', 'm='+maxintrons, 'outdir=DenoiseIntrons/']

    Denois = open('DenoiseI.sh','w')
    Denois.write(' '.join(x for x in denoised)+'\n')
    Denois.close()
    os.chmod('DenoiseI.sh',0o755)

    qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=40g,num_proc=4', '-binding', 'linear:4', 'DenoiseI.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)




def Extractor(homolis):
    if os.path.exists('Extractor'):
        os.system('rm -rf Extractor')
    else:
        pass

    qsublis,works = [],[]

    with open(homolis,'r') as f:
        for line in f:
            Extract = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'Extractor', 'p=true', 'c=true', 'r=true', 'd=false', 'f=false']
            line = line.strip().split()
            name, gff, fna = line[0], line[1], line[2]
            works.append(name)
            Extract.extend(['a='+gff, 'g='+fna, 'outdir=Extractor/'+name])
            extrt = open(name+'.extract.sh', 'w')
            extrt.write(' '.join(x for x in Extract)+'\n')
            extrt.write('echo "$?: Extractor has done!"\n')
            extrt.close()
            os.chmod(name+'.extract.sh', 0o755)
            qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=20g,num_proc=4', '-binding', 'linear:4', name+'.extract.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            qsublis.append(qsubs.stdout.strip().split()[2].decode())

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)

    for work in works:
        while True:
            if os.path.exists(currentdir+'/Extractor/'+work+'/assignment.tabular'):
                print(work+' extractor run done!')
                break
            else:
                qsubwork(work+'.extract.sh',20,4)
                time.sleep(30)


def searchhit(ref):

    if os.path.exists('Searchhit'):
        pass
    else:
        os.system('mkdir Searchhit')
    os.chdir('Searchhit')

    qsublis,works = [],[]
    for name in os.listdir(currentdir+'/Extractor/'):
        os.system('mkdir '+name)
        os.chdir(name)
        works.append(name)
        query = glob.glob(currentdir+'/Extractor/'+name+'/cds-parts.fasta')[0]

        refdb = [mmseqs, 'createdb', currentdir+'/'+ref, 'ref.DB']
        qrydb = [mmseqs, 'createdb', query,'query.DB']
        createindex = [mmseqs, 'convertalis', 'ref.DB', 'tmp']
        search = [mmseqs, 'search', 'query.DB', 'ref.DB', 'resultDB', 'align', '-a']
        convertalis = [mmseqs, 'convertalis', 'query.DB', 'ref.DB', 'resultDB', 'search.tabular', '--format-output', '"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen"']

        mapped = open('mmseqs.sh','w')
        mapped.write(' '.join(x for x in refdb)+'\n')
        mapped.write('echo "$?: make refDB"\n')
        mapped.write(' '.join(x for x in qrydb)+'\n')
        mapped.write('echo "$?: make queryDB"\n')
        mapped.write(' '.join(x for x in createindex)+'\n')
        mapped.write('echo "$?: create index"\n')
        mapped.write(' '.join(x for x in search)+'\n')
        mapped.write('echo "$?: search"\n')
        mapped.write(' '.join(x for x in convertalis)+'\n')
        mapped.write('echo "$?: convertalis"\n')

        mapped.close()
        os.chmod('mmseqs.sh',0o755)

        #qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=50g,num_proc=8', '-binding', 'linear:8', 'mmseqs.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=100g,num_proc=4', '-P', args.groupID+'_super', '-binding', 'linear:4', '-q', 'st_supermem.q', 'mmseqs.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        qsublis.append(qsubs.stdout.strip().split()[2].decode())
        os.chdir(currentdir+'/Searchhit')

    os.chdir(currentdir)
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)

    for work in works:
        while True:
            if os.path.exists(currentdir+'/Searchhit/'+work+'/search.tabular'):
                print(work+' of mapping was ran done!')
                break
            else:
                time.sleep(30)


def GeneModelMapper(ref,maxintrons):

    mapper = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'GeMoMa', 'm='+maxintrons, 'sort=TRUE', 'outdir=output','s=search.tabular', 't=genome.fa', 'c=cds-parts.fasta','a=assignment.tabular']

    if os.path.exists('GeneModelMapper'):
        os.system('rm -rf GeneModelMapper')

    os.system('mkdir GeneModelMapper')
    os.chdir('GeneModelMapper')

    qsublis,works = [],[]

    for name in os.listdir(currentdir+'/Extractor/'):
        works.append(name)
        os.system('mkdir '+name)
        os.chdir(name)
        os.symlink(currentdir+'/'+ref, currentdir+'/GeneModelMapper/'+name+'/genome.fa')
        os.symlink(currentdir+'/Searchhit/'+name+'/search.tabular', currentdir+'/GeneModelMapper/'+name+'/search.tabular')
        os.symlink(currentdir+'/Extractor/'+name+'/cds-parts.fasta', currentdir+'/GeneModelMapper/'+name+'/cds-parts.fasta')
        os.symlink(currentdir+'/Extractor/'+name+'/assignment.tabular', currentdir+'/GeneModelMapper/'+name+'/assignment.tabular')
        #os.symlink(currentdir+'/ERE/coverage.bedgraph', currentdir+'/GeneModelMapper/'+name+'/coverage.bedgraph')
        #os.symlink(currentdir+'/DenoiseIntrons/denoised_introns.gff', currentdir+'/GeneModelMapper/'+name+'/denoised_introns.gff')

        Gmapper = open('GeneModelMapper.sh', 'w')
        Gmapper.write(' '.join(x for x in mapper)+'\n')
        Gmapper.write('echo "$?: GeMoMe pipeline"\n')
        Gmapper.close()
        os.chmod('GeneModelMapper.sh', 0o755)

        #qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=200g,num_proc=8', '-binding', 'linear:8', 'GeneModelMapper.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=200g,num_proc=4', '-P', args.groupID+'_super', '-binding', 'linear:2', '-q', 'st_supermem.q', 'GeneModelMapper.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        qsublis.append(qsubs.stdout.strip().split()[2].decode())
        os.chdir(currentdir+'/GeneModelMapper')

    os.chdir(currentdir)
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis, qsublis)

    for work in works:
        while True:
            if os.path.exists(currentdir+'/GeneModelMapper/'+work+'/output/predicted_annotation.gff'):
                print(work+' was ran done!')
                break
            else:
                time.sleep(30)


def GMMfilter():
    Filter = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'GAF', 'outdir=output', 'tf=true']
    if os.path.exists('AnnotationFilter'):
        os.system('rm -rf AnnotationFilter')

    os.system('mkdir AnnotationFilter')

    os.chdir('AnnotationFilter')

    for name, gff in zip(os.listdir(currentdir+'/GeneModelMapper/'), glob.glob(currentdir+'/GeneModelMapper/*/output/predicted_annotation.gff')):
        GAF = ['p='+name, 'g='+gff, 'w=1']
        Filter.extend(GAF)

    #subprocess.run(Filter)
    qsublis = []
    clean = open('Combine_Filter.sh','w')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=1.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=2.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=3.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=4.0"'+'\n')
    clean.write(' '.join(x for x in Filter)+' f="score/aa>=5.0"'+'\n')
    clean.close()
    os.chmod('Combine_Filter.sh',0o755)
    qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=100g,num_proc=2', '-P', args.groupID+'_super', '-binding', 'linear:2', '-q', 'st_supermem.q', 'Combine_Filter.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())
    os.chdir(currentdir)

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis, qsublis)


def AnnotationFinalizer(ref,name):

    if os.path.exists('AnnotationFinalizer'):
        os.system('rm -rf AnnotationFinalizer')

    os.system('mkdir AnnotationFinalizer')
    os.chdir('AnnotationFinalizer')

    qsublis = []

    os.symlink(currentdir+'/'+ref, currentdir+'/AnnotationFinalizer/genome.fa') 
    #os.symlink(currentdir+'/DenoiseIntrons/denoised_introns.gff', currentdir+'/AnnotationFinalizer/denoised_introns.gff')
    #os.symlink(currentdir+'/ERE/coverage.bedgraph', currentdir+'/AnnotationFinalizer/coverage.bedgraph')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions.gff', currentdir+'/AnnotationFinalizer/filtered_predictions.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_1.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_1.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_2.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_2.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_3.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_3.gff')
    os.symlink(currentdir+'/AnnotationFilter/output/filtered_predictions_4.gff', currentdir+'/AnnotationFinalizer/filtered_predictions_4.gff')

    finalizes = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'AnnotationFinalizer', 'g=genome.fa', 'n=false', 'rename=SIMPLE', 'p='+name+'G']

    #subprocess.run(finalizes)
    Final = open('Finalizer.sh','w')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions.gff outdir=Result_score2aa_1.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_1.gff outdir=Result_score2aa_2.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_2.gff outdir=Result_score2aa_3.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_3.gff outdir=Result_score2aa_4.0/'+'\n')
    Final.write(' '.join(x for x in finalizes)+' a=filtered_predictions_4.gff outdir=Result_score2aa_5.0/'+'\n')
    Final.write('echo "$?: Finalizer done!"\n')
    Final.close()
    os.chmod('Finalizer.sh',0o755)
    qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=100g,num_proc=2', '-P', args.groupID+'_super', '-binding', 'linear:2', '-q', 'st_supermem.q', 'Finalizer.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis, qsublis)
    os.chdir(currentdir)


if __name__ == '__main__':

    version = "v1.0"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
======================================================================
This pipeline is used for call GeMoMa annotating pipeline.

Version: v1.0
Author: Panyouliang, panyouliang@genomics.cn
Date: 2024-01-06, yyyy-mm-dd
======================================================================''')
    parser.add_argument('-v', '--version', action='version', version=version)
    parser.add_argument('-ref', metavar='fasta', type=str, required=True, help='Please input the reference genome [.fa, .fasta].')
    parser.add_argument('-homolog_list', metavar='file', type=str, required=True, help='Please input the homolog_list.txt.')
    parser.add_argument('-bamdir', metavar='path', type=str, required=False, help='Please input the bamfile path.')
    parser.add_argument('-threads', metavar='int', type=str, default='4', required=False, help='Please input the threads number, default=4, [OPTIONAL].')
    parser.add_argument('-maxintrons', metavar='int', type=str, default='100000',required=False, help='Please input the max intron [int], default=100000, [OPTIONAL].')
    parser.add_argument('-strand', metavar='rna-strandness', type=str, default='FR_UNSTRANDED', required=False, help='Defines whether the reads are stranded, range={FR_UNSTRANDED, FR_FIRST_STRAND, FR_SECOND_STRAND}, default = FR_UNSTRANDED, [OPTIONAL]')
    parser.add_argument('-name', metavar='str', type=str,default='GeMoMa',required=False, help='Please input the specie name (OPTIONAL).')
    parser.add_argument('-groupID', metavar='str', type=str,default='P18Z10200N0107',required=False, help='Please assign the groups, default=P18Z10200N0107, (OPTIONAL).')
    parser.add_argument('-MEM', metavar='int', type=str,default='20',required=False, help='Please assign the Memory size, default=20G, (OPTIONAL).')
    args = parser.parse_args()

    #ERE(args.bamdir,args.strand)
    #DenoiseIntrons(args.maxintrons)
    Extractor(args.homolog_list)
    searchhit(args.ref)
    GeneModelMapper(args.ref, args.maxintrons)
    GMMfilter()
    AnnotationFinalizer(args.ref,args.name)

