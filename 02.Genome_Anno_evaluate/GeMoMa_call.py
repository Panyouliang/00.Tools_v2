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


def ERE(bamdir,strand):

    qsublis = []

    ERE_tools = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'ERE', 's='+strand, 'maximumcoverage=500', 'outdir=step01.ERE/']
    bamlis = glob.glob(bamdir+'/*bam')
    for bam in bamlis:
        ERE_tools.append('m='+bam)

    ERErun = open('ERE.sh', 'w')
    ERErun.write(' '.join(x for x in ERE_tools)+'\n')
    ERErun.close()
    os.chmod('ERE.sh',0o755)

    qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=50g,num_proc=4', '-binding', 'linear:4', 'ERE.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)

    while True:
        if os.path.exists('step01.ERE/introns.gff'):
            print('ERE was ran done!')
            break
        else:
            time.sleep(30)



def DenoiseIntrons(maxintrons):

    qsublis = []

    denoised = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'DenoiseIntrons', 'i=step01.ERE/introns.gff', 'coverage_unstranded=step01.ERE/coverage.bedgraph', 'm='+maxintrons, 'outdir=step02.DenoiseIntrons/']

    Denois = open('DenoiseI.sh','w')
    Denois.write(' '.join(x for x in denoised)+'\n')
    Denois.close()
    os.chmod('DenoiseI.sh',0o755)

    qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=20g,num_proc=8', '-binding', 'linear:4', 'DenoiseI.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    qsublis.append(qsubs.stdout.strip().split()[2].decode())

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)




def Extractor(homolis):
    if os.path.exists('step03.Extractor'):
        os.system('rm -rf step03.Extractor')
    else:
        pass

    qsublis = []

    with open(homolis,'r') as f:
        for line in f:
            Extract = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'Extractor', 'p=true', 'c=true', 'r=true', 'd=false', 'f=true', 'sefc=true']
            line = line.strip().split()
            name, gff, fna = line[0], line[1], line[2]
            Extract.extend(['a='+gff, 'g='+fna, 'outdir=step03.Extractor/'+name])
            extrt = open(name+'.extract.sh', 'w')
            extrt.write(' '.join(x for x in Extract)+'\n')
            extrt.write('echo "$?: Extractor has done!"\n')
            extrt.close()
            os.chmod(name+'.extract.sh', 0o755)
            qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=30g,num_proc=2', '-binding', 'linear:2', name+'.extract.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            qsublis.append(qsubs.stdout.strip().split()[2].decode())

    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)


def searchhit(ref):

    if os.path.exists('step04.Searchhit'):
        pass
    else:
        os.system('mkdir step04.Searchhit')
    os.chdir('step04.Searchhit')

    qsublis,works = [],[]
    for name in os.listdir(currentdir+'/step03.Extractor/'):
        os.system('mkdir '+name)
        os.chdir(name)
        works.append(name)
        query = glob.glob(currentdir+'/step03.Extractor/'+name+'/cds-parts.fasta')[0]

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

        #qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=100g,num_proc=4', '-binding', 'linear:4', 'mmseqs.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=100g,num_proc=4', '-P', args.groupID+'_super', '-binding', 'linear:4', '-q', 'st_supermem.q', 'mmseqs.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        qsublis.append(qsubs.stdout.strip().split()[2].decode())
        os.chdir(currentdir+'/step04.Searchhit')

    os.chdir(currentdir)
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis,qsublis)

    for work in works:
        while True:
            if os.path.exists(currentdir+'/step04.Searchhit/'+work+'/search.tabular'):
                print(work+' of mapping was ran done!')
                break
            else:
                time.sleep(30)


def GeneModelMapper(ref,maxintrons):

    mapper = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'GeMoMa', 'm='+maxintrons, 'sort=TRUE', 'outdir=output','s=search.tabular', 't=genome.fa', 'c=cds-parts.fasta','a=assignment.tabular','i=denoised_introns.gff','coverage=UNSTRANDED', 'coverage_unstranded=coverage.bedgraph']

    if os.path.exists('step05.GeneModelMapper'):
        os.system('rm -rf step05.GeneModelMapper')

    os.system('mkdir step05.GeneModelMapper')
    os.chdir('step05.GeneModelMapper')

    qsublis,works = [],[]

    for name in os.listdir(currentdir+'/step03.Extractor/'):
        works.append(name)
        os.system('mkdir '+name)
        os.chdir(name)
        os.symlink(currentdir+'/'+ref, currentdir+'/step05.GeneModelMapper/'+name+'/genome.fa')
        os.symlink(currentdir+'/step04.Searchhit/'+name+'/search.tabular', currentdir+'/step05.GeneModelMapper/'+name+'/search.tabular')
        os.symlink(currentdir+'/step03.Extractor/'+name+'/cds-parts.fasta', currentdir+'/step05.GeneModelMapper/'+name+'/cds-parts.fasta')
        os.symlink(currentdir+'/step03.Extractor/'+name+'/assignment.tabular', currentdir+'/step05.GeneModelMapper/'+name+'/assignment.tabular')
        os.symlink(currentdir+'/step01.ERE/coverage.bedgraph', currentdir+'/step05.GeneModelMapper/'+name+'/coverage.bedgraph')
        os.symlink(currentdir+'/step02.DenoiseIntrons/denoised_introns.gff', currentdir+'/step05.GeneModelMapper/'+name+'/denoised_introns.gff')

        Gmapper = open('GeneModelMapper.sh', 'w')
        Gmapper.write(' '.join(x for x in mapper)+'\n')
        Gmapper.write('echo "$?: GeMoMe pipeline"\n')
        Gmapper.close()
        os.chmod('GeneModelMapper.sh', 0o755)

        #qsubs = subprocess.run(['qsub', '-clear','-cwd','-q','st.q','-P', args.groupID, '-l', 'vf=120g,num_proc=8', '-binding', 'linear:8', 'GeneModelMapper.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        qsubs = subprocess.run(['qsub', '-cwd', '-l', 'vf=150g,num_proc=2', '-P', args.groupID+'_super', '-binding', 'linear:2', '-q', 'st_supermem.q', 'GeneModelMapper.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

        qsublis.append(qsubs.stdout.strip().split()[2].decode())
        os.chdir(currentdir+'/step05.GeneModelMapper')

    os.chdir(currentdir)
    qstatlis = get_all_qstat_job()
    moniter_group_completed(qstatlis, qsublis)

    for work in works:
        while True:
            if os.path.exists(currentdir+'/step05.GeneModelMapper/'+work+'/output/predicted_annotation.gff'):
                print(work+' was ran done!')
                break
            else:
                time.sleep(30)


def GMMfilter():
    Filter = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'GAF', 'outdir=output', 'tf=true']
    if os.path.exists('step06.AnnotationFilter'):
        os.system('rm -rf step06.AnnotationFilter')

    os.system('mkdir step06.AnnotationFilter')

    os.chdir('step06.AnnotationFilter')

    for name, gff in zip(os.listdir(currentdir+'/step05.GeneModelMapper/'), glob.glob(currentdir+'/step05.GeneModelMapper/*/output/predicted_annotation.gff')):
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

    if os.path.exists('step07.AnnotationFinalizer'):
        os.system('rm -rf step07.AnnotationFinalizer')

    os.system('mkdir step07.AnnotationFinalizer')
    os.chdir('step07.AnnotationFinalizer')

    qsublis = []

    os.symlink(currentdir+'/'+ref, currentdir+'/step07.AnnotationFinalizer/genome.fa') 
    os.symlink(currentdir+'/step02.DenoiseIntrons/denoised_introns.gff', currentdir+'/step07.AnnotationFinalizer/denoised_introns.gff')
    os.symlink(currentdir+'/step01.ERE/coverage.bedgraph', currentdir+'/step07.AnnotationFinalizer/coverage.bedgraph')
    os.symlink(currentdir+'/step06.AnnotationFilter/output/filtered_predictions.gff', currentdir+'/step07.AnnotationFinalizer/filtered_predictions.gff')
    os.symlink(currentdir+'/step06.AnnotationFilter/output/filtered_predictions_1.gff', currentdir+'/step07.AnnotationFinalizer/filtered_predictions_1.gff')
    os.symlink(currentdir+'/step06.AnnotationFilter/output/filtered_predictions_2.gff', currentdir+'/step07.AnnotationFinalizer/filtered_predictions_2.gff')
    os.symlink(currentdir+'/step06.AnnotationFilter/output/filtered_predictions_3.gff', currentdir+'/step07.AnnotationFinalizer/filtered_predictions_3.gff')
    os.symlink(currentdir+'/step06.AnnotationFilter/output/filtered_predictions_4.gff', currentdir+'/step07.AnnotationFinalizer/filtered_predictions_4.gff')

    finalizes = [java, '-Xms20G', '-Xmx400G', '-jar', GEMOMA, 'CLI', 'AnnotationFinalizer', 'g=genome.fa', 'u=YES', 'i=denoised_introns.gff', 'r=3', 'c=UNSTRANDED', 'coverage_unstranded=coverage.bedgraph', 'n=false', 'rename=SIMPLE', 'p='+name+'G']

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
    parser.add_argument('-bamdir', metavar='path', type=str, required=True, help='Please input the bamfile path.')
    parser.add_argument('-threads', metavar='int', type=str, default='4', required=False, help='Please input the threads number, default=4, [OPTIONAL].')
    parser.add_argument('-maxintrons', metavar='int', type=str, default='100000',required=False, help='Please input the max intron [int], default=100000, [OPTIONAL].')
    parser.add_argument('-strand', metavar='rna-strandness', type=str, default='FR_UNSTRANDED', required=False, help='Defines whether the reads are stranded, range={FR_UNSTRANDED, FR_FIRST_STRAND, FR_SECOND_STRAND}, default = FR_UNSTRANDED, [OPTIONAL]')
    parser.add_argument('-name', metavar='str', type=str,default='GeMoMa',required=False, help='Please input the specie name (OPTIONAL).')
    parser.add_argument('-groupID', metavar='str', type=str,default='P18Z10200N0107',required=False, help='Please assign the groups, default=P18Z10200N0107, (OPTIONAL).')
    parser.add_argument('-MEM', metavar='int', type=str,default='20',required=False, help='Please assign the Memory size, default=20G, (OPTIONAL).')
    args = parser.parse_args()

    #ERE(args.bamdir,args.strand)
    #DenoiseIntrons(args.maxintrons)
    #Extractor(args.homolog_list)
    #searchhit(args.ref)
    GeneModelMapper(args.ref, args.maxintrons)
    GMMfilter()
    AnnotationFinalizer(args.ref,args.name)

