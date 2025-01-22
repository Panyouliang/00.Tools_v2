bam=$1

export LD_LIBRARY_PATH="/share/app/gcc-5.2.0/lib64:$LD_LIBRARY_PATH"

date
SCRIPTDIR="/zfssz3/NASCT_BACKUP/RD_09A/GROUP2/ifs1/97.workflow/BWA-ReSeq"
$SCRIPTDIR/genomeCoverageBed -ibam $bam >$bam.coverage
echo "exit code: $?"
samtools view -q 20 -Sb $bam | $SCRIPTDIR/genomeCoverageBed -ibam - >$bam.mapQ20.coverage
echo "exit code: $?"

date
awk '{if(a[$1]!=1){print ;a[$1]=1}}' $bam.coverage|awk '{if($2!=0){print $1"\t"1}else{print $1"\t"(1-$5)}}' >$bam.coverage.per_seq
echo "exit code: $?"
date
awk 'BEGIN {pc=""} 
{
    c=$1;
    if (c == pc) {
        cov=cov+$2*$5;
    } else {
    print pc,cov;
    cov=$2*$5;
    pc=c}
} END {print pc,cov}' $bam.coverage | tail -n +2 > $bam.coverage.percontig
echo "exit code: $?"

awk '(ARGIND==1){len[$1]=$4}(ARGIND==2){cov[$1]=$2}(ARGIND==3){print $1"\t"len[$1]"\t"$2"\t"cov[$1]}' $bam.coverage $bam.coverage.percontig $bam.coverage.per_seq  >$bam.len.cov.depth
echo "exit code: $?"

rm -rf $bam.coverage.per_seq $bam.coverage.percontig

date

date
awk '{if(a[$1]!=1){print ;a[$1]=1}}' $bam.mapQ20.coverage|awk '{if($2!=0){print $1"\t"1}else{print $1"\t"(1-$5)}}' >$bam.mapQ20.coverage.per_seq
echo "exit code: $?"
date
awk 'BEGIN {pc=""} 
{
    c=$1;
    if (c == pc) {
        cov=cov+$2*$5;
    } else {
    print pc,cov;
    cov=$2*$5;
    pc=c}
} END {print pc,cov}' $bam.mapQ20.coverage | tail -n +2 > $bam.mapQ20.coverage.percontig
echo "exit code: $?"

awk '(ARGIND==1){len[$1]=$4}(ARGIND==2){cov[$1]=$2}(ARGIND==3){print $1"\t"len[$1]"\t"$2"\t"cov[$1]}' $bam.mapQ20.coverage $bam.mapQ20.coverage.percontig $bam.mapQ20.coverage.per_seq  >$bam.mapQ20.len.cov.depth
echo "exit code: $?"

rm -rf $bam.mapQ20.coverage.per_seq $bam.mapQ20.coverage.percontig

date
