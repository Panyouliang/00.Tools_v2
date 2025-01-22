bam=$1
depth=5

export LD_LIBRARY_PATH="/share/app/gcc-5.2.0/lib64:$LD_LIBRARY_PATH"

date
SCRIPTDIR="/zfssz3/NASCT_BACKUP/RD_09A/GROUP2/ifs1/97.workflow/BWA-ReSeq"
$SCRIPTDIR/genomeCoverageBed -ibam $bam >$bam.coverage
echo "exit code: $?"
samtools view -q 20 -Sb $bam | $SCRIPTDIR/genomeCoverageBed -ibam - >$bam.mapQ20.coverage
echo "exit code: $?"

date

python3 /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/03.AlignInfoStatis/depth.pre.py $bam.coverage $depth >$bam.len.cov.depth
echo "exit code: $?"

date

python3 /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/03.AlignInfoStatis/depth.pre.py $bam.mapQ20.coverage $depth >$bam.mapQ20.len.cov.depth
echo "exit code: $?"

date
