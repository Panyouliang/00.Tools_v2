#! /usr/bin/bash

genome=$1
bam=$2

echo date;
/share/app/samtools/1.11/bin/samtools depth $bam > $bam.depth.txt

echo date;
python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/03.AlignInfoStatis/depth_coverage.py $genome $bam.depth.txt > $bam.coverage_depth.txt

echo date;
