gff=$1

perl /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/FindOverlapAtCDSlevel.pl $gff $gff > overlap.txt

awk '$1!=$2' overlap.txt | awk '$6>0.7 || $7 > 0.7' | awk '{print $1"\t"$2}' > locus.overlap.ids
