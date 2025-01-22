gff=$1
species='IntronSize'

python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/genome_script/The_interval_intron.py $gff $species > intron.txt

python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/genome_script/sort_exon.py intron.txt $species 

