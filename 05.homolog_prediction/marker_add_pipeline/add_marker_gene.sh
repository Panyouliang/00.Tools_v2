
gwpep=$1
gwgff=$2
gwalg=$3
reblastdb=$4
refgff=$5
markerlis=$6


/hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/tanshangjin/software/miniconda3/envs/mmseqs2/bin/mmseqs easy-search $gwpep $reblastdb blast.m8 tmp --threads 16


python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/05.homolog_prediction/marker_add_pipeline/script/mmseq_get.best-hit.py blast.m8 $markerlis > blast.m8.id2id

perl /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/FindOverlapAtCDSlevel.pl $gwgff $refgff > overlap.CDS.lis

python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/05.homolog_prediction/marker_add_pipeline/script/get_blast2overlap.py blast.m8.id2id $gwalg overlap.CDS.lis > output.txt
