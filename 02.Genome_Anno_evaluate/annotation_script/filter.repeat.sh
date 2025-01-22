#! /usr/bin/bash

cds=$1
gff=$2

perl /home/xiongzj/bin/fastaDeal.pl --attr id:lc $cds > $cds.lc
awk '($2>=0.7)' $cds.lc > $cds.lc.0.7

python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/drop_ID.py $cds.lc.0.7 $gff  > $gff.filter.gff
#perl /home/xiongzj/bin/filter_gff_gene_lenght.pl --threshold 300 --exons 1 $gff.filter.gff > $gff.filter.gff.1.gff
#perl /home/xiongzj/bin/filter_gff_gene_lenght.pl --threshold 150 --exons 2 $gff.filter.gff.1.gff > $gff.filter.gff.2.gff

python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/filter_low_length.cdsgene.py $gff.filter.gff > $gff.detransposon.gff

rm $cds.lc $cds.lc.0.7 $gff.filter.gff
#python /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/02.Genome_Anno_evaluate/annotation_script/GetTopIsoform.py -gff $gff.filter.gff.relow.gff -sp $gff.filter.gff.relow.gff
