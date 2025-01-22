samtools merge -@ 16 -o sgs.sort.bam /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/12.genome_assemble/03.polish/test/round2/align/part_020.sort.bam /ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/00.Tools/12.genome_assemble/03.polish/test/round2/align/part_019.sort.bam
samtools markdup -@ 16 -r sgs.sort.bam sgs.sort.markdup.bam
rm sgs.sort.bam
samtools index -@ 16 sgs.sort.markdup.bam
