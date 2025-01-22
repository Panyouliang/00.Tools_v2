#! /usr/bin/bash

path=$1
genome=$2

echo "

/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin/BuildDatabase -engine ncbi -name mydb $path/$genome > $genome.repeatmodeler.log

/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin/RepeatModeler -engine ncbi -database mydb -threads 12  > run.out



" > $genome.sh

qsub -clear -cwd -q st.q -P P18Z10200N0107 -l vf=6g,p=12 -binding linear:12 $genome.sh
