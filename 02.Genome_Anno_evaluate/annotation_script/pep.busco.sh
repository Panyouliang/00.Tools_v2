echo start at `date`

export PATH="/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/bin:$PATH"
source /ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/panyouliang/miniconda3/etc/profile.d/conda.sh

conda activate busco


pep=$1

vertebrata_odb10='/jdfssz1/ST_EARTH/P18Z10200N0107/guoqunfei/02.genome_assembly/01.gammaridea/01.Hirondellea_gigas/08.Genome_evaluation/01.integrity/01.Hifiasm_v5/busco_set/vertebrata_odb10'
metazoa_odb10='/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work/02.E.granulosus_ST/10.busco/metazoa_odb10'
sauropsida_odb10='/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/guoqunfei/database/BUSCO/sauropsida_odb10'
arthropoda_odb10='/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/guoqunfei/database/BUSCO/arthropoda_odb10'
hymenoptera_odb10='/ldfssz1/ST_EARTH/Reference/ST_DIVERSITY/USER/guoqunfei/database/BUSCO/hymenoptera_odb10'
aves_odb10='/ldfssz1/ST_EARTH/P18Z10200N0107/panyouliang/my_work5_B10k/buscodb/aves_odb10'



busco -i ${pep} -l $sauropsida_odb10 -o ${pep}_sauropsida -m proteins -f --offline -c 20

echo end at `date`
