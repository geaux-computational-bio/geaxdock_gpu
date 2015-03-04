#!/bin/bash

#complex=$var1
#echo $complex

complex=1a07C1

bin=../src/dock
input_dir=.
parameters_dir=${input_dir}/parameters
complex_dir=${input_dir}/${complex}

pdb_file=${complex_dir}/$(echo ${complex:0:5}).pdb
sdf_file=${complex_dir}/${complex}.sdf
ff_file=${complex_dir}/${complex}-0.8.ff

nor_a_file=${parameters_dir}/08_nor_a
nor_b_file=${parameters_dir}/08_nor_b
para_file=${parameters_dir}/paras

TO_ADD_CENTER=sdf_center.py
python $TO_ADD_CENTER -l $sdf_file -f $ff_file


cmd="\
${bin} \
-id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
-para ${para_file} \
\
-ns 3000 \
-nc 10 \
-floor_temp 0.044f \
-ceiling_temp 0.036f \
-nt 1 \
-t 0.02f \
-r 0.08f \
"


echo ${cmd}
${cmd}




