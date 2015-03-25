#!/bin/bash

#complex=$var1
#echo $complex

complex=1a07C1

bin=../src/trace
input_dir=.
parameters_dir=${input_dir}/parameters
complex_dir=${input_dir}/${complex}

pdb_file=${complex_dir}/$(echo ${complex:0:5}).pdb
sdf_file=${complex_dir}/${complex}.sdf
ff_file=${complex_dir}/${complex}-0.8.ff

opt_file=${parameters_dir}/08ff_opt
nor_a_file=${parameters_dir}/08_nor_a
nor_b_file=${parameters_dir}/08_nor_b
para_file=${parameters_dir}/paras

trace_file=${parameters_dir}/trace.txt
conf_file=${parameters_dir}/${complex}_new.sdf


cmd="\
${bin} \
-id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
-opt ${opt_file} \
-na ${nor_a_file} \
-nb ${nor_b_file} \
-para ${para_file} \
\
-ns 3000 \
-nc 10 \
-floor_temp 0.044f \
-ceiling_temp 0.036f \
-nt 1 \
-t 0.02f \
-r 0.08f \
-tr ${trace_file} \
-lc ${conf_file}\
"


echo ${cmd}
${cmd}




