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

trace_file=./trajectory.sample.out
conf_file=${parameters_dir}/${complex}_new



cmd="\
${bin} \
--id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
--para ${para_file} \
--trace \
--conf $conf_file \
\
--csv $trace_file \
--nc 10 \
--floor_temp 0.04 \
--ceiling_temp 0.036 \
--nt 1 \
-t 0.02 \
-r 0.08 \
"

echo ${cmd}
${cmd}




