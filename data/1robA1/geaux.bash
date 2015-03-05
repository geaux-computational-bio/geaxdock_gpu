#!/bin/bash


################################################################################
# run GeauxDock
################################################################################

# user-defined
complex=1robA1

bin=../../src/dock
# TO_ADD_CENTER=/home/jaydy/Workspace/Bitbucket/geauxdock/data/sdf_center.py

ff_dir=.
parameters_dir=../parameters

pdb=$(echo ${complex:0:5}).pdb
ff=${complex}-0.8.ff
sdf=${complex}.sdf

pdb_file=${ff_dir}/$pdb
sdf_file=${ff_dir}/$sdf
ff_file=${ff_dir}/$ff


para_file=${parameters_dir}/paras

# running

cmd="\
${bin} \
-id ${complex} \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
-para ${para_file} \
\
-nc 10 \
-floor_temp 0.04f \
-ceiling_temp 0.036f \
-nt 1 \
-t 0.02f \
-r 0.08f \
"


echo ${cmd}
${cmd}
