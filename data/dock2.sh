#!/usr/bin/env bash


readonly bin=../src/dock

readonly paras=./parameters/paras

pdb_file=./10gs/10gsA00.pdb
sdf_file=./10gs/10gs_ligand.sdf
ff_file=./10gs/10gsA00.ff

csv_file=./10gs/10gsA00.csv


cmd="\
${bin} \
--id  10gsA00 \
-p ${pdb_file} \
-l ${sdf_file} \
-s ${ff_file} \
\
--para ${paras} \
\
--csv $csv_file \
 --nc 10 \
--floor_temp 0.004 \
--ceiling_temp 0.036 \
--nt 1 \
-t 0.02 \
-r 0.08 \
"


echo ${cmd}
${cmd}
