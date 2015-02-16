#!/bin/sh

BASEDIR=~/work/astex
CWD=`pwd`

#MODE=("lig")
#SUFFIX=("-0.4.ligands.sdf")

MODE=("lig" "prt")
SUFFIX=("-0.4.ligands.sdf" "-0.4.templates.pdb")



cd ${BASEDIR}


i=0
while [ $i -lt ${#MODE[@]} ]
do
	for fn in *
	do
#		echo -n ${MODE[i]}
#		ls ${BASEDIR}/${fn}/${fn}${SUFFIX[i]}

		 # echo ""
		 # echo "${CWD}/count ${MODE[i]} ${BASEDIR}/${fn}/${fn}${SUFFIX[i]}"
		${CWD}/count ${MODE[i]} ${BASEDIR}/${fn}/${fn}${SUFFIX[i]}
#		returnval=$?
#		printf "%s \t %s \t %d\n" ${fn} ${MODE[i]} ${returnval}
	done

	i=`expr $i + 1`
done






