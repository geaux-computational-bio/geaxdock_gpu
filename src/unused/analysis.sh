################################################################################
# use ./analysis to load the hdf5 record files and redirect the output to a csv format
hd_paths=$(ls output_*/*.h5)
rep=$1

for hd_path in $hd_paths; do
    echo -e "hdf file path:\t\t\t"$hd_path

    dir_path=$(dirname $hd_path)

    hd_fn=${hd_path##*/}

    base_name=${hd_fn%.h5}

    csv_path=$dir_path/${base_name}_$rep.csv
    echo -e "csv output path:\t\t"$csv_path"\n"

    # ./analysis -nl 20 -l 1 -p 2 -r 0 -e 1 $hd_path 
    ./analysis -rep $rep -nl 4500 -l 1 -p 2 -e $hd_path > $csv_path
done


################################################################################
# use pandas to laod the csv file and remove the duplicates
# echo $dir_path
# echo $base_name

python removeDup.py $dir_path ${rep}_${base_name}
