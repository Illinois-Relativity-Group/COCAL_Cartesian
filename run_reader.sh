#!/bin/bash

cd $reader_dir

echo "Datatype: $data"
echo "Root directory: $root"
echo "Grid size: $x_min to $x_max, $y_min to $y_max, $z_min to $z_max"
echo "Resolution: $res"

echo "Reading data from $output_dir"

sed -i "2s|.*|$output_dir|" PATH2ID.txt

if [ "$data" = "bht" ]; then
    gfortran coc2cac_bht_WL_cartesian.f90 -o coc2cac_bht_WL_cartesian
    echo "Running coc2cac_bht_WL_cartesian.f90"
    time ./coc2cac_bht_WL_cartesian $x_min $x_max $y_min $y_max $z_min $z_max $res

elif [ "$data" = "rns" ]; then
    gfortran coc2cac_rs_WL_v1_cartesian.f90 -o coc2cac_rs_WL_v1_cartesian
    echo "Running coc2cac_rs_WL_v1_cartesian.f90"
    time ./coc2cac_rs_WL_v1_cartesian $x_min $x_max $y_min $y_max $z_min $z_max $res

elif [ "$data" = "bns" ]; then
    gfortran coc2cac_bns_cartesian.f90 -o coc2cac_bns_cartesian
    echo "Running coc2cac_bns_cartesian.f90"
    time ./coc2cac_bns_cartesian $x_min $x_max $y_min $y_max $z_min $z_max $res
elif [ "$data" = "nst" ]; then
    gfortran coc2cac_nst_cartesian.f90 -o coc2cac_nst_cartesian
    echo "Running coc2cac_nst_cartesian.f90"
    time ./coc2cac_nst_cartesian $x_min $x_max $y_min $y_max $z_min $z_max $res
fi

mv data.txt $root

cd $root

sed -i -E 's/([0-9])([+-][0-9]{2,3})/\1E\2/g' data.txt

tail -n +2 data.txt | awk '{print $1, $2, $3, $4+$5}' > plot.3d

python3 convert_to_hdf5.py

cd $root
