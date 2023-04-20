#/bin/bash

base_dir=/home/bsorenson/OMI/arctic_comp/
cd $base_dir

make -f Make_omi_colocate

REDO_ALL=false
EXTRACT=true
REPROCESS=false
INCLUDE_NEW=false   # Turns on the pfile check

#base_dir=$(pwd)
#echo $base_dir

data_dir=${base_dir}comp_data/

cd $data_dir

if $EXTRACT; then
  filelist=(${data_dir}*.gz)
  #echo $filelist

  for ((i=0;i < ${#filelist[@]}; i++)); do
    #colocate_file=${data_dir}colocated_subset_${part2}.hdf5
   
    if ($REDO_ALL) ; then
      echo "Extracting ${filelist[$i]}"
      tar -xvzf ${filelist[$i]}
    else
      tester="${filelist[$i]##*_}"
      file_time="${tester%.tar*}"
      just_date=$(echo $file_time | head -c 8)
      #check_file1=${data_dir}${just_date}/${file_time}/omi_pfile_${file_time}.hdf5
      check_file2=${data_dir}${just_date}/${file_time}/omi_shawn_${file_time}.hdf5
      #check_file=${data_dir}${just_date}/${file_time}/omi_shawn_${file_time}.hdf5
      if !(test -f "$check_file2") ; then
        echo "Extracting ${filelist[$i]}"
        tar -xvzf ${filelist[$i]}
      fi
    fi
    #fi

  done
fi 

dirlist=$(ls -d */)
#echo $dirlist

# Loop over the data date directories
# -----------------------------------
for d1 in ${dirlist} ; do
  #echo "${dirlist[$i]}"
  echo $d1
  #cd $d1
  sublist=$(ls -d $d1*/)

  ##echo $sublist

  for d2 in ${sublist} ; do
    #echo $d2
    #part1=$(dirname "$d2")
    part2=$(basename "$d2")

    # Check if the colocated file exists
    # ----------------------------------
    colocate_file=${data_dir}colocated_subset_${part2}.hdf5
    #echo ${colocate_file}
    if !(test -f "$colocate_file") || $REPROCESS ; then
      echo "$colocate_file does NOT exist. Processing data"

      modis_ch2=$(ls ${d2}modis_ch2*)
      modis_ch7=$(ls ${d2}modis_ch7*)
      nsidc_file=$(ls ${d2}nsidc*)
      ceres_file=$(ls ${d2}ceres*)
      omi_file=$(ls ${d2}omi*)

      echo $modis_ch2 $modis_ch7 $nsidc_file $ceres_file $omi_file

      # Run the OMI colocation executable
      # ---------------------------------
      ${base_dir}omi_comp_exec $modis_ch2 $modis_ch7 $nsidc_file $ceres_file $omi_file

      echo
    else
      echo "$colocate_file exists. Not processing anything"
      echo
    fi

    #echo $part2
    #datalist=$(ls $d2)
    #echo $datalist
  done

  #for ((j=0; j < ${#sublist[@]}; j++))
  #do
  #  cd ${sublist[$j]}
  #  ls *
  #done

done
