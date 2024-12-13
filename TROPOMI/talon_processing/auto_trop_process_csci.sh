#/bin/bash

code_dir=/home/blake.sorenson/OMI/tropomi_colocate/
#make -f Make_trop_colocate
make -f Make_trop_colocate_thin

cd $base_dir
base_dir=/home/blake.sorenson/CSCI_544/final_project/

REDO_ALL=false
EXTRACT=true
REPROCESS=false

#base_dir=$(pwd)
#echo $base_dir

data_dir=${base_dir}coloc_data/
prep_dir=${base_dir}prep_data/

cd $prep_dir

if $EXTRACT; then
  filelist=(${prep_dir}*.gz)
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
      check_file=${prep_dir}${file_time}/omi_shawn_${file_time}.hdf5
      #check_file=${prep_dir}${just_date}/${file_time}/omi_shawn_${file_time}.hdf5
      if !(test -f "$check_file") ; then
        echo "Extracting ${filelist[$i]}"
        tar -xvzf ${filelist[$i]}
      fi
    fi
    #fi

  done
fi 

pwd

dirlist=$(ls -d */)
#echo $dirlist

# Loop over the data date directories
# -----------------------------------
for d1 in ${dirlist} ; do
  #echo "${dirlist[$i]}"
  echo $d1
  #cd $d1
  #sublist=$(ls -d $d1*/)

  #part1=$(dirname "$d2")
  part2=$(basename "$d1")


  # Check if the colocated file exists
  # ----------------------------------
  coloc_trop_file=${data_dir}colocated_tropomi_${part2}.hdf5
  #echo ${colocate_file}
  if !(test -f "$coloc_trop_file") || $REPROCESS ; then
    echo "$coloc_trop_file does NOT exist. Processing data"

    trop_file=$(ls ${d1}tropomi_*)
    omi_prep_file=$(ls ${d1}omi_shawn_*)

    # NOTE: ADD A CHECK FOR A MATCHING OMI COMP FILE IN arctic_comp
    #       IF IT EXISTS, RUN TROP_COLOC_THIN
    # -------------------------------------------------------------
    dv1=$(echo $d1 | head -c 8)
    dv2=$(echo $d1 | head -c 12)
    #omi_prep_file=$(ls /home/blake.sorenson/CSCI_/arctic_comp/comp_data/${dv1}/${dv2}/omi_shawn_${dv2}.hdf5)
    #omi_prep_file=$(ls /home/bsorenson/OMI/arctic_comp/comp_data/${dv1}/${dv2}/omi_shawn_${dv2}.hdf5)

    # Run the OMI colocation executable
    # ---------------------------------
    echo "time ${code_dir}trop_coloc_thin_exec $trop_file $omi_prep_file"
    time ${code_dir}trop_coloc_thin_exec $trop_file $omi_prep_file
   
    ##!#if (test -f "$omi_prep_file"); then
    ##!#  echo "OMI prep file found. Processing thin data."

    ##!#  # Run the OMI colocation executable
    ##!#  # ---------------------------------
    ##!#  echo $trop_file $omi_prep_file
    ##!#  time ${base_dir}trop_coloc_thin_exec $trop_file $omi_prep_file

    ##!#else 
    ##!#  echo "No OMI prep file found. Processing full data"

    ##!#  trop_file=$(ls ${d1}tropomi_*)
    ##!#  omi_name_file=$(ls ${d1}omi_filename.txt)

    ##!#  echo $trop_file $omi_name_file

    ##!#  # Run the OMI colocation executable
    ##!#  # ---------------------------------
    ##!#  time ${base_dir}trop_coloc_exec $trop_file $omi_name_file
    ##!#fi 

    echo
  else
    echo "$coloc_trop_file exists. Not processing anything"
    echo
  fi

  #echo $part2
  #datalist=$(ls $d2)
  #echo $datalist

  #for ((j=0; j < ${#sublist[@]}; j++))
  #do
  #  cd ${sublist[$j]}
  #  ls *
  #done

done
