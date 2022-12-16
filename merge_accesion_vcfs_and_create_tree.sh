#!/bin/bash

########################
# Initialize variables #
########################
merge_input=""
dir_name=""
vcf_output=$1
tree_output=$2
# echo "vcf_output: $vcf_output";
# echo "tree_output: $tree_output";

########################
# Files pre-processing #
########################
i=0
for vcf_file in "$@"
do
  ((i++))
  if [[ $i -lt 3 ]]; then
    continue
  fi
  echo "Started pre-processing $vcf_file"

	if [ -z "$dir_name" ];
	then
		dir_name=$(dirname ${vcf_file})
	fi

  temp_file=$(basename ${vcf_file})
	temp_file=${temp_file/#//temp_}
	temp_file="$dir_name$temp_file"
  cp ${vcf_file} ${temp_file}
  # echo "Created $temp_file"

  if [[ $temp_file == *.gz ]];
  then
    gzip -d ${temp_file}
    temp_file=${temp_file%".gz"}
		# echo "Decompressed $temp_file"
  fi

  if [[ $temp_file != *.vcf ]];
  then
    printf "$temp_file is not a VCF file" >&2
		rm temp_*
        exit 1
  fi

  bgzip ${temp_file}
  cmprsd_temp_file=${temp_file/%/.gz}
  # echo "Compressed $temp_file into $cmprsd_temp_file"

  tabix -f -p vcf ${cmprsd_temp_file}
  # echo "Indexed $cmprsd_temp_file"

  merge_input="$merge_input $cmprsd_temp_file"
	echo "Finished pre-processing $vcf_file"
done

###############
# VCF merging #
###############
echo "Starting merging..."
bcftools merge --merge none ${merge_input} > ${vcf_output}
echo "Merged all vcf files into $vcf_output"
rm ${dir_name}/temp_*
# echo "Removed all temp files"

###################
# Tree generation #
###################
echo "Starting tree generation..."
vk phylo tree upgma ${vcf_output} > ${tree_output}
echo "Created tree $tree_output"
