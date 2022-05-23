#!/bin/bash

merge_input=""
dir_name=""
for vcf_file in "$@"
do
	if [ -z "$dir_name" ];
	then
		dir_name=$(dirname ${vcf_file})
	fi
	
    # echo "Started processing $vcf_file"
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
	echo "Finished processing $vcf_file"
done

echo "Starting merging..."
merge_output="$dir_name/all_accessions.vcf"
bcftools merge --merge none ${merge_input} > ${merge_output}
echo "Merged all vcf files into $merge_output"
rm ${dir_name}/temp_*
# echo "Removed all temp files"

echo "Starting tree creation..."
tree_output="$dir_name/tree.newick"
vk phylo tree upgma ${merge_output} > ${tree_output}
echo "Created tree $tree_output"
