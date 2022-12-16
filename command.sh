# !/bin/bash
#PBS -S /bin/bash
#PBS -j oe
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N simulation_00
#PBS -e /groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/simulation_00/logs
#PBS -o /groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/simulation_00/logs
#PBS -l select=ncpus=4:mem=10gb

###############
# Parameters: #
###############
code_dir=/groups/itay_mayrose/adisivan/PanGenomeSimulator
job_dir=/groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/$PBS_JOBNAME
accessions_vcf_files=$code_dir/vcf_inputs/intersection_*
vcf_file=$code_dir/vcf_inputs/all_accessions.vcf
#vcf_file=$job_dir/all.snp.filter.ref.variable.rename.vcf
tree_file=$code_dir/vcf_inputs/tree.newick
#tree_file=$job_dir/all.snp.filter.ref.variable.rename.min9.phy.treefile
gff_file=/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53_1000.gff3
#gff_file=/groups/itay_mayrose/adisivan/arabidopsis/ensemblgenomes/gff3/Arabidopsis_thaliana.TAIR10.53.gff3
annotations_file=$code_dir/test_arabidopsis/test_all_chroms_annotations.csv
#annotations_file=$job_dir/all_chroms_annotations.csv
fasta_ref=$code_dir/test_arabidopsis/Arabidopsis_3_chroms_1000.f
#fasta_ref=/groups/itay_mayrose/adisivan/PanGenomeSimulator/test_arabidopsis/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
model_file=$job_dir/mut_model.p

########
# Job: #
########

## 0. Initialization
source ~/.bashrc
hostname
conda activate PanGenomeSimulator
export PATH="$CONDA_PREFIX/bin:$PATH"
cd $code_dir
echo "Job Name: $PBS_JOBNAME"
touch $job_dir/logs/out.txt
if [[ ! -e $job_dir/logs/out.txt ]]; then
    mkdir -p $job_dir/logs
    touch $job_dir/logs/out.txt
fi
date >> $job_dir/logs/out.txt

## 1. Process VCF files to infer tree and build a merged VCF file
start=`date +%s`
sh ./merge_accesion_vcfs_and_create_tree.sh $vcf_file $tree_file $accessions_vcf_files
end=`date +%s`
runtime=$((end-start))
echo "VCF merging and tree generation took $runtime seconds" >> $job_dir/logs/out.txt 2>&1

## 2. Process annotation file
python process_annotations.py --gff $gff_file --csv $annotations_file

## 3. Build model
start=`date +%s`
python /groups/itay_mayrose/adisivan/PanGenomeSimulator/mutation_model.py \
-r $fasta_ref -v $vcf_file -o $model_file -a $annotations_file
end=`date +%s`
runtime=$((end-start))
echo "Model building took $runtime seconds" >> $job_dir/logs/out.txt 2>&1

## 4. Simulate pan-genome
start=`date +%s`
python /groups/itay_mayrose/adisivan/PanGenomeSimulator/pangenome_simulator.py \
-r $fasta_ref -m $model_file -o $job_dir/output -a $annotations_file -t $tree_file --vcf --max-threads 3 \
-R 150 --pe 250 20 -c 30
end=`date +%s`
runtime=$((end-start))
echo "Simulation took $runtime seconds" >> $job_dir/logs/out.txt 2>&1
