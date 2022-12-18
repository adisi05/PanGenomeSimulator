# !/bin/bash
#PBS -S /bin/bash
#PBS -j oe
#PBS -r y
#PBS -q itaym
#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH
#PBS -N simulation_05_x60
#PBS -e /groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/simulation_05_x60/logs
#PBS -o /groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/simulation_05_x60/logs
#PBS -l select=ncpus=1:mem=20gb

###############
# Parameters: #
###############
job_dir=/groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/$PBS_JOBNAME
art_path=/groups/itay_mayrose/adisivan/PanGenomeSimulator/Simulator/ART/art_bin_MountRainier/art_illumina
fasta_dir=/groups/itay_mayrose/adisivan/PanGenomeSimulator/Panoramic/simulation_05
fasta_files=$fasta_dir/output_[A-Z]*.fasta
paired=-p
read_length=150
coverage=60
insert_size=250
insert_std=20

########
# Job: #
########

source ~/.bashrc
hostname
conda activate PanGenomeSimulator
export PATH="$CONDA_PREFIX/bin:$PATH"
echo "Job Name: $PBS_JOBNAME"
touch $job_dir/logs/out.txt
if [[ ! -e $job_dir/logs/out.txt ]]; then
    mkdir -p $job_dir/logs
    touch $job_dir/logs/out.txt
fi
date >> $job_dir/logs/out.txt
cd $job_dir

echo "Starting reads simulation" >> $job_dir/logs/out.txt 2>&1
start=`date +%s`
for fasta_file in $fasta_files
do
  echo "Started simulation for $fasta_file"
  fastq_file="${fasta_file%.fasta}_read"
  fastq_file="${fastq_file#$fasta_dir/}"
  $art_path $paired -i $fasta_file -l $read_length -f $coverage -o $fastq_file -m $insert_size -s $insert_std
	echo "Finished simulation for $fasta_file"
done
end=`date +%s`
runtime=$((end-start))
echo "reads simulation took $runtime seconds" >> $job_dir/logs/out.txt 2>&1
