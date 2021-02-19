#!/bin/bash
#SBATCH --partition=pi_kaminski
#SBATCH --job-name=my_immcantation_job
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1-00:00:00


start=$(date "+%s")

# define input/output data folder
input_folder="/gpfs/loomis/pi/kaminski/xy48/MISC/lucas_data"
output_folder="/gpfs/loomis/pi/kaminski/xy48/MISC/lucas_data_changeo"

# make output folder
# if folder exist, delete it first
[ -d "${output_folder}" ] && rm -r $output_folder
mkdir $output_folder

# each folder represent a sample; loop all folders
echo $input_folder
for path in $input_folder/*
do
    [ -d "${path}" ] || continue # if not a directory, skip
    dirname="$(basename "${path}")"
    
    # split string to get sample name and type
    # the folder name pattern: sample name-TCR/BCR_VHT_cellranger
    s_name=${dirname:0:(-19)}
    s_type=${dirname:(-18):(-15)}
    
    # only process TCR data
    if [[ "$s_type" == "TCR" ]]
    then
        singularity exec -B ${path}:/data immcantation-latest_dev.sif \
        changeo-10x -s /data/filtered_contig.fasta \
                    -a /data/filtered_contig_annotations.csv \
                    -o ${output_folder}/${s_name} \
                    -g human \
                    -t tr
    fi
done

stop=$(date "+%s")
time=$(echo "scale=2;($stop - $start) / 3600.0" | bc)
echo "Whole pipeline completed. Elapsed time: $time hours."