1, Download metadata files and Accession lists from GSE100058 for PN and GSE143038 for ORN

2, download SRR files. It doesn't work on Xshell so download SRR in windows then transfer to Xshell

acc="path_to_acc_list"
ourdir1="path_to_output_directory"
cat $acc | xargs -n 1 -P 4 prefetch -O ${outdir}

3, fasterq-dump. This time inside Xshell

#!/bin/bash
#SBATCH --job-name=Fasterq-dump
#SBATCH --partition=cu
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --output=/public/home/tongchaogroup/yijielin/err/Fasterq-dump.out
#SBATCH --error=/public/home/tongchaogroup/yijielin/err/Fasterq-dump.err

ulimit -n 4096

conda init bash
eval "$(/public/home/tongchaogroup/yijielin/miniconda3/bin/conda shell.bash hook)"
conda activate general

sra_dir="/public/home/tongchaogroup/yijielin/Acc_List/Mz19_srr_files"
out_dir="/public/home/tongchaogroup/yijielin/FastQ/PN_Mz19" 

for srr in $(cat /public/home/tongchaogroup/yijielin/Acc_List/Mz19_Acc_List.txt); do
  fasterq-dump --outdir $out_dir ${sra_dir}/${srr}/${srr}.sra
done

4, STAR.

#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --partition=cu
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --output=/public/home/tongchaogroup/yijielin/err/STAR.out
#SBATCH --error=/public/home/tongchaogroup/yijielin/err/STAR.err

ulimit -n 4096

conda init bash
eval "$(/public/home/tongchaogroup/yijielin/miniconda3/bin/conda shell.bash hook)"
conda activate general

acc_list="/public/home/tongchaogroup/yijielin/Acc_List/Mz19_Acc_List.txt"
fastq_dir="/public/home/tongchaogroup/yijielin/FastQ/PN_Mz19"
output_dir="/public/home/tongchaogroup/yijielin/star_result/PN_Mz19"
genome_dir="/home/lilab/ten_a_isoform_project/genome_directory_dm6_100bp"  
sj_output_dir="/public/home/tongchaogroup/yijielin/star_sj/PN_Mz19"
BAM_output_dir="/public/home/tongchaogroup/yijielin/star_BAM/PN_Mz19"

while read -r srr; do
    fq1="$fastq_dir/${srr}_1.fastq"
    fq2="$fastq_dir/${srr}_2.fastq"

    # Check if both FASTQ files exist
    if [[ -f "$fq1" && -f "$fq2" ]]; then

        # Run STAR for this SRR sample
        STAR --runThreadN 32 \
             --genomeDir "$genome_dir" \
             --readFilesIn "$fq1" "$fq2" \
             --outFileNamePrefix "$output_dir/${srr}_" \
             --outSAMtype BAM SortedByCoordinate
             
        mv "${output_dir}/${srr}_SJ.out.tab" "$sj_output_dir/"
        mv "${output_dir}/${srr}Aligned.sortedByCoord.out.bam" "$BAM_output_dir/"

    else
        echo "Warning: FASTQ files for $srr not found. Skipping."
    fi
done < "$acc_list"













