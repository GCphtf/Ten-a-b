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

for srr in $(cat DC3_Acc.txt); do
  fasterq-dump --outdir out_dir ${sra_dir}/${srr}/${srr}.sra
done



