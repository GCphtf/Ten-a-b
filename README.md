1, Download metadata files and Accession lists from GSE100058 for PN and GSE143038 for ORN----------------------------------------------------------------------------------------------------------------------------------------------------------

2, download SRR files. It doesn't work on Xshell so download SRR in windows then transfer to Xshell------------------------------------------------------------------------------------------------------------------------------------------------------------

acc="path_to_acc_list"
ourdir1="path_to_output_directory"
cat $acc | xargs -n 1 -P 4 prefetch -O ${outdir}

3, fasterq-dump. This time inside Xshell----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

4, STAR.----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

5, detect junction-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd

jb1 = 12194224
jb2 = 12195768
jb3 = 12196624
jb4 = 12197185
jb5 = 12044350
jb6 = 12077154
jb7 = 12077242
jb8 = 12145794

ja1 = 12230371
ja2 = 12273934
ja3 = 12203879
ja4 = 12204138
ja5 = 12203985
ja6 = 12204138
ja7 = 12203976
ja8 = 12204138
ja9 = 12204410
ja10 = 12229665
ja11 = 12194224
ja12 = 12195768
ja13 = 12196624
ja14 = 12197185
ja15 = 12197256
ja16 = 12204138
ja17 = 12204410
ja18 = 12229665


# Define paths for input and output
acc_list_path = "/public/home/tongchaogroup/yijielin/Acc_List/Mz19_Acc_List.txt"
outpath = "/public/home/tongchaogroup/yijielin/find_junction_from_sj_file/PN_Mz19_all_junction_b_andor.csv"

# Initialize lists to store SRR IDs and junction counts
SRR_list = []
Junction_count_a = []
Junction_count_b = []
processed = 0

# Open the accession list and loop through each SRR accession
with open(acc_list_path, 'r') as file:
    for line in file:
        srr = line.strip()  # Clean the SRR ID
        SRR_list.append(srr)

        try:
            # Load the SJ.out.tab file for the current SRR
            srr_file = f"/public/home/tongchaogroup/yijielin/star_sj/PN_Mz19/{srr}_SJ.out.tab"
            df = pd.read_csv(srr_file, sep="\t", header=None)

            # Assign column names based on STAR's SJ.out.tab structure
            df.columns = [
                "chromosome", "start", "end", "strand", "intron_motif",
                "annotated", "unique_reads", "multi_reads", "max_overhang"
            ]

            # Filter the DataFrame for the specified junction(s)
            jc_dfb1 = df[(df["start"] == jb1) | (df["end"] == jb2)]
            jc_dfb2 = df[(df["start"] == jb3) | (df["end"] == jb4)]
            jc_dfb3 = df[(df["start"] == jb5) & (df["end"] == jb6)]
            jc_dfb4 = df[(df["start"] == jb7) & (df["end"] == jb8)]
            jc_total_b_df = pd.concat([jc_dfb1,jc_dfb2,jc_dfb3,jc_dfb4], ignore_index=True)

            # For a: Filter the DataFrame for the specified junction(s)
            jc_dfa1 = df[(df["start"] == ja1) | (df["end"] == ja2)]
            jc_dfa2 = df[(df["start"] == ja3) & (df["end"] == ja4)]
            jc_dfa3 = df[(df["start"] == ja5) & (df["end"] == ja6)]
            jc_dfa4 = df[(df["start"] == ja7) & (df["end"] == ja8)]
            jc_dfa5 = df[(df["start"] == ja9) & (df["end"] == ja10)]
            jc_dfa6 = df[(df["start"] == ja11) & (df["end"] == ja12)]
            jc_dfa7 = df[(df["start"] == ja13) & (df["end"] == ja14)]
            jc_dfa8 = df[(df["start"] == ja15) & (df["end"] == ja16)]
            jc_dfa9 = df[(df["start"] == ja17) & (df["end"] == ja18)]
            jc_total_a_df = pd.concat([jc_dfa1,jc_dfa2,jc_dfa3,jc_dfa4,jc_dfa5,jc_dfa6,jc_dfa7,jc_dfa8,jc_dfa9], ignore_index=True)

            # Check if the junction exists and get the unique read count
            if not jc_total_df_a.empty:
                unique_reads_a = jc_total_a_df.iloc[0]["unique_reads"]
                Junction_count_a.append(unique_reads_a)
                print(f"{srr} found with {unique_reads} reads for a.")
            else:
                Junction_count_a.append(0)
                if processed == 0:
                    processed = 1
                    print("code is normal, a file is processed with 0 reads mapped to the junction")
            if not jc_total_df_b.empty:
                unique_reads_b = jc_total_b_df.iloc[0]["unique_reads"]
                Junction_count_b.append(unique_reads_b)
                print(f"{srr} found with {unique_reads} reads for b.")
            else:
                Junction_count_b.append(0)
                if processed == 0:
                    processed = 1
                    print("code is normal, a file is processed with 0 reads mapped to the junction")
        except Exception as e:
            print(f"Error processing {srr}: {e}")
            Junction_count.append(0)  # Assign 0 for problematic files

# Create a DataFrame from the results and save to CSV
jc_df_a = pd.DataFrame({
    'SRR': SRR_list,
    'Junction_count': Junction_count_a
})

jc_df_b = pd.DataFrame({
    'SRR': SRR_list,
    'Junction_count': Junction_count_b
})

final_df = jc_df_a[jc_df_a["Junction_count"] >= 0]

final_df.to_csv(outpath, index=False)
print(f"Results saved to: {outpath}")

5, detect junction, within have a see b-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import pandas as pd

jb1 = 12194224
jb2 = 12195768
jb3 = 12196624
jb4 = 12197185
jb5 = 12044350
jb6 = 12077154
jb7 = 12077242
jb8 = 12145794

ja1 = 12230371
ja2 = 12273934
ja3 = 12203879
ja4 = 12204138
ja5 = 12203985
ja6 = 12204138
ja7 = 12203976
ja8 = 12204138
ja9 = 12204410
ja10 = 12229665
ja11 = 12194224
ja12 = 12195768
ja13 = 12196624
ja14 = 12197185
ja15 = 12197256
ja16 = 12204138
ja17 = 12204410
ja18 = 12229665


# Define paths for input and output
acc_list_path = "/public/home/tongchaogroup/yijielin/Acc_List/Mz19_Acc_List.txt"
outpath1 = "/public/home/tongchaogroup/yijielin/find_junction_from_sj_file/PN_Mz19_all_jc_have_a_see_b_xiaoyu10.csv"
outpath2 = "/public/home/tongchaogroup/yijielin/find_junction_from_sj_file/PN_Mz19_all_jc_have_a_see_b_dayu10.csv"
# Initialize lists to store SRR IDs and junction counts
SRR_list = []
Junction_count_b = []
processed = 0

# Open the accession list and loop through each SRR accession
with open(acc_list_path, 'r') as file:
    for line in file:
        srr = line.strip()  # Clean the SRR ID
        SRR_list.append(srr)

        try:
            # Load the SJ.out.tab file for the current SRR
            srr_file = f"/public/home/tongchaogroup/yijielin/star_sj/PN_Mz19/{srr}_SJ.out.tab"
            df = pd.read_csv(srr_file, sep="\t", header=None)

            # Assign column names based on STAR's SJ.out.tab structure
            df.columns = [
                "chromosome", "start", "end", "strand", "intron_motif",
                "annotated", "unique_reads", "multi_reads", "max_overhang"
            ]

            # Filter the DataFrame for the specified junction(s)
            jc_dfb1 = df[(df["start"] == jb1) | (df["end"] == jb2)]
            jc_dfb2 = df[(df["start"] == jb3) | (df["end"] == jb4)]
            jc_dfb3 = df[(df["start"] == jb5) & (df["end"] == jb6)]
            jc_dfb4 = df[(df["start"] == jb7) & (df["end"] == jb8)]
            jc_total_b_df = pd.concat([jc_dfb1,jc_dfb2,jc_dfb3,jc_dfb4], ignore_index=True)

            # For a: Filter the DataFrame for the specified junction(s)
            jc_dfa1 = df[(df["start"] == ja1) | (df["end"] == ja2)]
            jc_dfa2 = df[(df["start"] == ja3) & (df["end"] == ja4)]
            jc_dfa3 = df[(df["start"] == ja5) & (df["end"] == ja6)]
            jc_dfa4 = df[(df["start"] == ja7) & (df["end"] == ja8)]
            jc_dfa5 = df[(df["start"] == ja9) & (df["end"] == ja10)]
            jc_dfa6 = df[(df["start"] == ja11) & (df["end"] == ja12)]
            jc_dfa7 = df[(df["start"] == ja13) & (df["end"] == ja14)]
            jc_dfa8 = df[(df["start"] == ja15) & (df["end"] == ja16)]
            jc_dfa9 = df[(df["start"] == ja17) & (df["end"] == ja18)]
            jc_total_a_df = pd.concat([jc_dfa1,jc_dfa2,jc_dfa3,jc_dfa4,jc_dfa5,jc_dfa6,jc_dfa7,jc_dfa8,jc_dfa9], ignore_index=True)

            # Check if the junction exists and get the unique read count
            if not jc_total_a_df.empty:
                if not jc_total_b_df.empty:
                    unique_reads_b = jc_total_b_df.iloc[0]["unique_reads"]
                    Junction_count_b.append(unique_reads_b)
                    print(f"{srr} found with {unique_reads_b} reads for b.")
                else:
                    Junction_count_b.append(0)
                    if processed == 0:
                        processed = 1
                        print("code is normal, a file is processed with 0 reads mapped to the junction")
            else:
                Junction_count_b.append(10000000000000)
        except Exception as e:
            print(f"Error processing {srr}: {e}")
            Junction_count_b.append(1000000000000)  # Assign 0 for problematic files

# Create a DataFrame from the results and save to CSV
jc_df_b = pd.DataFrame({
    'SRR': SRR_list,
    'Junction_count': Junction_count_b
})

final_df1 = jc_df_b[jc_df_b["Junction_count"] == 0]

final_df1.to_csv(outpath1, index=False)
print(f"Results saved to: {outpath1}")

final_df2 = jc_df_b[(jc_df_b["Junction_count"] >= 10) & (jc_df_b["Junction_count"] <= 10000)]

final_df2.to_csv(outpath2, index=False)
print(f"Results saved to: {outpath1}")







