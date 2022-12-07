#!/bin/bash

. /tools/soft/conda3_v1/etc/profile.d/conda.sh

#Starting from base env
#/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022/SD001/bin3c_clust/fasta/*.fna

DATA_PATH="/store/bioinf/analysis/hicmicrobiome/run_2022_01/SDD_2022"
CL_FILES="bin3c_clust/fasta/*.fna"


for name in SD002 SD003 SD005 SD006 SD008 SD009;
do
mkdir ${name}
echo "${DATA_PATH}/${name}/${CL_FILES}"

cat ${DATA_PATH}/${name}/${CL_FILES} > "${name}/all_CL.fna"

#------------------------Plasflow_0.95t--------------------------------------------------
#conda activate plasflow
#mkdir ${name}/plsflow_out
#PlasFlow/PlasFlow.py --input "${name}/all_CL.fna" --output "${name}/plsflow_out/plasflow_predictions.tsv" --threshold 0.95
#cp ${name}/plsflow_out/plasflow_predictions.tsv ${name}/
#conda deactivate

#------------------------Mob-Suite--------------------------------------------------
conda activate mobsuite
mkdir ${name}/mobr_out
mob_recon --infile "${name}/all_CL.fna" --outdir "${name}/mobr_out" -c -f -n 40
less ${name}/mobr_out/contig_report.txt|cut -f5,21 > ${name}/ALL_MOB_REC_${name}.txt

mkdir ${name}/mobt_out
mob_typer --multi --infile "${name}/all_CL.fna" --out_file "${name}/mobt_out/sample_mobtyper_results.txt" -n 40 
less ${name}/mobt_out/sample_mobtyper_results.txt |grep -v "non-mobilizable"| cut -f14,1,17,15,18|awk '{print $5" "$1" "$2" "$7" "$8}'> ${name}/ALL_MOB_TYP_${name}.txt
conda deactivate 

#------------------------Viralverify_-p--------------------------------------------------
conda activate viralverify
mkdir ${name}/viral_v
viralverify -f "${name}/all_CL.fna" --hmm ../databases/pfam/Pfam-A.hmm -o "${name}/viral_v"
#less ${name}/viral_v/all_CL_result_table.csv |awk -F"," '{print $1 "\t" $2 "\t" $3}'> FINAL_PREDICTION_VIRALV_${name}.csv
less ${name}/viral_v/all_CL_result_table.csv|awk -F"," '{print $1","$2","$3","$4}'>${name}/FINAL_PREDICTION_VIRALV_${name}.csv 

conda deactivate
#------------------------blast_plasmids--------------------------------------------------
conda activate hicmag_py37
mkdir ${name}/PLSDB_res
#
blastn -query "${name}/all_CL.fna" -db /store/bioinf/data/PLSDB/plsdb.fna -out ${name}/PLSDB_res/PLSDB_plasmids.txt -outfmt "6 qseqid  sseqid qlen slen evalue pident qcovs length salltitles" -num_threads 40 -evalue 10E-5
cp ${name}/PLSDB_res/PLSDB_plasmids.txt ${name}/PLSDB_plasmids_${name}.txt 
conda deactivate 
done

