#!/bin/bash
. /tools/soft/conda3_v1/etc/profile.d/conda.sh

DATA_PATH="../"
GRAPH_FHALF_PATH="../obolensky"
GRAPH_SHALF_PATH="assembly_graph.fastg"

conda activate scapp

#for name in B-9817;
for name in SD001;
do
#in scapp folder!
mkdir "${name}_new"

#create assemblies of putative plasmids
scapp -g "${GRAPH_FHALF_PATH}/${name}/${GRAPH_SHALF_PATH}" -o "${name}_new/" -r1 "${DATA_PATH}/${name}_R1.fastq.gz" -r2 "${DATA_PATH}/${name}_R2.fastq.gz" -p 40

#conda deactivate
conda activate hicmag_py37


#get which contigs(cl) were mapped on the assembly and select no duplicated ones
cp /home/alexclear/Pl_V_identify/${name}/all_CL.fna ${name}_new/


makeblastdb -in ${name}_new/assembly_graph.confident_cycs.fasta -dbtype nucl

blastn -query ${name}_new/all_CL.fna -db ${name}_new/assembly_graph.confident_cycs.fasta -out ${name}_new/hzfm6cl.txt -outfmt "6 qseqid  sseqid qlen evalue  pident length qcovs salltitles" -num_threads 30
#chmod +x ${name}_new/hzfm6cl.txt

#get cl blasting with assembly with pid=99 and qcov>60 and get first better hit
less ${name}_new/hzfm6cl.txt |awk '{if ($5>90 && $7>20){print $0}}'|awk '{$(NF--)=""; print}'|awk -v OFS="\t" '{$1=$1; print}'|sort -u -k1,1>${name}_new/results_cl_blast_filtered.txt 

#for joining
less ${name}_new/results_cl_blast_filtered.txt |awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'|sort -n >${name}_new/results_cl_blast_filtered_col_swapted.txt 

#what organisms are from the assembly:
blastn -query ${name}_new/assembly_graph.confident_cycs.fasta -db /store/bioinf/data/blast_DBs/nt/nt.fasta -out ${name}_new/hzfm6ass.txt -outfmt "6 qseqid sseqid qlen slen evalue  pident length qcovs salltitles" -num_threads 60 -evalue 10E-5


#get results with the best cov/id and qlen/slen; 5-pid;7-qcov
less ${name}_new/hzfm6ass.txt |awk 'function abs(v) {return v < 0 ? -v : v}{if ($6>90 && $8>60 && abs(($3-$4)/$4)<=0.09){print $0" "abs(($3-$4)/$4)}}'|awk '_[$1]++ < 1'|sort -n >${name}_new/results_ass_blast_filtered.txt 

join ${name}_new/results_cl_blast_filtered_col_swapted.txt ${name}_new/results_ass_blast_filtered.txt >${name}_new/joined_ass_cl.txt 
less ${name}_new/joined_ass_cl.txt |awk '{print $1"\t"$2"\t"$14"\t"$15" "$16" "$17" "$18" "$19" "$20}'>${name}_new/${name}_final_joined_ass_cl.txt 
conda deactivate

#echo "error"
#exp - seqtk
conda activate exp
less ${name}_new/${name}_final_joined_ass_cl.txt|cut -f1|uniq>${name}_new/good_rnode.txt
seqtk subseq ${name}_new/assembly_graph.confident_cycs.fasta ${name}_new/good_rnode.txt>good_${name}_assembly_graph.confident_cycs.fasta 
conda deactivate

done

#conda deactivate

