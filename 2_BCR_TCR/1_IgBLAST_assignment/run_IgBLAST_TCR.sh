#!/bin/bash 

# LSF options 

#BSUB -o igblast_tcr-%J_%I-output.log
#BSUB -e igblast_tcr-%J_%I-error.log
#BSUB -q "normal"
#BSUB -G open_target_OTAR2064-grp

## Set PE 

#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M8000
#BSUB -n 1 

## Task range. 
#BSUB -J "igblast_tcr[1-69]" 

# Get list of target samples 
samples=$(awk '{print $1}' /nfs/users/nfs_c/cs54/OTAR2064_cs54/2_IgBLAST_assignment/final_data/slemap_cellranger_to_pool.txt)
pools=$(awk '{print $2}' /nfs/users/nfs_c/cs54/OTAR2064_cs54/2_IgBLAST_assignment/final_data/slemap_cellranger_to_pool.txt)

# Get sample to be processed by *this* task 
# extract the Nth sample in the list of samples, $samples, where N == $LSB_JOBINDEX 
this_sample=$(echo "${samples}" | sed -n ${LSB_JOBINDEX}p)
this_pool=$(echo "${pools}" | sed -n ${LSB_JOBINDEX}p)

input_fasta=/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/singlecell/cellranger/${this_sample}/per_sample_outs/${this_sample}/vdj_t/filtered_contig.fasta
input_annotations=/lustre/scratch126/opentargets/opentargets/OTAR2064/working/slemap/singlecell/cellranger/${this_sample}/per_sample_outs/${this_sample}/vdj_t/filtered_contig_annotations.csv

outdir=/nfs/users/nfs_c/cs54/OTAR2064_cs54/2_IgBLAST_assignment/final_data/IgBLAST_output_TCR

# . ~/my_software/immcantation-venv/bin/activate

module load HGI/softpack/users/cs54/immcantation_igblast/3

AssignGenes.py igblast -s ${input_fasta} -b /nfs/users/nfs_c/cs54/igblast_databases_nov24/igblast \
   --organism human --loci tr --format blast --outdir ${outdir} --outname ${this_pool} --exec ~/my_software/ncbi-igblast-1.22.0/bin/igblastn


MakeDb.py igblast -i ${outdir}/${this_pool}_igblast.fmt7 -s ${input_fasta} \
   -r /nfs/users/nfs_c/cs54/igblast_databases_nov24/germlines/imgt/human/vdj --10x ${input_annotations} --extended --outdir ${outdir} --outname ${this_pool}

module unload HGI/softpack/users/cs54/immcantation_igblast/3
