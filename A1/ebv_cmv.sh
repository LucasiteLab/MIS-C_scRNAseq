#!/bin/bash

#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=96G
#SBATCH --time=3-00:00:00
#SBATCH --job-name ebv_cmv
#SBATCH --output ebv_cmv_%J.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<ngravindra@gmail.com>

# gw92 fpaths
MMT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-mm10-3.0.0/"
MNT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/mm10-3.0.0_GFP_premrna/"
HHT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-GRCh38-3.0.0/"
HNT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/GRCh38-3.0.0.premrna/"
HMT="/ycga-gpfs/apps/bioinfo/genomes/10xgenomics/refdata-cellranger-hg19_and_mm10-1.2.0"
RMT="/gpfs/ycga/apps/bioinfo/genomes/10xgenomics/refdata-cellranger-macaca-mulatta/macaca_mulatta"
GMT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/mm10-3.0.0_GFP/"
CHK="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/Gallus_gallus_NCBI_build3.1/gallus_ncbi_3.1"
TBB="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/TBB927"
TCO="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/TcongolenseIL3000/"
DRT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/Danio_rerio.GRCz11"
CEG="/ycga/sequencers/pacbio/gw92/10x/reference/c_elegans_ws268/"
HRV="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/hg19_HRV/"
HRP="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/GRCh38_RPVIRUS/"
HSV="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/GRCh38_HSV/"
RHT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/Rattus_norvegicus.Rnor_6.0.95/"
VHT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0/"
VMT="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-vdj-GRCm38-alts-ensembl-3.1.0"
ZSG="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/mm10-3.0.0_ZSG/"
H2B="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/mm10-3.0.0_H2B/"
ZSG1="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/mm10-3.0.0_ZSG1/"
HCO="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/hum_scv2_exons/"
HCU="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/hum_scv2_exons_3utr/"
HSU="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/hum_scv2_S-3utr/"
HCR="/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/hum_scv2_orf1-10/"

# ngr specified files
human_fasta=/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa
human_gtf=/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

mouse_fasta=/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-mm10-3.0.0/fasta/genome.fa
mouse_gtf=/gpfs/ycga/sequencers/pacbio/gw92/10x/reference/refdata-cellranger-mm10-3.0.0/genes/genes.gtf

scv2_fasta=/home/ngr4/apps/references/MT020880.1.fasta
scv2_gtf=/home/ngr4/apps/references/scv2_orf1-10.gtf

hACE2_fasta=/home/ngr4/apps/references/hACE2.fa
hACE2_gtf=/home/ngr4/apps/references/hACE2.gtf

mNG_fasta=/home/ngr4/apps/references/mNG.fa
mNG_gtf=/home/ngr4/apps/references/mNG.gtf

ebv_fasta=/home/ngr4/apps/references/EBV.fa
ebv_gtf=/home/ngr4/apps/references/EBV.gtf

cmv_fasta=/home/ngr4/apps/references/CMV.fa
cmv_gtf=/home/ngr4/apps/references/CMV.gtf

output_dir=/home/ngr4/scratch60/

# concatenate fastas and gtfs
cat ${human_fasta} ${ebv_fasta} ${cmv_fasta} > ${output_dir}hs_ebv_cmv.fa
cat ${human_gtf} ${ebv_gtf} ${cmv_gtf}> ${output_dir}hs_ebv_cmv.gtf

# cellranger
module load cellranger/3.1.0

cellranger mkref --genome=hs_ebv_cmv --fasta=${output_dir}hs_ebv_cmv.fa --genes=${output_dir}hs_ebv_cmv.gtf