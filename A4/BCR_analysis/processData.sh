# Kenneth B. Hoehn
# 2/1/2021
# Run BCR data through immcantation pipeline
# First, open immcantation Docker image to run data
# docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.0.0 bash
# toggle which sets to run by commenting out dirs lines
 
dirs=("raw/32M-BCR_VHT_cellranger" "raw/34F-BCR_VHT_cellranger" "raw/35F-BCR_VHT_cellranger" "raw/36M-BCR_VHT_cellranger" "raw/111-BCR_VHT_cellranger" "raw/112-BCR_VHT_cellranger" "raw/113-BCR_VHT_cellranger" "raw/CITE1_BCR_VHT_cellranger" "raw/CITE2_BCR_VHT_cellranger" "raw/CITE3_BCR_VHT_cellranger" "raw/CITE4_BCR_VHT_cellranger" "raw/CITE5_BCR_VHT_cellranger" "raw/CITE6_BCR_VHT_cellranger" "raw/HA5876BLD_Bcell_VHT" "raw/HA5877BLD_Bcell_VHT" "raw/HA5894BLD-BCR_VHT" "raw/HA5952BLD-BCR_VHT" "raw/HA5953BLD-BCR_VHT" "raw/HA5957BLD-BCR_VHT" "raw/NS0A_BCR_VHT_cellranger" "raw/NS0B_BCR_VHT_cellranger" "raw/NS1A_BCR_VHT_cellranger" "raw/NS1B_BCR_VHT_cellranger" "raw/TP6A_BCR_VHT_cellranger" "raw/TP6B_BCR_VHT_cellranger" "raw/TP7A_BCR_VHT_cellranger" "raw/TP7B_BCR_VHT_cellranger" "raw/TP8B_BCR_VHT_cellranger" "raw/TP9B_BCR_VHT_cellranger" "raw/TS2A_BCR_VHT_cellranger" "raw/TS2B_BCR_VHT_cellranger" "raw/TS3A_BCR_VHT_cellranger" "raw/TS3B_BCR_VHT_cellranger" "raw/TS4A_BCR_VHT_cellranger" "raw/TS4B_BCR_VHT_cellranger" "raw/TS5A_BCR_VHT_cellranger" "raw/TS5B_BCR_VHT_cellranger")

for wdir in "${dirs[@]}"
do
	echo $wdir
	changeo-10x -s $wdir/filtered_contig.fasta -a $wdir/filtered_contig_annotations.csv -o . \
    	-g human -t ig -p 3 -o $wdir
done

