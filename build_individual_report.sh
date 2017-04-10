#!/bin/bash
#

OUTPUT_DATA_FOLDER=$1
PROJECT_INPUT_DATA_FOLDER=$2
BLAST_RESULTS_FILE=$3
DATA_SET_NAME=$4

# genscan data retrieval was not run on biocluster, but the data and screen shots from ucsc table browser
# are in the folder below
GENSCAN_DATA_FOLDER="$PROJECT_INPUT_DATA_FOLDER"/genscan_genes/2017-01-24-get-stickleback-genscan-data/whole_gene
GENSCAN_BED_FILE="$GENSCAN_DATA_FOLDER"/stickleback_genscan_01_24_2017.bed
# for our report file we need a file with the genscan column names
GENSCAN_BED_HDR="$GENSCAN_DATA_FOLDER"/stickleback_genscan_bed_cols_header.txt
# we also need to convert the chrI, chrII to groupI, groupII, etc chrom names
# and change the chrUn chrom coordinates to scaffolds
# the converted genscan bed file, the one to be merged into the report is defined below
# it is written to back to the GENSCAN_DATA_FOLDER folder.
GENSCAN_BED_FILE_CONVERTED="$GENSCAN_DATA_FOLDER"/stickleback_genscan_01_24_2017_converted.bed
# we also store the conversion script with the genscan data
GENSCAN_CONV_SCRPT="$GENSCAN_DATA_FOLDER"/stickleback_genscan_convert_chr_to_group_and_unchr_to_scaffolds.sh
# now call the script with the input and output 
WRK_DIR="$OUTPUT_DATA_FOLDER"/temp_wrk_dir_genscan
## NO NEED TO RUN AGAIN # sh $GENSCAN_CONV_SCRPT $GENSCAN_BED_FILE $GENSCAN_BED_FILE_CONVERTED $WRK_DIR

# NSCAN data retrieval was not run on biocluster, but the data and screen shots from ucsc table browser
# are in the folder below
NSCAN_DATA_FOLDER="$PROJECT_INPUT_DATA_FOLDER"/N-SCAN/2017-01-24-get-stickle-back-NSCAN-data/whole_gene
NSCAN_BED_FILE="$NSCAN_DATA_FOLDER"/stickleback_NSCAN_01_24_2017.bed
# for our report file we need a file with the NSCAN column names
NSCAN_BED_HDR="$NSCAN_DATA_FOLDER"/stickleback_NSCAN_bed_cols_header.txt
# we also need to convert the chrI, chrII to groupI, groupII, etc chrom names
# and change the chrUn chrom coordinates to scaffolds
# the converted NSCAN bed file, the one to be merged into the report is defined below
# it is written to back to the NSCAN_DATA_FOLDER folder.
NSCAN_BED_FILE_CONVERTED="$NSCAN_DATA_FOLDER"/stickleback_NSCAN_01_24_2017_converted.bed
# we also store the conversion script with the NSCAN data
NSCAN_CONV_SCRPT="$NSCAN_DATA_FOLDER"/stickleback_NSCAN_convert_chr_to_group_and_unchr_to_scaffolds.sh
# now call the script with the input and output 
WRK_DIR="$OUTPUT_DATA_FOLDER"/temp_wrk_dir_NSCAN
## NO NEED TO RUN AGAIN # sh $NSCAN_CONV_SCRPT $NSCAN_BED_FILE $NSCAN_BED_FILE_CONVERTED $WRK_DIR

# We will need a tmp file for some file manipulation below
TMP="$OUTPUT_DATA_FOLDER"/tmp.txt

# now that we have all the data, use bedtools to decorate the blastx results with TSS, genscan and N-SCAN data
BED_FILE="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME".bed
echo "$BED_FILE"
# remove comment lines beginning with '#' with grep
# remove the header line with tail -n +2
# now use awk to put the 12 standard bed file columns if front 
# for now we have have data for the first 6 bed columns chrom, start, stop, name, score(fpkm), strand and the rest are "."
grep '^#' -v $BLAST_RESULTS_FILE | tail -n +2 | awk -F"\t" -v OFS="\t" -v SET="$DATA_SET_NAME" '{print $1, $2, $3, $4,$5,$6,".",".",".",".",".",".", SET, $7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > "$BED_FILE"

# remove comment lines beginning with '#' with grep
# get the header row with head -n 1
# now use awk to put the 12 standard bed file columns if front of original data
# for now we have have data for the first 4 bed columns chrom, start, stop, name, the rest are "."
HDR_FILE="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr.txt
grep '^#' -v $BLAST_RESULTS_FILE | head -n 1 | awk -F"\t" -v OFS="\t" '{print "chrom","chromStart","chromEnd","name","fpkm","strand","thickStart","thickEnd","itemRbg","blockCount","blockSizes","blockStarts","data_source",$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > $HDR_FILE

module load bedtools/2.21.0
# sort the blastx results
SORTED_BED1="$OUTPUT_DATA_FOLDER"/sorted_bed1.bed
SORTED_BED2="$OUTPUT_DATA_FOLDER"/sorted_bed2.bed
bedtools sort -i $BED_FILE > $SORTED_BED1
# sort the genscan data and combine it with the blastx & tss data
BED_w_GENSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_w_genscan.bed
bedtools sort -i $GENSCAN_BED_FILE_CONVERTED > $SORTED_BED2
bedtools intersect -wao -a $SORTED_BED1 -b $SORTED_BED2 -s > $BED_w_GENSCAN 2>bedtools_error.txt
# combine the genscan header with the blastx & tss header
HDR_w_GENSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr_w_genscan.txt
paste $HDR_FILE $GENSCAN_BED_HDR > $HDR_w_GENSCAN
awk -F"\t" -v OFS="\t" '{print $0,"overlap_[genscan]"}' $HDR_w_GENSCAN > $TMP && mv $TMP $HDR_w_GENSCAN

# sort the N-SCAN data and combine it with the blastx,tss & genscan data
BED_w_GENSCAN_NSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_w_genscan_nscan.bed
bedtools sort -i $NSCAN_BED_FILE > $SORTED_BED2
bedtools intersect -wao -a $BED_w_GENSCAN -b $SORTED_BED2 -s > $BED_w_GENSCAN_NSCAN 2>>bedtools_error.txt
# combine the N-SCAN header with the blastx & tss & genscan header
HDR_w_GENSCAN_NSCAN="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_hdr_w_genscan_nscan.txt
paste $HDR_w_GENSCAN $NSCAN_BED_HDR > $HDR_w_GENSCAN_NSCAN
awk -F"\t" -v OFS="\t" '{print $0,"overlap_[N-SCAN]"}' $HDR_w_GENSCAN_NSCAN > $TMP && mv $TMP $HDR_w_GENSCAN_NSCAN



#
# TEST AUGUSTUS EXONS
#
# AUGUSTUS data retrieval was not run on biocluster, but the data and screen shots from ucsc table browser
# are in the folder below with the data file and header file.
BED_FILE_CONVERTED=/home/groups/simons/Joe/sb_blast_project/input_data/AUGUSTUS/2017-02-12-get-stickle-back-augustus-data/exon/stickleback_augustus_02-13-17_exons_converted.bed
# for our report file we need a file with the NSCAN column names
BED_HDR=/home/groups/simons/Joe/sb_blast_project/input_data/AUGUSTUS/2017-02-12-get-stickle-back-augustus-data/exon/stickleback_AUGUSTUS_exons_bed_cols_header.txt
sh add_data_to_results.sh $BED_w_GENSCAN_NSCAN $BED_FILE_CONVERTED $HDR_w_GENSCAN_NSCAN $BED_HDR "Augustus_exons" $OUTPUT_DATA_FOLDER

#
# TEST Human Protein Block Exons
#
# Human Protein data retrieval was not run on biocluster, but the data and screen shots from ucsc table browser
# are in the folder below with the data file and header file.
BED_FILE_CONVERTED=/home/groups/simons/Joe/sb_blast_project/input_data/Human_Blast_Track/2017-02-13-get-stickleback-human-protein-track/stickleback_human_proteins_02-13-17_blocks_converted.bed
# for our report file we need a file with the NSCAN column names
BED_HDR=/home/groups/simons/Joe/sb_blast_project/input_data/Human_Blast_Track/2017-02-13-get-stickleback-human-protein-track/stickleback_HG18_blastp_blocks_bed_cols_header.txt
sh add_data_to_results.sh $BED_w_GENSCAN_NSCAN $BED_FILE_CONVERTED $HDR_w_GENSCAN_NSCAN $BED_HDR "HG18_blastp" $OUTPUT_DATA_FOLDER

#
# TEST genscan Exons
#
# genscan exon data retrieval was not run on biocluster, but the data and screen shots from ucsc table browser
# are in the folder below with the data file and header file.
BED_FILE_CONVERTED=/home/groups/simons/Joe/sb_blast_project/input_data/genscan_genes/2017-01-24-get-stickleback-genscan-data/exon/stickleback_genscan_02-13-17_exons_converted.bed
# for our report file we need a file with the NSCAN column names
BED_HDR=/home/groups/simons/Joe/sb_blast_project/input_data/genscan_genes/2017-01-24-get-stickleback-genscan-data/exon/stickleback_genscan_exons_bed_cols_header.txt
sh add_data_to_results.sh $BED_w_GENSCAN_NSCAN $BED_FILE_CONVERTED $HDR_w_GENSCAN_NSCAN $BED_HDR "genscan_exon" $OUTPUT_DATA_FOLDER

#
# TEST N-scan Exons
#
# N-SCAN exon data retrieval was not run on biocluster, but the data and screen shots from ucsc table browser
# are in the folder below with the data file and header file.
BED_FILE_CONVERTED=/home/groups/simons/Joe/sb_blast_project/input_data/N-SCAN/2017-01-24-get-stickle-back-NSCAN-data/exon/stickleback_NSCAN_02-13-17_exons_converted.bed
# for our report file we need a file with the NSCAN column names
BED_HDR=/home/groups/simons/Joe/sb_blast_project/input_data/N-SCAN/2017-01-24-get-stickle-back-NSCAN-data/exon/stickleback_NSCAN_exons_bed_cols_header.txt
sh add_data_to_results.sh $BED_w_GENSCAN_NSCAN $BED_FILE_CONVERTED $HDR_w_GENSCAN_NSCAN $BED_HDR "NSCAN_exons" $OUTPUT_DATA_FOLDER


# put the final header and final data together and write it out.
FINAL_FILE="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_file_to_review.txt
cat $HDR_w_GENSCAN_NSCAN $BED_w_GENSCAN_NSCAN > $FINAL_FILE

# Finally we will add the data for each set together
ALL_DATA_SET_REPORT="$OUTPUT_DATA_FOLDER"/all_data_set_report.txt
if [ ! -f $ALL_DATA_SET_REPORT ];
then
	# if the all data set report does not exist, first add then header
	cat $HDR_w_GENSCAN_NSCAN > $ALL_DATA_SET_REPORT
fi
# now add the data to the report
cat $BED_w_GENSCAN_NSCAN >> $ALL_DATA_SET_REPORT


echo "end of individual data script"