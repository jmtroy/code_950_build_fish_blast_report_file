#!/bin/bash
#

# uses module load bedtools/2.21.0 - but assume its already loaded
# we also assume MAIN_RESULTS_BED is already sorted.

MAIN_RESULTS_BED=$1
DATA_TO_ADD_BED=$2
MAIN_RESULTS_HEADER=$3
DATA_TO_ADD_HEADER=$4
DATA_TO_ADD_NAME=$5
OUTPUT_DATA_FOLDER=$6

# sort the DATA_TO_ADD_BED and combine it with the MAIN_RESULTS_BED data
TEMP_SORTED_BED="$OUTPUT_DATA_FOLDER"/temp_sorted.bed
TEMP_OUT="$OUTPUT_DATA_FOLDER"/temp_out.bed
# use bedtools to sort file
bedtools sort -i $DATA_TO_ADD_BED > $TEMP_SORTED_BED

# Now do a bedtools intersect to add data to the main_results_bed if it overlaps
# -wao option:  Write the original A and B entries plus the number of base pairs of 
# overlap between the two features. However, A features w/o overlap are also reported 
# with a NULL B feature and overlap = 0
# -s option: Force “strandedness”. That is, only report hits in B that overlap A on the 
#    same strand. By default, overlaps are reported without respect to strand.
bedtools intersect -wao -a $MAIN_RESULTS_BED -b $TEMP_SORTED_BED -s > $TEMP_OUT 2>>bedtools_error.txt
# remove main_results_bed and move contents of temp_out into main_results_bed
rm $MAIN_RESULTS_BED
mv $TEMP_OUT $MAIN_RESULTS_BED

# combine the N-SCAN header with the blastx & tss & genscan header
TEMP_HDR="$OUTPUT_DATA_FOLDER"/tmp_hdr.txt
OVERLAP_COL_HDR="$DATA_TO_ADD_NAME"_genomic_overlap
paste $MAIN_RESULTS_HEADER $DATA_TO_ADD_HEADER | awk -F"\t" -v OFS="\t" -v X=$OVERLAP_COL_HDR '{print $0,X}' > $TEMP_HDR
# remove main_results_bed and move contents of temp_out into main_results_bed
rm $MAIN_RESULTS_HEADER
mv $TEMP_HDR $MAIN_RESULTS_HEADER

# remember to remove TEMP_SORTED_BED
rm $TEMP_SORTED_BED
