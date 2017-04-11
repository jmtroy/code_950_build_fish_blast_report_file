#!/bin/bash

# decorate "file to update" with closest peak
# assume that file to update is t a .bed file
# create another file with the closest peaks

# This script finds the closest H3K4ME3 peak to each exon.
# The 1st input is just the output folder path
# The 2nd input is the FILE_TO_UPDATE with mouse exons 
# The 3rd input is the HDR_TO_UPDATE with the column headers of FILE_TO_UPDATE
# The 4th input is the PEAK file from homer with h3k4me3 peaks
# The 5th input is the data set name

# The peak file was created by merge_anno_peak_pos_files.sh called earlier in the process
# (merge_anno_peak_pos_files.sh is called in the main script and merges 2 homer files into one and does formatting)

# The format of the H3K4ME3_PEAK_FILE file is
# 1 chrom  (1,2,3,...)
# 2 start
# 3 end
# 4 homer id/name
# 5 empty (".")
# 6 strand / meaningless
# 7 annotation from homer
# 8 findPeaks score from homer
# 9 Fold Change vs Control (from Homer)
# 10 p-value vs Control (from Homer)

OUTPUT_DATA_FOLDER=$1
FILE_TO_UPDATE=$2 
HDR_TO_UPDATE=$3 
H3K4ME3_PEAK_FILE=$4
DATA_SET_NAME=$5


echo "closest_peak.sh OUTPUT_DATA_FOLDER $OUTPUT_DATA_FOLDER"
echo "closest_peak.sh FILE_TO_UPDATE     $FILE_TO_UPDATE"
echo "closest_peak.sh HDR_TO_UPDATE      $HDR_TO_UPDATE"
echo "closest_peak.sh H3K4ME3_PEAK_FILE  $H3K4ME3_PEAK_FILE"
echo "closest_peak.sh DATA_SET_NAME      $DATA_SET_NAME"

module load bedtools/2.21.0

tmp1_bed="$OUTPUT_DATA_FOLDER"/closest_peak_tmp1.bed
tmp1_hdr="$OUTPUT_DATA_FOLDER"/closest_peak_tmp1.hdr

#
# use bedtools closest to find the closest H3K4ME3 peak to each exon following the rules below
# 1) if the exon is on the "-" strand the closest peak is the peak on the 3` end
# 2) if the exon is on the "+" strand the closest peak is the peak on the 5` end
# options used...
# -a is the file with exons
# -b is the file with the H3K4ME3 peaks
# -D a include distance in output, upstream features are reported as negative, 
# and use strand of feature in file a to determine what is upstream and downstream
# -id ignore downstream
bedtools closest -a $FILE_TO_UPDATE -b $H3K4ME3_PEAK_FILE -t first -D a -id > $tmp1_bed
mv $tmp1_bed $FILE_TO_UPDATE

# add columns to the header using awk
# use print 0 for the original header and then add the 10 other column headers
# do not use "bedScore.H3K4ME3_PEAK" as column name, use "." for col 5 (jmt 3/28/17)
awk -F $'\t'  'BEGIN {OFS=FS} {print $0, "chrom.H3K4ME3_PEAK", "chromStart.H3K4ME3_PEAK", "chromEnd.H3K4ME3_PEAK", "name.H3K4ME3_PEAK", ".", "Strand.H3K4ME3_PEAK", "homerAnnotation.H3K4ME3_PEAK","homerFindPeaksScore.H3K4ME3_PEAK","homeFoldChangeVsControl.H3K4ME3_PEAK","homerPvalueVsControl.H3K4ME3_PEAK","distance.H3K4ME3_PEAK"}' $HDR_TO_UPDATE > $tmp1_hdr
mv $tmp1_hdr $HDR_TO_UPDATE

# create file with only the peaks that were closest. which is the last 11 columns in 
# the updated header file  (but omit the very last one which is the distance from the exon to the peak)
# and the last 11 columns in the updated exon file ($FILE_TO_UPDATE) (also omit the very last one which is the distance from the exon to the peak)
# added sort -u to get only unique data.
CLOSEST_PEAKS="$OUTPUT_DATA_FOLDER"/"$DATA_SET_NAME"_closest_peaks.txt
awk -F $'\t' 'BEGIN {OFS=FS} {print $(NF-11),$(NF-10),$(NF-9), $(NF-4), $(NF-7),$(NF-6),$(NF-5),$(NF-8),$(NF-3),$(NF-2),$(NF-1),".","dataset"}' $HDR_TO_UPDATE > $CLOSEST_PEAKS
awk -F $'\t' -v dataset=$DATA_SET_NAME 'BEGIN {OFS=FS} {print $(NF-11),$(NF-10),$(NF-9), $(NF-4), $(NF-7),$(NF-6),$(NF-5),$(NF-8),$(NF-3),$(NF-2),$(NF-1), ".",dataset}' $FILE_TO_UPDATE | sort -u >> $CLOSEST_PEAKS
