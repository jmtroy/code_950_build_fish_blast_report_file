#!/bin/bash
#

#
# Take a stickleback bed file with groupI, groupII ... scaffold_xxxx, ....
# and convert to bed file with chrI, chrII, ...
#


INPUT_BED=$1
OUTPUT_BED=$2
WRK_DIR=$3
CHRUN_TO_SCAFFOLDS_TXT=/home/groups/simons/Joe/sb_blast_project/input_data/scaffolds_to_unchr_frm_abbas/ChrUn_to_scaffolds.txt

# test if the WRK_DIR variable is unused or empty and if so create a value with the date and time
if [ -z $WRK_DIR ]
then 
	dt=`date +%Y%m%d`
	tm=`date +%H%M%S`
	RUN_ID=
	WRK_DIR=tmp_wrk_dir_"$dt"_"$tm"
fi

# create the WRK_DIR if it does not exist
if ! [ -d "$WRK_DIR" ]
then
  mkdir -p "$WRK_DIR"
fi


# First step, convert the CHRUN_TO_SCAFFOLDS_TXT to a bed file, CHRUN_TO_SCAFFOLDS_BED
SCAFFOLDS_TO_CHRUN_BED="$WRK_DIR"/scaffolds_to_ChrUn.bed

# convert mapping file to a bed file that will be used to aid in the mapping of the chrUN repeats to scaffold locations
# 1) use tail -n +2 "$CHRUN_TO_SCAFFOLDS_TXT"  to remove the header from the original mapping file
# 2) then pipe the results to the awk command to reformat into a bed file.
# 2a) the    -F $'\t'    tells the awk command the file has a tab delimiter
# 2b) the    '{print $3"\t"$4-1"\t"$5"\t"$1"\t"$2}'   tells the awk command how to order the orginal file columns to create the bed file.
# 2c) note the  "\t"  between each column is to tell awk the resulting bed file is tab separated.
# 2d) also note    $4-1   to convert the start location from "1-based" to "0-based" as required for a bed file.  for more info on 1-based and 0-based see: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
tail -n +2 "$CHRUN_TO_SCAFFOLDS_TXT" | awk -F $'\t' '{print $1"\t0\t"$2"\t"$3"\t"$4}'  > "$SCAFFOLDS_TO_CHRUN_BED"

# Step 2 sort the input file to be converted
SORTED_INPUT_BED="$WRK_DIR"/sorted_input.bed
bedtools sort -i $INPUT_BED > $SORTED_BED2



























