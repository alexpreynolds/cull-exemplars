#!/bin/bash

set -euxo pipefail

left_window_size=${1}
right_window_size=${2}
bin_size=${3}
build_bed=${4}
input_gz=${5}
input_bed=${6}
output_bed=${7}

#
# the total width of a padded bin is the bin interval itself, plus
# the width of windows upstream (left) and downstream (right) of the bin 
#
total_width=$((${left_window_size} + ${right_window_size} + ${bin_size}))

#
# extract intervals to sorted bed file
#
gunzip -c ${input_gz} \
    | sort-bed - \
    | ${PWD}/fix_bins.py ${bin_size} \
    > ${input_bed}

#
# print nine-column data with insignificant priority/p-value scores
#
bedops --chop ${bin_size} ${build_bed} \
    | awk -v FS="\t" -v OFS="\t" '{ print $1, $2, $3, ".", "0", ".", "1", "1", "." }' \
    > .placeholder_intervals

#
# write a bed file containing all bins over genome
# keeping only those which are of total desired width
#
bedmap --count --echo --echo-map .placeholder_intervals ${input_bed} \
    | awk -v FS="|" -v OFS="\t" '($1 == 0){ print $2 }($1 == 1){ print $3 }' \
    | bedops --range -${left_window_size}:${right_window_size} --everything - \
    | awk -v FS="\t" -v OFS="\t" -v t=${total_width} '(($3 - $2) == t)' \
    > ${output_bed}

#
# cleanup
#
rm -f .placeholder_intervals