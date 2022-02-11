SHELL=/bin/bash

DESIRED_HITS=15181508
MAXIMUM_POSSIBLE_HITS=15181508 # -> wc -l ${INTERVALS_BED}
# WINDOW_SIZE=2000
WINDOW_SIZE=12400
BIN_SIZE=200
# WINDOW_SPAN_COUNT=10 # -> WINDOW_SIZE / BIN_SIZE
WINDOW_SPAN_COUNT=62 # 124 # 62 * 2
DATA_DIR=${PWD}/../data
RESULTS_DIR=${PWD}/../results
# INTERVALS_GZ=${DATA_DIR}/intervals.txt.gz
# INTERVALS_BED=${DATA_DIR}/intervals.bed
# INTERVALS_GZ=${DATA_DIR}/pairwiseMetrics_cancer_non-cancer_s1.txt.gz
# INTERVALS_BED=${DATA_DIR}/pairwiseMetrics_cancer_non-cancer_s1.bed
# INTERVALS_GZ=${RESULTS_DIR}/highValueData.txt.gz
# INTERVALS_BED=${RESULTS_DIR}/highValueData.bed
# INTERVALS_GZ=${RESULTS_DIR}/flattenedData.txt.gz
# INTERVALS_BED=${RESULTS_DIR}/flattenedData.bed
INTERVALS_GZ=${DATA_DIR}/All_833_biosamples.S1.scores.txt.gz
INTERVALS_BED=${DATA_DIR}/All_833_biosamples.S1.scores.bed
ASSEMBLY=hg19
BUILD_BED=${DATA_DIR}/${ASSEMBLY}.bed
WINDOWS_BED=${DATA_DIR}/windows.bed
PQ_RESULTS_BED=${RESULTS_DIR}/pq.bed
PQ_MAX_RESULTS_BED=${RESULTS_DIR}/pq.maxhits.bed
PQ_MAX_WITH_KTH_PERC_RESULTS_BED=${RESULTS_DIR}/pq.maxhitsWithKthPercentile.bed
WIS_RESULTS_BED=${RESULTS_DIR}/wis.bed
MAX_MEAN_RESULTS_BED=${RESULTS_DIR}/maxMean.bed
 
KTH_PERC := $(shell cut -f5 ${INTERVALS_BED} | Rscript -e 'quantile(as.numeric(readLines("stdin")), c(0.99))' | tail -1)

all:

#
# for disjoint-windows
#

# INPUT_FN=${DATA_DIR}/windows.fixed.25k.bed
LEFT_WINDOW_WIDTH=12000
RIGHT_WINDOW_WIDTH=12000
EXTEND_WIDTH=12000
WINDOW_WIDTH=25000
INTERVAL_STEP=1000
# WINDOWS_BED=${INPUT_FN}
# WINDOW_SPAN_COUNT=60

extend_windows:
	${PWD}/extend_windows.sh ${LEFT_WINDOW_WIDTH} ${RIGHT_WINDOW_WIDTH} ${WINDOW_WIDTH} ${DATA_DIR}/removeOverlap.bed.starch ${INPUT_FN}

assembly:
	fetchChromSizes ${ASSEMBLY} \
		| grep -v ".*_.*" \
		| grep -v chrM \
		| awk -v FS="\t" -v OFS="\t" '{ print $$1, "0", $$2 }' \
		| sort-bed - \
		> ${BUILD_BED}

#	${PWD}/prep.sh ${WINDOW_SIZE} ${WINDOW_SIZE} ${BIN_SIZE} ${BUILD_BED} ${INTERVALS_GZ} ${INTERVALS_BED} ${WINDOWS_BED}
prep:
	${PWD}/prep.4column.sh ${WINDOW_SIZE} ${WINDOW_SIZE} ${BIN_SIZE} ${BUILD_BED} ${INTERVALS_GZ} ${INTERVALS_BED} ${WINDOWS_BED}

intervals_summary:
	cut -f5 ${INTERVALS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

max_mean_window_sweep:
	mkdir -p ${RESULTS_DIR}
	${PWD}/max_mean_window_sweep.py --input ${WINDOWS_BED} --k ${DESIRED_HITS} --window-span ${WINDOW_SPAN_COUNT} \
		| sort-bed --tmpdir ${PWD} --max-mem 2G - \
		> ${MAX_MEAN_RESULTS_BED}

priority_queue:
	mkdir -p ${RESULTS_DIR}
	${PWD}/max_heap.py --input ${WINDOWS_BED} --k ${DESIRED_HITS} --window-span ${WINDOW_SPAN_COUNT} \
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		> ${PQ_RESULTS_BED}

priority_queue_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${PQ_RESULTS_BED}

priority_queue_summary:
	cut -f5 ${PQ_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

priority_queue_max_hits:
	mkdir -p ${RESULTS_DIR}
	${PWD}/max_heap.py --input ${WINDOWS_BED} --k ${MAXIMUM_POSSIBLE_HITS} --window-span ${WINDOW_SPAN_COUNT} \
		> ${PQ_MAX_RESULTS_BED}
#	${PWD}/max_heap.py --input ${WINDOWS_BED} --k ${MAXIMUM_POSSIBLE_HITS} --window-span ${WINDOW_SPAN_COUNT} \
#		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
#		> ${PQ_MAX_RESULTS_BED}

priority_queue_max_hits_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${PQ_MAX_RESULTS_BED}

priority_queue_max_hits_summary:
	cut -f5 ${PQ_MAX_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

priority_queue_max_hits_with_kth_percentile:
	mkdir -p ${RESULTS_DIR}
	${PWD}/max_heap_with_kth_percentile.py --input ${WINDOWS_BED} --k ${MAXIMUM_POSSIBLE_HITS} --window-span ${WINDOW_SPAN_COUNT} --kth-percentile ${KTH_PERC} > .priority_queue_max_hits_with_kth_percentile.tmp

	awk -v FS="\t" -v OFS="\t" '($$6=="+")' .priority_queue_max_hits_with_kth_percentile.tmp \
		| bedmap --delim '\t' --echo-map-range --echo-map-id --echo-map-score - \
		| sort-bed - \
		| awk -v FS="\t" -v OFS="\t" '{ print $$0, "+" }' \
		| ${PWD}/max_heap_with_kth_percentile_get_max_state_for_max_score.py \
		> .priority_queue_max_hits_with_kth_percentile.tmp.for.map
	
	bedops --merge .priority_queue_max_hits_with_kth_percentile.tmp.for.map \
		| bedmap --delim '\t' --multidelim '|' --echo --echo-map-id --echo-map-score - .priority_queue_max_hits_with_kth_percentile.tmp.for.map \
		| awk -v FS="\t" -v OFS="\t" '{ print $$0, "+" }' \
		| ${PWD}/max_heap_with_kth_percentile_reduce_map.py ${WINDOW_SIZE} ${BIN_SIZE} \
		> .priority_queue_max_hits_with_kth_percentile.tmp.for.map.merge.map

	awk -v FS="\t" -v OFS="\t" '($$6=="-")' .priority_queue_max_hits_with_kth_percentile.tmp \
		| bedmap --delim '\t' --echo-map-range --echo-map-id --echo-map-score - \
		| sort-bed - \
		| awk -v FS="\t" -v OFS="\t" '{ print $$0, "+" }' \
		| ${PWD}/max_heap_with_kth_percentile_get_max_state_for_max_score.py \
		> .priority_queue_max_hits_with_kth_percentile.tmp.rev.map
	
	bedops --merge .priority_queue_max_hits_with_kth_percentile.tmp.rev.map \
		| bedmap --delim '\t' --multidelim '|' --echo --echo-map-id --echo-map-score - .priority_queue_max_hits_with_kth_percentile.tmp.rev.map \
		| awk -v FS="\t" -v OFS="\t" '{ print $$0, "-" }' \
		| ${PWD}/max_heap_with_kth_percentile_reduce_map.py ${WINDOW_SIZE} ${BIN_SIZE} \
		> .priority_queue_max_hits_with_kth_percentile.tmp.rev.map.merge.map
	
	sort-bed .priority_queue_max_hits_with_kth_percentile.tmp.for.map.merge.map .priority_queue_max_hits_with_kth_percentile.tmp.rev.map.merge.map \
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		| gzip -c \
		> .priority_queue_max_hits_with_kth_percentile.tmp.union.gz

	${PWD}/prep_kth_percentile_final_pass.sh ${WINDOW_SIZE} ${WINDOW_SIZE} ${BIN_SIZE} ${BUILD_BED} .priority_queue_max_hits_with_kth_percentile.tmp.union.gz .priority_queue_max_hits_with_kth_percentile.tmp.union.bed .priority_queue_max_hits_with_kth_percentile.tmp.windows

	${PWD}/max_heap_with_kth_percentile.py --is-final-pass --input .priority_queue_max_hits_with_kth_percentile.tmp.windows --k ${MAXIMUM_POSSIBLE_HITS} --window-span ${WINDOW_SPAN_COUNT} \
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		> ${PQ_MAX_WITH_KTH_PERC_RESULTS_BED}

	rm -f .priority_queue_max_hits_with_kth_percentile.tmp*

priority_queue_max_hits_with_kth_percentile_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${PQ_MAX_WITH_KTH_PERC_RESULTS_BED}

priority_queue_max_hits_with_kth_percentile_summary:
	cut -f5 ${PQ_MAX_WITH_KTH_PERC_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

wis:
	mkdir -p ${RESULTS_DIR}
	${PWD}/weighted_interval_scheduling_via_iteration.py --input ${WINDOWS_BED} --k ${DESIRED_HITS} --window-span ${WINDOW_SPAN_COUNT} \
		> ${WIS_RESULTS_BED}
#	${PWD}/weighted_interval_scheduling_via_iteration.py --input ${WINDOWS_BED} --k ${DESIRED_HITS} --window-span ${WINDOW_SPAN_COUNT} \
#		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
#		> ${WIS_RESULTS_BED}

wis_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${WIS_RESULTS_BED}

wis_summary:
#	cut -f5 ${WIS_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'
	cut -f4 ${WIS_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

clean:
	rm -f ${WINDOWS_BED} ${INTERVALS_BED} ${PQ_RESULTS_BED}
