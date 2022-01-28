SHELL=/bin/bash

DESIRED_HITS=15181508
MAXIMUM_POSSIBLE_HITS=15181508 # -> wc -l ${INTERVALS_BED}
WINDOW_SIZE=2000
BIN_SIZE=200
WINDOW_SPAN_COUNT=10 # -> WINDOW_SIZE / BIN_SIZE
DATA_DIR=${PWD}/../data
RESULTS_DIR=${PWD}/../results
# INTERVALS_GZ=${DATA_DIR}/intervals.txt.gz
# INTERVALS_BED=${DATA_DIR}/intervals.bed
# INTERVALS_GZ=${RESULTS_DIR}/highValueData.txt.gz
#INTERVALS_BED=${RESULTS_DIR}/highValueData.bed
INTERVALS_GZ=${RESULTS_DIR}/flattenedData.txt.gz
INTERVALS_BED=${RESULTS_DIR}/flattenedData.bed
ASSEMBLY=hg19
BUILD_BED=${DATA_DIR}/${ASSEMBLY}.bed
WINDOWS_BED=${DATA_DIR}/windows.bed
PQ_RESULTS_BED=${RESULTS_DIR}/pq.bed
PQ_MAX_RESULTS_BED=${RESULTS_DIR}/pq.maxhits.bed
WIS_RESULTS_BED=${RESULTS_DIR}/wis.bed

all:

assembly:
	fetchChromSizes ${ASSEMBLY} \
		| grep -v ".*_.*" \
		| grep -v chrM \
		| awk -v FS="\t" -v OFS="\t" '{ print $$1, "0", $$2 }' \
		| sort-bed - \
		> ${BUILD_BED}

prep:
	${PWD}/prep.sh ${WINDOW_SIZE} ${WINDOW_SIZE} ${BIN_SIZE} ${BUILD_BED} ${INTERVALS_GZ} ${INTERVALS_BED} ${WINDOWS_BED}

intervals_summary:
	cut -f5 ${INTERVALS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

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
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		> ${PQ_MAX_RESULTS_BED}

priority_queue_max_hits_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${PQ_MAX_RESULTS_BED}

priority_queue_max_hits_summary:
	cut -f5 ${PQ_MAX_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

wis:
	mkdir -p ${RESULTS_DIR}
	${PWD}/weighted_interval_scheduling_via_iteration.py --input ${WINDOWS_BED} --k 1007638 \
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		> ${WIS_RESULTS_BED}

wis_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${WIS_RESULTS_BED}

wis_summary:
	cut -f5 ${WIS_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

clean:
	rm -f ${WINDOWS_BED} ${INTERVALS_BED} ${PQ_RESULTS_BED}
