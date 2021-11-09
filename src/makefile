SHELL=/bin/bash

DESIRED_HITS=250000
WINDOW_SIZE=2000
BIN_SIZE=200
WINDOW_SPAN_COUNT=10 # -> WINDOW_SIZE / BIN_SIZE
DATA_DIR=${PWD}/../data
INTERVALS_GZ=${DATA_DIR}/intervals.txt.gz
INTERVALS_BED=${DATA_DIR}/intervals.bed
ASSEMBLY=hg19
BUILD_BED=${DATA_DIR}/${ASSEMBLY}.bed
WINDOWS_BED=${DATA_DIR}/windows.bed
RESULTS_DIR=${PWD}/../results
PQ_RESULTS_BED=${RESULTS_DIR}/pq.bed
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
	${PWD}/max_heap.py --input ${WINDOWS_BED} --k ${DESIRED_HITS} --window-span ${WINDOW_SPAN_COUNT} \
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		> ${PQ_RESULTS_BED}

priority_queue_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${PQ_RESULTS_BED}

priority_queue_summary:
	cut -f5 ${PQ_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

wis:
	${PWD}/weighted_interval_scheduling_via_iteration.py --input ${WINDOWS_BED} --k ${DESIRED_HITS} \
		| bedops --range ${WINDOW_SIZE}:-${WINDOW_SIZE} --everything - \
		> ${WIS_RESULTS_BED}

wis_validate_no_overlaps:
	${PWD}/validate_no_overlaps.py ${WINDOW_SIZE} < ${WIS_RESULTS_BED}

wis_summary:
	cut -f5 ${WIS_RESULTS_BED} | Rscript -e 'summary (as.numeric (readLines ("stdin")))'

clean:
	rm -f ${WINDOWS_BED} ${INTERVALS_BED} ${PQ_RESULTS_BED}