#!/usr/bin/env python

import sys

FIELD_DELIM = ';'

for line in sys.stdin:
    (chrom, start, stop, ids_str, scores_str, strand) = line.rstrip().split('\t')
    ids = ids_str.split(FIELD_DELIM)
    scores = [float(x) for x in scores_str.split(FIELD_DELIM)]
    max_score = max(scores)
    max_score_index = scores.index(max_score)
    max_id = ids[max_score_index]
    id = max_id
    score = str(max_score)
    sys.stdout.write('{}\n'.format('\t'.join([chrom, start, stop, id, score, strand])))