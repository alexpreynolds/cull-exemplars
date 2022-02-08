#!/usr/bin/env python

import sys

bin_size = int(sys.argv[1])

for line in sys.stdin:
    elements = line.rstrip().split('\t')
    chrom = elements[0]
    start = int(elements[1])
    stop = int(elements[2])
    remainder = elements[3:]
    if stop - start > bin_size and stop % bin_size == 0 and start % bin_size == 0:
        for new_start in range(start, stop, bin_size):
            new_stop = new_start + bin_size
            sys.stdout.write('{}\t{}\t{}\t{}'.format(chrom, str(new_start), str(new_stop), '\t'.join(remainder)))
    elif stop - start == bin_size:
        sys.stdout.write(line)