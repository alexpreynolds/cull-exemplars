# cull-exemplars
Filter lower-scoring, neighboring genomic intervals within a specified distance

## Overview

We apply two methods to filter intervals and meet the following criteria:

1. Find a good, potentially optimal distribution of scored intervals over the input set, which should be some specified distance away from each other (i.e., "mutually compatible").

2. We would like both methods to prioritize high-scoring intervals.

3. We would like 250k such elements out of the full superset, which meet the first two criteria.

## Data

### Input

We start with `data/intervals.txt.gz`, which contains ~4.5M intervals within `hg19` space.

Each interval is a 200nt "bin" with a priority score in the fifth column (other columns are allowed and are ignored during processing).

Not all of `hg19` is covered. All starting bins have a positive, non-zero score. We take advantage of this property to fill in gaps in `hg19` space with "placeholder" or "sentinel" bins, each with a zero score. This property will be exploited to know when the priority queue method is "exhausted", i.e., when there are no more valid bins to filter.

### Output

The output from each of the two methods (discussed below) contains a filtered subset of the input bins, which meet the desired criteria.

### Targets

Running the `assembly` and `prep` targets sets up files required for running the two filtering methods (`priority_queue` and `wis`).

## Methods

### Priority queue

In the `priority_queue` and `priority_queue_max_hits` targets, we put padded bins into a [priority queue](https://en.wikipedia.org/wiki/Priority_queue), ordered by score. 

We pop the highest-scoring element off the queue, and keep it if it does not overlap any other popped elements. If it does, we still keep it, but all lower-scored overlaps within 2kb are marked as rejected.

We pop the next highest-scoring element and ask if it has already been rejected. If it hasn't, we keep it and we again mark any lower-scored overlaps as rejected, if there are any. If it was previously rejected, we skip it and keep popping and testing, until there are no elements left in the queue to pop-and-test, or until we hit some preselected number of intervals (250k).

Intervals we keep are written to standard output.

In the `priority_queue_max_hits` target, we specify a larger `k` value of 4453117 (~4.5M). This is the number of intervals we started with and therefore the maximum number we can get back from any method. In the optimal case, all intervals would be 2kb or more away from each other, but in reality there are overlaps and so some fraction of these will get filtered. In this case, we will therefore eventually pop an element off the queue that has a score of zero. This is an artificial "placeholder" or sentinel bin, so we cease testing at this point and write whatever intervals we have found to standard output.

### Weighted-interval scheduling

In the `wis` target, we use a [dynamic programming](https://en.wikipedia.org/wiki/Interval_scheduling#Weighted) approach to trace a path of high-scoring elements which do not overlap one another within 2kb.

This locally optimal path contains elements which do not overlap within 2kb. If we have more than 250k such elements in this path, we use a priority queue to pull out the top 250k such elements.

Intervals retrieved via this path are written to standard output.

### Validation

The `priority_queue_validate_no_overlaps` target tests that the outputted bins from the `priority_queue` method are not within 2kb of each other. Likewise, the `wis_validate_no_overlaps` target tests that the bins from the `wis` method are not within 2kb of each other.

## Results

Here is a summary of the baseline score distribution of the input elements. This is obtained via the `intervals_summary` target:

```
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.3819   0.9664   2.2576   4.5391   5.3447 102.5106 
```

The filtered subsets should meet or beat these values. 

Here is the score profile of `priority_queue` output, via the `priority_queue_summary` target:

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.962   3.566   6.455  10.689  13.111 102.511 
```

And likewise for the output from the `wis` method, via the `wis_summary` target:

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.3943  1.3949  3.9383  8.4296 10.4911 99.7358 
```

In addition to giving a better score result, the `priority_queue` method ran much faster on this data than the `wis` method.

### Priority queue (max-hits)

In the case of the `priority_queue_max_hits` target output, we get 395649 elements before we exhaust the queue. 

These elements have the following score profile:

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.382   1.206   3.236   7.119   8.480 102.511 
```

Exhausting the queue performs worse than the `wis` method, score-wise, but this method returns more elements (~400k) as opposed to `wis`, which can return at most only ~251k elements when it is run to completion.
