#!/usr/bin/env python

'''

Assume a genome-wide signal vector S, and a window of size k.

1. In a continuous sliding window, calculate max() over S within size k
   windows -- make sure the estimator is in the center, i.e. k/2 to the left
   and k/2 to the right. This will result in a step-like signal, indicating
   regions containing high local maxima, i.e. peaks.

2. In the same way, calculate mean() over S within size k windows. This
   constitutes a smooth estimator for regions with generally higher signals
   and will be used to break ties.

3. Order the genome by 1) (primary) and 2) (secondary) -- this way, 2) breaks
   ties in 1), since there will be many of those. Walk through the ordered
   list in a descending manner, and remove elements/windows that overlap
   higher scoring ones.

'''

import sys
import numpy as np
import pandas as pd
from io import StringIO
import click

@click.command()
@click.option('--input', help='input filename')
@click.option('--k', type=int, help='samples')
@click.option('--window-span', type=int, default=10, help='number of windows used for overlap/rejection testing')
def main(input, k, window_span):
    column_names = ['Chromosome', 'Start', 'End', 'ID', 'Score', 'Strand', 'PVal', 'CorrectedPVal', 'Signif']
    df = pd.read_csv(input, sep='\t', header=None, names=column_names, dtype={"ID": "string", "Signif": "string"})
    # df = df.iloc[0:500]
    n = len(df.index)
    r = np.zeros(n, dtype=np.bool_)
    w = df.loc[:, 'Score'].values
    est_max = np.zeros(n)
    np.put(est_max, np.arange(0, window_span), np.nan)
    np.put(est_max, np.arange(n - window_span, n), np.nan)
    est_mean = np.zeros(n)
    np.put(est_mean, np.arange(0, window_span), np.nan)
    np.put(est_mean, np.arange(n - window_span, n), np.nan)
    for idx in range(window_span, len(w) - window_span + 1):
        window = w[idx - window_span:idx + window_span]
        est_max[idx] = np.amax(window)
        est_mean[idx] = np.mean(window)
    df["EstMax"] = est_max
    df["EstMean"] = est_mean
    df["OriginalIndex"] = range(n)
    df.sort_values(["EstMax", "EstMean"], ascending=[False, False], inplace=True)
    df = df.reset_index(drop=True)

    '''
    We don't need a heap queue as the sort order is already set.
    We just need to walk through the data frame row by row, test that
    the interval is not excluded (and append it, if so, and exclude
    all intervals over the range). We omit intervals with a zero score,
    and we omit those edge intervals that had a NaN max or mean score.
    '''
    q = []
    for i, v in enumerate(df.index):
        if df.iloc[v]["Score"] == 0 or np.isnan(df.iloc[v]["EstMax"]) or np.isnan(df.iloc[v]["EstMean"]):
            continue
        original_i = df.iloc[v]["OriginalIndex"]
        start_i = (original_i - window_span) if (original_i - window_span) > 0 else 0
        stop_i = (original_i + window_span + 1) if (original_i + window_span + 1) <= n else n
        if not np.any(r[start_i : stop_i]):
            for i in range(start_i, stop_i):
                r[i] = True
            q.append(v)
    
    '''
    Print indices in raw order.
    '''
    o = StringIO()
    df = df.drop(columns=["EstMax", "EstMean", "OriginalIndex"])
    df.iloc[q].to_csv(o, sep='\t', index=False, header=False)
    sys.stdout.write('{}'.format(o.getvalue()))
    
if __name__ == '__main__':
    main()