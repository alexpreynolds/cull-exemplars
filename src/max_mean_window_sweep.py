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
import bottleneck as bn
import numpy as np
import pandas as pd
from io import StringIO
import click

@click.command()
@click.option('--input', help='input filename')
@click.option('--k', type=int, help='samples')
@click.option('--window-span', type=int, default=10, help='number of windows used for overlap/rejection testing')
def main(input, k, window_span):
    
    sys.stderr.write('A\n')
    # column_names = ['Chromosome', 'Start', 'End', 'ID', 'Score', 'Strand', 'PVal', 'CorrectedPVal', 'Signif']
    column_names = ['Chromosome', 'Start', 'End', 'Score']
    
    
    # df = pd.read_csv(input, sep='\t', header=None, names=column_names, dtype={"ID": "string", "Signif": "string"})
    df = pd.read_csv(input, sep='\t', header=None, names=column_names)
    # df = df.iloc[0:500]
    
    
    n = len(df.index)
    r = np.zeros(n, dtype=np.bool_)
    df['IndexOfMaxVal'] = -1
    
    # w = df.loc[:, 'Score'].values
    w = df.loc[:, 'Score'].to_numpy()
        
    sys.stderr.write('B\n')
    fws = 2 * window_span
    est_max_v2 = bn.move_max(w[window_span-1:], window=fws, min_count=window_span+1)
    est_max_v2 = np.append(est_max_v2, np.repeat(np.nan, window_span - 1))
    est_mean_v2 = bn.move_mean(w[window_span-1:], window=fws, min_count=window_span+1)
    est_mean_v2 = np.append(est_mean_v2, np.repeat(np.nan, window_span - 1))

    df["EstMax"] = est_max_v2
    df["EstMean"] = est_mean_v2
    
    df["OriginalIndex"] = range(n)

    df_copy = df.copy()
    
    df.sort_values(["EstMax", "EstMean"], ascending=[False, False], inplace=True)
    
    df = df.reset_index(drop=True)

    '''
    We don't need a heap queue as the sort order is already set.
    We just need to walk through the data frame row by row, test that
    the interval is not excluded (and append it, if so, and exclude
    all intervals over the range). We omit intervals with a zero score,
    and we omit those edge intervals that had a NaN max or mean score.
    '''
    sys.stderr.write('C [{}]\n'.format(len(df.index)))
    q = []
    for i, v in enumerate(df.index):
        if i % 10000 == 0: sys.stderr.write('C -> {}\n'.format(i))
        if df.iloc[v]["Score"] == 0 or np.isnan(df.iloc[v]["EstMax"]) or np.isnan(df.iloc[v]["EstMean"]):
            continue
        original_i = df.iloc[v]["OriginalIndex"]
        start_i = (original_i - window_span) if (original_i - window_span) > 0 else 0
        stop_i = (original_i + window_span + 1) if (original_i + window_span + 1) <= n else n
        if not np.any(r[start_i : stop_i]):
            # for i in range(start_i, stop_i):
            #     r[i] = True
            np.put(r, np.arange(start_i, stop_i), True)
            q.append(v)
            vals_over_span = df_copy[start_i : stop_i]['Score'].values
            df.at[i, 'IndexOfMaxVal'] = np.argmax(vals_over_span)
            # print(start_i, stop_i)
            # print(vals_over_span)
            # print(np.argmax(vals_over_span))
            # sys.exit(0)
    
    '''
    Print indices in raw order.
    '''
    sys.stderr.write('D\n')
    o = StringIO()
    df["Score"] = df["EstMax"]
    df = df.drop(columns=["EstMax", "EstMean", "OriginalIndex"])
    df.iloc[q].to_csv(o, sep='\t', index=False, header=False)
    sys.stdout.write('{}'.format(o.getvalue()))
    
if __name__ == '__main__':
    main()