#!/usr/bin/env python

import sys
from math import nan
from unicodedata import name
from xml.dom import ValidationErr
import pandas as pd
import numpy as np
import scipy.stats as st

group = "WIS"
middle_bin = 62
fn = "../results/wis.bed"

# group = "MM"
# middle_bin = 62
# fn = "../results/maxMean.bed"

df = pd.read_csv(fn, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score', 'IndexMaxBin'])
df = df.sort_values(by=['Score'], ascending=False)

print('{}'.format('\t'.join(["Group", "TopNStart", "TopNEnd", "Mean", "Median", "Mode", "CI_left", "CI_right", "FreqCenter", "FracFreqCenter"])))

start = 0
end = 100000
bin_size = 1000
for subset_size in range(start, end, bin_size):
    subset = df[subset_size:subset_size+bin_size]['IndexMaxBin'].values
    raw_subset = df[subset_size:subset_size+bin_size]
    raw_subset_centered = raw_subset[raw_subset['IndexMaxBin'] == middle_bin]
    raw_subset_centered_values = raw_subset_centered['IndexMaxBin'].values
    # print(len(raw_subset_centered_values), len(raw_subset_centered_values)/subset_size)
    mean = nan
    median = nan
    mode = nan
    ci = [nan, nan]
    freq_center = nan
    frac_freq_center = nan
    try:
        mean = round(np.mean(subset))
        median = int(np.median(subset))
        mode = st.mode(subset)[0][0]
        ci = [round(x) for x in st.t.interval(alpha=0.95, df=len(subset)-1, loc=np.mean(subset), scale=st.sem(subset))]
        freq_center = len(raw_subset_centered_values)
        freq_center_denominator = bin_size
        frac_freq_center = len(raw_subset_centered_values) / freq_center_denominator
    except (ValueError, IndexError):
        pass    
    
    # print(subset_size, mean, median, mode, ci)
    print('{}'.format('\t'.join([group, str(subset_size), str(subset_size + bin_size), str(mean), str(median), str(mode), str(ci[0]), str(ci[1]), str(freq_center), str(frac_freq_center)])))