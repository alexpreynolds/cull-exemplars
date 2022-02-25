#!/usr/bin/env python

import sys
from math import nan
from unicodedata import name
from xml.dom import ValidationErr
import pandas as pd
import numpy as np
import scipy.stats as st

# group = "WIS"
# middle_bin = 62
# fn = "../results/wis.bed"

group = "MM"
middle_bin = 62
fn = "../results/maxMean.bed"

df = pd.read_csv(fn, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score', 'IndexMaxBin'])
df = df.sort_values(by=['Score'], ascending=False)

print('{}'.format('\t'.join(["Group", "TopNStart", "TopNEnd", "Mean", "Median", "Mode", "CI_left", "CI_right"])))

start = 0
end = 100000
bin_size = 1000
for subset_size in range(start, end, bin_size):
    subset = df[subset_size:subset_size+bin_size]['Score'].values
    # print(len(raw_subset_centered_values), len(raw_subset_centered_values)/subset_size)
    mean = nan
    median = nan
    mode = nan
    ci = [nan, nan]
    try:
        mean = np.mean(subset)
        median = np.median(subset)
        mode = st.mode(subset)[0][0]
        ci = [x for x in st.t.interval(alpha=0.95, df=len(subset)-1, loc=np.mean(subset), scale=st.sem(subset))]
    except (ValueError, IndexError):
        pass    
    
    # print(subset_size, mean, median, mode, ci)
    print('{}'.format('\t'.join([group, str(subset_size), str(subset_size + bin_size), str(mean), str(median), str(mode), str(ci[0]), str(ci[1])])))