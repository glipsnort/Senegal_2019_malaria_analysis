#!/usr/local/bin/python3
#  Simply applies my current definition of 'related' and produces a filtered IBD file

from collections import defaultdict, Counter
from itertools import combinations 
import numpy as np
import sys

def main() :
  if len(sys.argv) != 2 : sys.exit('usage: filter_ibd.py <file base>')
  infile = 'output/' + sys.argv[1] + '.hmm_fract.txt'
  skimfile = 'output/' + sys.argv[1] + '_hmm_filtered.txt'
  sampfile = 'output/' + sys.argv[1] + '_good_mono_samples.txt'
  hetfile = 'output/' + sys.argv[1] + '_het.txt'

  n_rel = n_high = 0
  depth = {}
  with open(hetfile, 'r') as hetf :
    head = hetf.readline().rstrip().split()
    idx = {col: i for i, col in enumerate(head)}
    for line in hetf :
      pieces = line.rstrip().split()
      samp = pieces[idx['samp']]
      depth[samp] = float(pieces[idx['depth_gt5x']])

  samp_set = set()
  with open(infile, 'r') as inf, open(skimfile, 'w') as skimf :
    idx = {}
    head = inf.readline()
    print(head, file=skimf, end='')
    head = head.rstrip().split()
    for i, col in enumerate(head) :
      idx[col] = i
    for line in inf :
      pieces = line.rstrip().split()
      samp1 = pieces[ idx['sample1'] ]
      samp2 = pieces[ idx['sample2'] ]
      year1 = int(samp1.split('_')[2])
      year2 = int(samp2.split('_')[2])
      samp_set.add(samp1)
      samp_set.add(samp2)
      #      ngen = float(pieces[ idx['N_generation'] ])
#      ntrans = float(pieces[ idx['N_state_transition'] ])
      fIBD = float(pieces[ idx['fract_sites_IBD'] ])
      thresh = 0.05
      if depth[samp1] >= 0.5 and depth[samp2] >= 0.5 : 
        thresh = 0.04
      if fIBD > thresh : 
        # My current definition of related pairs
        print(line, file=skimf, end='')
        if year1 == 2019 and year2 == 2019 : 
          n_rel += 1
          if fIBD > 0.5 : n_high += 1

  print(n_high, 'pairs with IBD > 0.5, out of', n_rel, 'related 2019 pairs')
  nsamp = len(samp_set)
  print(nsamp, 'samples, {:.0f} pairs in complete data set'.format(nsamp * (nsamp-1) / 2))

main()
