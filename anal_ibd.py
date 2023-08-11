#!/usr/local/bin/python3
# for now, look at spatial distribution of relatedness for a given year
#  hmm results should already be restricted to good mono samples
#  outputs list of pairs that pass thresholds for relatedness

from ci import ci_calc
from collections import defaultdict, Counter
from itertools import combinations 
import math
import numpy as np
import random
import seaborn as sns
import sys

def main() :
  if len(sys.argv) != 3 : sys.exit('usage: anal_ibd.py <file base> <target year>')
  infile = 'output/' + sys.argv[1] + '_hmm_filtered.txt'
  target_year = int(sys.argv[2])
  outfile1 = 'results/' + sys.argv[1] + '_hmm_space1.txt'
  outfile2 = 'results/' + sys.argv[1] + '_hmm_space2.txt'
  outfile_bcode = 'output/' + sys.argv[1] + '_bcode_added.txt'
  sampfile = 'output/' + sys.argv[1] + '_good_mono_nonclone_samples.txt'
  polyfile = 'output/' + sys.argv[1] + '_all_poly_samples.txt'
  target_sites = ('SLP', 'SES', 'KDG', 'TBA', 'DEG', 'SMS', 'VLG', 'KLD', 'TAM')
  cols = sns.color_palette("colorblind", n_colors=3)
  max_bcode_diffs = 0     # to accept as highly related
  all_sites = set()       # all sites with either barcode or seq data
  all_samps = defaultdict(set)
  all_samps_w_rel = defaultdict(set)    # any kind of relative, within same site
  all_mono_samps = defaultdict(set)     # bcode or seq, anything that can be used to measure relatedness
  
  # Get polygenomic samples for pruning barcodes
  poly_samps = set()  
  with open(polyfile, 'r') as polyf :
    for line in polyf :
      samp = line.rstrip()
      year = int(samp.split('_')[2])
      if year != target_year : continue
      poly_samps.add(samp)

  # get list of all good, mono samples for the target year and count the number per site
  sitecounts = Counter()       # number of samples by site -- from sequencing
  sampsite = {}
  nsamp = 0
  global_seq_samps = set()
  seq_samp_sets = defaultdict(lambda: set())   # keyed by site, set of sequenced samples, for CI calc -- I may store only for selected sites
  with open(sampfile, 'r') as sampf :
    for line in sampf :
      samp = line.rstrip()
      pieces = samp.split(sep='_')
      site = pieces[1]
      year = int(pieces[2])
      if year != target_year : continue
      sampsite[samp] = site
      nsamp += 1
      sitecounts[site] += 1  
      if site in target_sites : seq_samp_sets[site].add(samp)
      all_sites.add(site)
      global_seq_samps.add(samp)
      all_mono_samps[site].add(samp)
      
  print(nsamp, 'samples, {:.0f} pairs in complete data set'.format(nsamp * (nsamp-1) / 2))

  # Read in the barcode data for all target year samples
  barcodes = {}
  site_bcode_samp_sets = defaultdict(lambda : set())
  barfile = 'data/mono_barcodes_filtered_2019.csv'
  barf = open(barfile, 'r')
  head = barf.readline().rstrip().split(sep=',')
  idx = {col: i for i, col in enumerate(head)}
  nskipped = 0
  for line in barf :
    pieces = line.rstrip().split(sep=',')
    samp = pieces[ idx['Sample_Name'] ]
    sub = samp.split(sep='_')
    site = sub[1]
    sampsite[samp] = site
    year = int(sub[2])
    if year != target_year : continue
    if samp in poly_samps :
      #print('skipping barcode because poly by seq', samp)
      nskipped += 1
      continue
    bcode = pieces[ idx['Barcode_String'] ]
    nN = bcode.count('N')
    nX = bcode.count('X')
#    if nN + nX > 1 : continue
    site_bcode_samp_sets[site].add(samp)
    all_sites.add(site)
    barcodes[samp] = bcode
    all_samps[site].add(samp)
    all_mono_samps[site].add(samp)
  barf.close()
  print('skipped', nskipped, 'barcodes because seq says poly')

  # Do the pairwise comparisons for barcodes
  count = 0
  bcode_samps_used = defaultdict(lambda: set())
  bcode_samps = barcodes.keys()
  site_bcode_ident = Counter()    # pairs for given site in target year with identical barcodes
  site_bcode_tried = Counter()  # all barcode pairs tried for given site
  site2_bcode_ident = Counter()   # same for pairs of sites, indexed by tuple of sites
  site2_bcode_tried = Counter()
  site_bcode_ident_sets = defaultdict(lambda : set())    # keyed by site, set of sample pairs that are identical
  site_bcode_tried_sets = defaultdict(lambda : set())    # keyed by site, set of sample pairs that were tried
  with open(outfile_bcode, 'w') as bcoutf :
    print('samp1\tsamp2\tndiff', file=bcoutf)
    for samp1, samp2 in combinations(bcode_samps, 2) :
      count += 1
      bcode1 = barcodes[samp1]
      bcode2 = barcodes[samp2]
      site1 = sampsite[samp1]
      site2 = sampsite[samp2]
      ndiff = 0
      ncomp = 0
      bcode_samps_used[site1].add(samp1)
      bcode_samps_used[site2].add(samp2)
      for ci,cj in zip( bcode1, bcode2 ) :
        if ci == 'N' or ci == 'X' or cj == 'N' or cj == 'X' :
          continue
        ncomp += 1
        if ci != cj : ndiff += 1
      print(samp1, samp2, str(ndiff), sep='\t', file=bcoutf)
      if site1 == site2 :
        site_bcode_tried_sets[site1].add( (samp1, samp2) )
        site_bcode_tried_sets[site1].add( (samp2, samp1) )
        site_bcode_tried[site1] += 1
        if ndiff <= max_bcode_diffs :
          site_bcode_ident[site1] += 1
          site_bcode_ident_sets[site1].add( (samp1, samp2) )
          site_bcode_ident_sets[site1].add( (samp2, samp1) )
          all_samps_w_rel[site1].add(samp1)
          all_samps_w_rel[site1].add(samp2)
      else :
        site2_bcode_tried[(site1, site2)] += 1
        site2_bcode_tried[(site2, site1)] += 1
        if ndiff <= max_bcode_diffs :
          if site1 in ('TAO', 'TAM') and site2 in ('TAO', 'TAM') : print(samp1, samp2)
          if site1 in ('ROB', 'KAF') and site2 in ('ROB', 'KAF') : print(samp1, samp2)
          if site1 in ('ROB', 'TBA') and site2 in ('ROB', 'TBA') : print(samp1, samp2)
          site2_bcode_ident[(site1, site2)] += 1
          site2_bcode_ident[(site2, site1)] += 1

  print('bcode pairs:', count)
  
  outf1 = open(outfile1, 'w')
  outf2 = open(outfile2, 'w')

  nident = 0
  nrelated = 0
  site_n_related = Counter()  # all pairs (all of which should be related) with ibd  <= .98
  site2_n_related = Counter()
  site_partrel_set = defaultdict(lambda : set())     # keyed by site, set of samp pairs that are part related
  
  with open(infile, 'r') as inf :
    idx = {}
    head = inf.readline().rstrip().split()
    for i, col in enumerate(head) :
      idx[col] = i
    for line in inf :
      pieces = line.rstrip().split()
      samp1 = pieces[ idx['sample1'] ]
      samp2 = pieces[ idx['sample2'] ]
      if samp1 not in global_seq_samps : continue    # this will get rid of anything not from target year
      if samp2 not in global_seq_samps : continue
      ngen = float(pieces[ idx['N_generation'] ])
      ntrans = float(pieces[ idx['N_state_transition'] ])
      fIBD = float(pieces[ idx['fract_sites_IBD'] ])

      site1 = sampsite[samp1]
      site2 = sampsite[samp2]
      if site1 == site2 :
        if fIBD > 0.98 :
          print('huh?')
        else :
          site_partrel_set[site1].add( (samp1, samp2) )
          site_partrel_set[site1].add( (samp2, samp1) )
          site_n_related[site1] += 1
      else :
        if fIBD > 0.98 :
          print('skipping clones between', site1, site2)
          continue
        else :
          site2_n_related[(site1, site2)] += 1
          site2_n_related[(site2, site1)] += 1
  
  sites = sorted(list(all_sites))
  print(*sites)
  # Related pairs within each site
  print('site\tN_samples\tN_pairs\tN_pairs_part_rel\tfract_pairs_part_rel\tf_samp_samesite_anyrel', end='', file=outf1)
  print('\tN_bcode_samps\tN_bcode_pairs\tN_bcode_ident\tfract_ident\tci_lo_fpartrel\tci_hi_fpartrel\tci_lo_fident', end='', file=outf1)
  print('\tci_hi_fident\ttotal_mono_samps', file=outf1)
  for site in sites :      
    n = sitecounts[site]
    npair = round(n * (n-1) / 2)
    frelat = 0
    if npair > 0 :
      frelat = site_n_related[site] / npair
    nbc_pairs = site_bcode_tried[site]
    ident_frac_bc = 0
    if nbc_pairs > 0 :
      ident_frac_bc = site_bcode_ident[site] / nbc_pairs
    ci_lo_partrel = ci_hi_partrel = frelat
    ci_lo_ident = ci_hi_ident = ident_frac_bc
    ident_mean = ident_frac_bc
    if site in target_sites :    # calculate CI for these sites
      ci_lo_partrel, ci_hi_partrel = ci_calc(site_partrel_set[site], seq_samp_sets[site], frelat, 0.55)
      ci_lo_ident, ci_hi_ident = ci_calc(site_bcode_ident_sets[site], site_bcode_samp_sets[site], ident_frac_bc, 0.6)
    print(site, n, npair, site_n_related[site], '{:.3e}'.
          format(frelat), '{:.3e}'.format(len(all_samps_w_rel[site]) / len(all_samps[site])),
          len(bcode_samps_used[site]), nbc_pairs, site_bcode_ident[site],
          '{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:d}'.
          format(ident_mean, frelat-ci_lo_partrel, ci_hi_partrel-frelat, ident_mean-ci_lo_ident, \
                 ci_hi_ident-ident_mean, len(all_mono_samps[site])), file=outf1, sep='\t')

  # Related/identical pairs between sites
  print('site1\tsite2\tN1\tN2\tNpair\tfract_pairs_part_related\tN_bcode_pairs\tN_bcode_ident\tfract_pairs_ident', file=outf2)
  for sitei, sitej in combinations(sites, 2) :
    ni = sitecounts[sitei]
    nj = sitecounts[sitej]
    npair = ni * nj
    nbc_pairs = site2_bcode_tried[(sitei,sitej)]
    frelat = fident_bc = 0
    if npair > 0 :
#      fident = site2_n_ident[(sitei, sitej)] / npair
      frelat = site2_n_related[(sitei, sitej)] / npair
    if nbc_pairs > 0 :
      fident_bc = site2_bcode_ident[(sitei,sitej)] / nbc_pairs
      
    print(sitei, sitej, ni, nj, npair, '{:.3e}'.format(frelat), nbc_pairs, site2_bcode_ident[(sitei,sitej)],'{:.3e}'.
          format(fident_bc), file=outf2, sep='\t')

main()
