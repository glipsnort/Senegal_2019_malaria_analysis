from itertools import combinations 
import math
import numpy as np
import random

def ci_calc(related_set, samp_set, mean_val, sub_fract) :
  nsub = 1000
  bounds = (.16, .84)
  #bounds = (.05, .95)
  nsamp = len(samp_set)
  subsize = round(sub_fract * nsamp)
  relrate_list = []
  npairs = round(subsize * (subsize-1) / 2)
  for isub in range(nsub) :
    nrel = 0
    subset = random.sample(list(samp_set), subsize)
    for samp1, samp2 in combinations(subset, 2) :
      if (samp1, samp2) in related_set :
        nrel += 1
    relrate_list.append(nrel/npairs)
  return np.quantile(relrate_list, bounds)

def ci_calc_ident(related1, samps1, related2, pairs2, mean_val1, mean_val2) :
  # This part is the same as the previous function: just subsample the samples (by half) repeatedly and build up a distribution of values
  bounds = (.16, .84)
#  bounds = (.05, .95)
  nsub = 400
  nsamp1 = len(samps1)
  sub_fract1 = .6
  if mean_val1 < 0.1 : sub_fract1 = .5
  if mean_val1 < .05 : sub_fract1 = .4
  subsize1 = round(sub_fract1 * nsamp1)
  npairs1 = 0
  relrate_list1 = np.empty(nsub)
  if subsize1 > 2 :
    npairs1 = round(subsize1 * (subsize1-1) / 2)
    for isub in range(nsub) :
      nrel = 0
      subset = random.sample(list(samps1), subsize1)
      for samp1, samp2 in combinations(subset, 2) :
        if (samp1, samp2) in related1 : nrel += 1
        relrate_list1[isub] = nrel / npairs1

  # This part is a little complicated, because not all pairs are being used for the calculation of the value. Extract the  samples from
  #  the pairs, subsample them, and only examine pairs that were in the original set of pairs
  samps2 = set()
  for pair in pairs2 :
    samps2.add(pair[0])
    samps2.add(pair[1])
   sub_fract2 = .36 
  if mean_val2 < .1 : sub_fract2 = .25
  if mean_val2 < .05 : sub_fract2 = .16
  relrate_list2 = np.zeros(nsub)
  target_pairs = round(sub_fract2 * len(pairs2))
  samplist = list(samps2)
  for isub in range(nsub) :
    npairs2 = 0
    nrel = 0
    nsampled = 0
    random.shuffle(samplist)
    for samp1, samp2 in combinations(samplist, 2) :
      if (samp1, samp2) not in pairs2 : continue    # (both sample orders were entered into pairs2)
      npairs2 += 1
      if (samp1, samp2) in related2 : nrel += 1 
      nsampled += 1
      if nsampled >= target_pairs : break
    rr = 0
    if npairs2 > 0 :
      relrate_list2[isub] = nrel / npairs2

  if npairs1 == 0 :
    relrate_list = relrate_list2
  else :
    w1 = math.sqrt(npairs1)
    w2 = math.sqrt(npairs2)
    relrate_list = (relrate_list1 * w1 + relrate_list2 * w2) / (w1 + w2)
  print(npairs1, npairs2)
  return np.quantile(relrate_list, bounds), np.mean(relrate_list)
