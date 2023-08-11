#!/usr/local/bin/python3

from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import re
import seaborn as sns
import sys

def main() :

  if len(sys.argv) != 3 :
    sys.exit('usage: hetrate.py <file name base> <"yes"=compare to barcode>')
  het_thresh = .0024   # current definition of polygenomic
  cover_thresh = 0.25   # current def of good sample -- threshold on fraction of genome >= 5x cover
  cols = sns.color_palette("colorblind", n_colors=3)

  base = sys.argv[1]
  do_barcode = False
  if (sys.argv[2] == 'yes') : 
    do_barcode = True
  toutfile = 'output/' + base + '_every_sample.txt'
  toutf = open(toutfile, 'w')
  
  goutfile = 'output/' + base + '_good_mono_samples.txt'
  goutf = open(goutfile, 'w')

  aoutfile = 'output/' + base + '_all_mono_samples.txt'
  aoutf = open(aoutfile, 'w')

  bmoutfile = 'output/' + base + '_bad_mono_samples.txt'
  bmoutf = open(bmoutfile, 'w')

  apoutfile = 'output/' + base + '_all_poly_samples.txt'
  apoutf = open(apoutfile, 'w')
  
  poutfile = 'output/' + base + '_good_poly_samples.txt'
  poutf = open(poutfile, 'w')

  pall_outfile = 'output/' + base + '_good_samples.txt'
  pall_outf = open(pall_outfile, 'w')
  
  hetfile = 'output/' + base + '_het.txt'
  hetf = open(hetfile, 'r')
  head = hetf.readline().rstrip().split()
  idx = {}
  for i in range(len(head)) :
    idx[head[i]] = i

  nsamp = 0
  nsamp_mono = 0
  nsamp_good = 0
  nsamp_mono_good = 0
  nsamp_mono_bad = 0
  hetrate = {}
  hetrate_good = {}
  for line in hetf :
    pieces = line.rstrip().split()
    samp = pieces[0]
    year = samp.split('_')[2]
#    if year != '2019' : continue
    nhet = int(pieces[ idx['N_hetcall'] ])
    ntot = int(pieces[ idx['N_call'] ])
    f_cover = float(pieces[ idx['depth_gt5x'] ])
    het_fract = nhet / ntot
    samp = pieces[ idx['samp'] ]
    if re.search(r'Control', samp) or re.search(r'Ctrl', samp) :
      print('skipping', samp)
      continue
    hetrate[samp] = het_fract
    if samp == 'SEN_SLP_2014_001' : print(het_fract)
    nsamp += 1
    print(samp, file=toutf)
    if het_fract <= het_thresh : 
      print(samp, file=aoutf)
      nsamp_mono += 1
    if f_cover >= cover_thresh :
      hetrate_good[samp] = het_fract
      nsamp_good += 1
      # coverage > threshold and either poly or mono
      print(samp, file=pall_outf)
      
      if het_fract <= het_thresh : 
        # monogenomic and coverage > threshold
        print(samp, file=goutf)
        nsamp_mono_good += 1
    else :    # bad samples
      if het_fract <= het_thresh : 
        # monogenomic and coverage < threshold
        print(samp, file=bmoutf)
        nsamp_mono_bad += 1
    if het_fract > het_thresh : 
      print(samp, file=apoutf)
      if f_cover >= cover_thresh : 
        print(samp, file=poutf)
    
  print('N samples:', nsamp)
  print('N mono samples:', nsamp_mono)
  print('N good samples:', nsamp_good)
  print('N good mono samples:', nsamp_mono_good)
  print('N bad mono samples:', nsamp_mono_bad)

  pdf_file = 'results/' + base + '_hetrate.pdf'
  pp = PdfPages(pdf_file)
  fig, ax = plt.subplots()
  hets = [float(x) for x in hetrate.values()]
  print(len(hets))
  ax.hist(x=hets, bins=np.arange(0,0.03,.0006), rwidth=0.85, color=cols[0])
  ax.set_xlabel('Fraction of heterozygous calls')
  ax.set_ylabel('Number of samples')
  ax.set_title('Distribution of het call rate per sample (all samples)')
#  ax.set_ylim(top=40)
  fig.savefig(pp, format='pdf', bbox_inches='tight')

  fig1, ax1 = plt.subplots()
  hets_good = [float(x) for x in hetrate_good.values()]
  print(len(hets_good))
  ax1.hist(x=hets_good, bins=np.arange(0,0.03,.0005), rwidth=0.85, color=cols[0])
  ax1.set_xlabel('Het call rate')
  ax1.set_ylabel('Number of samples')
  ax1.set_title('Distribution of het call rate per sample (>50% high coverage)')
#  ax1.set_ylim(top=40)
  fig1.savefig(pp, format='pdf', bbox_inches='tight')

  # compare to barcode data
  if do_barcode :
    barfile = 'data/Senegal_AllBarcodes_20200701.txt'
    barf = open(barfile, 'r')
    head = barf.readline().split()
    idx = {}
    for i, col in enumerate(head) :
      idx[col] = i
    
    hetlist = []
    misslist = []
    nlist = []
    nsimple = Counter()
    ncomplex = Counter()
    for line in barf :
      pieces = line.rstrip().split(sep='\t')
      samp = pieces[idx['Sample_Name']]
      bcode = pieces[idx['Barcode_String']]
      nN = bcode.count('N')
      nX = bcode.count('X')
      if samp in hetrate_good :
        hetlist.append(hetrate_good[samp])
        misslist.append(nX)
        nlist.append(nN)
        if hetrate[samp] > het_thresh :
          ncomplex[nN] += 1
        else :
          nsimple[nN] += 1
      else :
        pass
    
    print('N het barcode calls\tN COI=1\tN COI>1\tfrac correct')
    for iN in range(3) :
      print(iN, nsimple[iN], ncomplex[iN], '{:.3f}'.format(nsimple[iN] / (nsimple[iN] + ncomplex[iN])), sep='\t')

#      print('did not find', samp)

    fig2, ax2 = plt.subplots()
    ax2.scatter(hetlist, rand_jitter(nlist), s=5, alpha=.8, color=cols[0])
    ax2.set_title('Het rates, barcode vs seq')
    ax2.set_xlabel('Het call rate in sequencing')
    ax2.set_ylabel('Number of het assays in barcode')
    ax2.set_xscale('log')
    ax2.set_xlim(left=.0001)
    fig2.savefig(pp, format='pdf', bbox_inches='tight')

    fig2, ax2 = plt.subplots()
    ax2.scatter(hetlist, rand_jitter(misslist), s=5, alpha=.8, color=cols[0])
    ax2.set_title('Barcode miss rate vs seq het rate')
    ax2.set_xlabel('Het call rate in sequencing')
    ax2.set_ylabel('Number of failed assays in barcode')
    fig2.savefig(pp, format='pdf', bbox_inches='tight')

  pp.close()

def rand_jitter(l):
  stdev = .01*(max(l)-min(l))
#  stdev = .15
  return l + np.random.randn(len(l)) * stdev



main()
