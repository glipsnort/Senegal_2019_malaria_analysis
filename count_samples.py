#!/usr/local/bin/python3

from collections import Counter

def main() :
  clonefile = 'output/sequenced_clones.txt'
  infiles = ['output/all_good_mono_samples.txt', 'output/all_good_poly_samples.txt']
  outfile = 'output/site_counts.txt'
  prunedfile = 'output/all_good_mono_nonclone_samples.txt'

  include_clones = False
  if include_clones : print('including clones')
  clones = set()
  with open(clonefile, 'r') as clonef :
    for line in clonef :
      clones.add(line.rstrip())
  bysite = {}
  all_seq_samps = []
  all_seq_samps.append( set() )
  all_seq_samps.append( set() )
  for itype, infile in enumerate(infiles) :
    bysite[infile] = Counter()
    with open(infile, 'r') as inf :
      for line in inf :
        samp = line.rstrip()
        pieces = samp.split('_')
        site = pieces[1]
        year = pieces[2]
        if year != '2019' : continue
        if not include_clones and samp in clones : continue
        bysite[infile][site] += 1
        all_seq_samps[itype].add(samp)
  print('n seq samps, mono/poly', len(all_seq_samps[0]), len(all_seq_samps[1]))

  mono_bar_bysite = Counter()
  poly_bar_bysite = Counter()
  mono_bar_samps = set()
  poly_bar_samps = set()
  
  barfile = 'data/mono_barcodes_filtered_2019.csv'
  barf = open(barfile, 'r')
  head = barf.readline().rstrip().split(sep=',')
  idx = {col: i for i, col in enumerate(head)}
  for line in barf :
    pieces = line.rstrip().split(sep=',')
    samp = pieces[ idx['Sample_Name'] ]
    sub = samp.split(sep='_')
    year = int(sub[2])
    site = sub[1]
    if year != 2019 : continue
    bcode = pieces[ idx['Barcode_String'] ]
    nN = bcode.count('N')
    nX = bcode.count('X')
    if nX > 1 : continue
    if samp in all_seq_samps[1] :
      poly_bar_samps.add(samp)
      poly_bar_bysite[site] += 1
    else :
      mono_bar_samps.add(samp)
      mono_bar_bysite[site] += 1
  barfile = 'data/poly_barcodes_filtered_2019.csv'
  barf = open(barfile, 'r')
  head = barf.readline().rstrip().split(sep=',')
  idx = {col: i for i, col in enumerate(head)}
  for line in barf :
    pieces = line.rstrip().split(sep=',')
    samp = pieces[ idx['Sample_Name'] ]
    sub = samp.split(sep='_')
    year = int(sub[2])
    site = sub[1]
    if year != 2019 : continue
    bcode = pieces[ idx['Barcode_String'] ]
    nN = bcode.count('N')
    nX = bcode.count('X')
    if nX > 1 : continue
    if samp in all_seq_samps[0] :
      mono_bar_samps.add(samp)
      mono_bar_bysite[site] += 1
    else :
      poly_bar_samps.add(samp)
      poly_bar_bysite[site] += 1
    
  merged_seq_samps = all_seq_samps[0] | all_seq_samps[1]
  merged_bar_samps = poly_bar_samps | mono_bar_samps
  all_samps = merged_bar_samps | all_seq_samps[0]
  
  sites = sorted(mono_bar_bysite.keys())
  outf = open(outfile, 'w')
  ndead = 0
  ndead_seq = ndead_bcode = 0
  print('site\tmonogenomic sequence\tpolygenomic sequence\tmono barcode\tpoly barcode', file=outf)
  for site in sites :
    if mono_bar_bysite[site] < 2 and bysite[infiles[0]][site] < 2 :
      ndead += 1
      if mono_bar_bysite[site] == 1 : ndead_bcode += 1
      if bysite[infiles[0]][site] == 1 : ndead_seq += 1
    print(site, bysite[infiles[0]][site], bysite[infiles[1]][site], mono_bar_bysite[site], poly_bar_bysite[site], sep='\t', file=outf)
  with open(prunedfile, 'w') as prunedf :
    for samp in all_seq_samps[0] :
      print(samp, file=prunedf)

  print('dead sites', ndead)
  print('seq samps:', len(merged_seq_samps) - ndead_seq)
  print('barcode samps:', len(merged_bar_samps) - ndead_bcode)
  print('total samps with both:', len( merged_bar_samps & merged_seq_samps ))
  print('mono seq samps:', len(all_seq_samps[0]) - ndead_seq)
  print('mono bar samps:', len(mono_bar_samps) - ndead_bcode)
  print('mono samps with both:', len( mono_bar_samps & all_seq_samps[0] ))
  print('total mono samps:', len( mono_bar_samps | all_seq_samps[0] )  - ndead_bcode)
  print('seq only samps:', len( merged_seq_samps - merged_bar_samps ))
  print('seq only mono samps:', len( all_seq_samps[0] - mono_bar_samps ))
  print('barcode only samps:', len( merged_bar_samps - all_seq_samps[0] ) - ndead_bcode)
  if include_clones : 
    print('total samps used (all barcode + mono seq, with clones):', len( all_samps ) - ndead_bcode)
  else :
    print('total samps used (all barcode + mono seq, no seq clones):', len( all_samps ) - ndead_bcode)

main()
