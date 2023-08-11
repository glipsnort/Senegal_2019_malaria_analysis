#!/usr/local/bin/python3

def main() :
  cover_file = 'output/all_het.txt'
  ibd_file = 'output/all.hmm_fract.txt'
  out_file = 'output/sequenced_clones.txt'
  
  coverage = {}
  outset = set()     # samples to output for removal
  with open(cover_file, 'r') as coverf :
    head = coverf.readline().rstrip().split()
    idx = {col : x for x, col in enumerate(head)}
    for line in coverf :
      pieces = line.rstrip().split()
      samp = pieces[idx['samp']]
      coverage[samp] = float(pieces[idx['depth_gt5x']])

  with open(ibd_file, 'r') as inf :
    head = inf.readline().rstrip().split()
    idx = {col : x for x, col in enumerate(head)}
    for line in inf :
      pieces = line.rstrip().split()
      samp1 = pieces[idx['sample1']]
      samp2 = pieces[idx['sample2']]
      ibd = float(pieces[idx['fract_sites_IBD']])
      if ibd < 0.98 : continue
      p1 = samp1.split('_')
      p2 = samp2.split('_')
      year1 = p1[2]
      year2 = p2[2]
      site1 = p1[1]
      site2 = p2[1]
      if year1 != year2 : continue
      if site1 != site2 : continue
      cover1 = coverage[samp1]
      cover2 = coverage[samp2]
      if samp1 in outset :
        continue
      if samp2 in outset :
        continue
      which = samp1 if cover1 >= cover2 else samp2
      outset.add(which)

  with open(out_file, 'w') as outf :
    for samp in outset :
      print(samp, file=outf)
main()
