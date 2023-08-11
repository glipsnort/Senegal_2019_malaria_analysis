#!/usr/local/bin/python3

from collections import Counter
import gzip
import re
import sys

def main() :
  kill_indel = False
  # Mininum genotyping call rate to accept variant
  min_call = 0.80
  # Option: minimum number of alt allele copies to keep variant
  #  (1 = kill monomorphic, 2 = kill singletons)
  min_copy = 0
  min_depth = 0
  
  if len(sys.argv) < 2 : sys.exit('Usage: vcf2hmm.ph <input vcf file name> [<files of samples to use>]')
  infile = sys.argv[1]
  seqfile = 'seq/seq.txt'
  allfile = 'seq/all.txt'
  all_samps = True
  good_samps = set()
  if len(sys.argv) == 3 :
    all_samps = False
    sampf = open(sys.argv[2], 'r')
    for line in sampf :
      samp = line.rstrip()
      print(samp)
      good_samps.add(samp)
    
  if re.search(r'\.gz$', infile) :
    inf = gzip.open(infile, 'rt')
  else :
    inf = open(infile, 'r')

  samples = []
  keep_samps = []
  nonchrom = Counter()
  depthsum = Counter()    # by sample
  ndepth = Counter()
  nhet = Counter()
  ncall = Counter()
  nocall = Counter()
  nhomo = 0
  nindel = 0
  with open(seqfile, 'w') as seqf, open(allfile, 'w') as allf :
    for line in inf :
      if re.match(r'\#CHROM', line) :
        samples = line.rstrip().split('\t')[9:]
        for samp in samples :
          if all_samps or samp in good_samps : keep_samps.append(samp)
        print('chrom', 'pos', '\t'.join(keep_samps), file=seqf, sep='\t')
        continue
      elif re.match(r'\#', line) : continue
      
      # data line
      pieces = line.rstrip().split('\t')
      chrom = str(pieces[0])
      # Handle falciparum chrom names
      match = re.search(r'^Pf3D7_(\d+)_v3', chrom)
      if match :
        chrom = int(match.group(1))
      pos = int(pieces[1])
      ref_all = pieces[3]
      alt_alls = pieces[4].split(',')
      kill = False
      if kill_indel and not re.match(r'[ACGT\.]$', ref_all) :
        nindel += 1
        continue
      for alt_all in alt_alls :
        if kill_indel and not re.match(r'[ACGT\.]$', alt_all) :
          nindel += 1
          kill = True
      if kill : continue

      formats = pieces[8].split(':')
      genotypes = pieces[9:]
      outline = '{:d}\t{:d}'.format(chrom, pos)
      nassay = 0

      for isamp, samp in enumerate(samples) :
        if not all_samps and samp not in good_samps : continue
        genotype = genotypes[isamp].split(':')
        call = ''
        for iform, form in enumerate(formats) :
          if form == 'GT' :
            call = genotype[iform]
          elif form == 'DP' :
            if genotype[iform] == '.' : depth = 0
            else : depth = int(genotype[iform])
            depthsum[samp] += depth
            ndepth[samp] += 1
        allele = '-2'
        if re.match(r'\d+', call) and depth > min_depth :
          match = re.match(r'(\d+)\/(\d+)', call)
          g1 = match.group(1)
          g2 = match.group(2)
          ncall[samp] += 1
          if g1 != g2 :
            allele = '-1'
            nhet[samp] += 1
          else :
            allele = g1
        nassay += 1
        outline += ('\t' + allele)
      print(outline, file=seqf)
      print('{:d}\t{:d}'.format(chrom, pos), ref_all, '\t'.join(alt_alls), sep='\t', file=allf)
                           
main()
