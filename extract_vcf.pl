#!/usr/bin/perl
use warnings;
use strict;

# The sole function of this script is to detect whether the input file has been compressed or not, 
#  and feed the data in the appropriate way to the script that actually does something.

my $base = shift;
my $min_depth = shift;
if (!defined $base) {die "Usage: extract_vcf.pl <file base name (e.g. 2019)> [min read depth]\n";}
my $vcf_file = "data/$base.vcf.gz";
my $out_file = "seq/$base" . '_seq.txt';
my $samp_file = "output/$base" . "_good_mono_samples.txt";
my $allele_file = "seq/$base" . "_allele.txt";

if (!-e $vcf_file) {die "File $vcf_file does not exist\n";}
if (!defined $min_depth) {$min_depth = 0;}

my $com;
# Check if vcf is gzipped
if ($vcf_file =~ m/\.gz$/) {
  $com = "zcat < $vcf_file | ./vcf_genotypes.pl $out_file $samp_file $allele_file $min_depth";
}
else {
  $com = "cat $vcf_file | ./vcf_genotypes.pl $out_file $samp_file $allele_file $min_depth";
}
my $resp = `$com`;
print $resp;
