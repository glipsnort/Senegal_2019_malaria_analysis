#!/usr/bin/perl
use warnings;
use strict;

# The sole function of this script is to detect whether the input file has been compressed or not, 
#  and feed the data in the appropriate way to the script that actually does something.

my $vcf_file = "data/all.vcf.gz";
my $out_file = "output/all_het.txt";
my $min_depth = 5
print "$vcf_file $out_file\n";
if (!-e $vcf_file) {die "File $vcf_file does not exist\n";}

my $com;
# Check if vcf is gzipped
if ($vcf_file =~ m/\.gz$/) {
  $com = "zcat < $vcf_file | ./find_het.pl $out_file $min_depth";
}
else {
  $com = "cat $vcf_file | ./find_het.pl $out_file $min_depth";
}
my $resp = `$com`;
print $resp;
