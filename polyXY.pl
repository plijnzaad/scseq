#!/bin/env perl

### find stretches of polyXY 's (e.g. AAAACCCCC but also CCCCCAAAAA, for
### all combinations of two diff nucleotides Needs (on stdin) sequences in
### tab format, i.e.  ID \t ACACACTCAGCGAGCTCCAC...ACACCACTGTT\n
### 
### This script is/was used to find out what a reasonable trimvalue is for the -xytrim option
### in add_bc_to_R2.pl
### (Has since been abandoned, results in miniscule mapping improvements)

use strict;

use Math::Combinatorics;
use Regexp::Optimizer;

my $ndiff=2;
my @combs=combine($ndiff, qw(A C T G));
@combs =  map { permute @$_ } @combs;

my $REs={};

my $o=Regexp::Optimizer->new;

my ($min)=shift;

die "Usage: $0:  ... | $0   N " unless $min > 0;

for my $comb (@combs) { 
   my $name=join("_",@$comb)."_";
  my $quant="{$min,}"; 
  my $RE=(join($quant, @$comb)). "$quant";
  $REs->{$name}=  $o->optimize(qr/$RE/i);
}

my @sids=sort keys %$REs;
my $hits={};

while(<>) {
  my ($id, $seq)=split("\t", $_);
  s/N/A/g;
  for my $sid (@sids) { 
    $hits->{$sid} += ($seq =~ $REs->{$sid});
  }
}

for my $sid (@sids) { 
  my $n=$hits->{$sid};
  $n=0 unless $n;
  print "$sid\t$n\n";
}
