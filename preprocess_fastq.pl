#!/usr/bin/perl -w -s
## In CELSeq2, read1 contains (in that order) CBC, UMI, polyT, and read2 contains
## the mRNA (starts of fragments). 
## This script puts the CBC and UMI from read1 in front of read2 and prints it
## to stdout.
##
## Original written by Lennart Kester.
##
## The current protocols have an artefact that tends to produces long
## ranges of polyA (and to a lesser extent polyT) Specifying
## e.g. -trim=A=12,T=18 will delete any occurrence of AAAAAAAAAAAA.*$
## TTTTTTTTTTTTTTTTTT.*$ from the read (the quality lines are trimmed in
## the same way). (The numbers suggested in the usage message correspond
## to roughly 0.1% of the actual occurrences in the human transcriptome)
##
## A more sophisticated trimming is provided by the -xytrim=N option.
## This trims any occurrence of polyX-polyY, where X and Y are any combination
## of different homopolymers of length N. E.g. -xytrim=3 would get rid of all
## AAACCC.*, CCCAAA.*, AAAGGG.*, GGGAAA.*, etc. (N=9 is a more reasonable value).
##
## By default CELSeq2 is used, i.e. UMI precedes the cell bar code. Use -protocol=1
## to swap them.

use tools;
use strict;

my $version=getversion($0);
warn "Running $0, version $version\n";

our($fastq, $umi_len, $cbc_len, $trim, $xytrim, $protocol);

if (!($fastq)){
  die "Usage: $0 -fastq=s_R1.fastq[.gz],s_R2.fastq[.gz] -umi_len=6 -cbc_len=8 [-trim=A18,T18] [-xytrim=9] [ -protocol=1 ] | gzip >  s_cbc.fastq.gz ";
}

$protocol=2 if !defined($protocol);

my $regexps ={};

if (defined($trim)) { 
  warn "finding -trim and -xytrim, -trim will be done first" if defined($xytrim);
  my @nucs=split(',', $trim);
  for my $nt (@nucs) { 
    my($nuc, $num)= ($nt =~ /^([ACGT])(\d+)/);
    die "expected string like -trim=A18,T18" unless $nuc && $num;
    $regexps->{$nuc}=$num;
  }

  for my $nuc ( keys  %$regexps ) { 
    my $re = '(' . $nuc x $regexps->{$nuc} . ".*)";
    $regexps->{$nuc}= qr/$re/;
  }
}

if(defined($xytrim)) { 
  use Math::Combinatorics;
  use Regexp::Optimizer;

  my $ndiff=2;
  my @combs=combine($ndiff, qw(A C T G));
  @combs =  map { permute @$_ } @combs;
  
  my $o=Regexp::Optimizer->new;
  
  for my $comb (@combs) { 
    my $name=join("_",@$comb)."_";
    my $quant="{$xytrim,}"; 
    my $re=(join($quant, @$comb)). "$quant";
    $regexps->{$name}=  $o->optimize(qr/($re.*)/i);
  }
}

my $ntrimmed={};
my $ntrimmedtotal={};
my @all=sort keys %$regexps;
my @regexpids=( grep(/^[ATCG]$/, @all ), grep(/_/, @all )); # first -trim, then -xytrim

for my $rid (keys @regexpids) {              # rid=regexp-id
  $ntrimmed->{$rid}=0;
  $ntrimmedtotal->{$rid}=0;
}

die "no -umi_len specified" unless $umi_len > 0; # length of the UMI
die "no -cbc_len specified" unless $cbc_len > 0; # length of the cell bar code

my $prefix_len = $cbc_len + $umi_len;

my @fastq = split(/\,/,$fastq);

# open fastq file
my $cat = "cat ";
$cat = "zcat " if $fastq =~ /\.gz/;

my($IN1, $IN2);
open($IN1, "$cat $fastq[0] |") || die "$fastq[0]: $!";
open($IN2, "$cat $fastq[1] |") || die "$fastq[1]: $!";

my ($line1, $line2, $bar);
my $trimmedlen={};

my (@lines1, @lines2);

READ:
while( not eof $IN1 and not eof $IN2) {
  for(my $i=0; $i<4;$i++) {             # 4 lines at a time
    $lines1[$i] = <$IN1>;
    $lines2[$i] = <$IN2>; 
  }

### id line:
  chomp($lines2[0]);
  my($id, $rest)=split(' ',$lines2[0]);

### sequence line:
  $bar = substr($lines1[1], 0, $prefix_len);
  my $umi=substr($bar,0, $umi_len);
  my $cbc=substr($bar, $umi_len, $cbc_len);
  if ($protocol == 1) { 
    $umi=substr($bar,0, $cbc_len);
    $cbc=substr($bar, $cbc_len, $umi_len);
  }

  $lines2[0] = "$id:cbc=$cbc:umi=$umi $rest\n";

### do trimming, if any:
  my $line2=$lines2[1];
  for my $rid (@regexpids) { 
    if( $line2 =~ $regexps->{$rid} ) { 
      my $newlen=length($line2) - length($1);
      $trimmedlen->{$rid}=$newlen;    # remember for the qual line
      $line2= substr($line2,0, $newlen);
      $ntrimmed->{$rid}++;
      $ntrimmedtotal->{$rid} += length($1)
    }
  }
  $lines2[1]=$line2;

### line with Phred qualities:
  $line2=$lines2[3];
  for my $rid (@regexpids) {               # trim qual line if seqline was
    if(exists($trimmedlen->{$rid})) { 
      $line2= substr($line2,0, $trimmedlen->{$rid});
    }
  }
  $lines2[3]=$line2;

  print  join("", @lines2);
  $trimmedlen={};
}                                       # READ

close $IN1 || die "$fastq[0]: $!";
close $IN2 || die "$fastq[1]: $!";

for my $rid (@regexpids) { 
  warn "trimmed $ntrimmed->{$rid} poly${rid}'s from the reads (totalling $ntrimmedtotal->{$rid} nucleotides)\n"
      if exists($ntrimmed->{$rid}) && $ntrimmed->{$rid} > 0;
}
