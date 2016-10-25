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

use strict;

our($fastq, $rb_len, $cbc_len, $trim, $xytrim);

if (!($fastq)){
  die "Usage: $0 -fastq=s_R1.fastq[.gz],s_R2.fastq[.gz] -rb_len=6  -cbc_len=8 [-trim=A18,T18] [-xytrim=9] | gzip >  s_cbc.fastq.gz ";
}

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

die "no -rb_len specified" unless $rb_len > 0; # length of the UMI
die "no -cbc_len specified" unless $cbc_len > 0; # length of the cell bar code

my $prefix_len = $cbc_len + $rb_len;
my $barcode_quality='F';                # i.e. 37
## if the quality is too low, bwa makes the BC:Z:<barcodesequence> all lowercase,
## and it is not mapped anyway due to -B N flag.

my @fastq = split(/\,/,$fastq);

# open fastq file
my $cat = "cat ";
$cat = "zcat " if $fastq =~ /\.gz/;

my($IN1, $IN2);
open($IN1, "$cat $fastq[0] |") || die "$fastq[0]: $!";
open($IN2, "$cat $fastq[1] |") || die "$fastq[1]: $!";

my $i = 0; 
my ($line1, $line2, $bar);
my $trimmedlen={};

LINE:
while( not eof $IN1 and not eof $IN2) {
  $line1 = <$IN1>;
  $line2 = <$IN2>; 
  if ($i == 0){                   # id-line
    print  $line2;
    $i++;
    next LINE;
  }
  if ($i == 1){                   # sequence line
    $bar = substr($line1, 0, $prefix_len);
    chomp($line2);
    # do trimming, if any
    for my $rid (@regexpids) { 
      if( $line2 =~ $regexps->{$rid} ) { 
        my $newlen=length($line2) - length($1);
        $trimmedlen->{$rid}=$newlen;    # remember for the qual line
        $line2= substr($line2,0, $newlen);
        $ntrimmed->{$rid}++;
        $ntrimmedtotal->{$rid} += length($1)
      }
    }
    print  "$bar$line2\n";
    $i++;
    next LINE;
  }
  if ($i == 2){                   # the '+'-line
    print $line2;
    $i++;
    next LINE;
  }
  if ($i == 3){                   # line with Phred qualities
    my $qual= $barcode_quality x $prefix_len;
    chomp($line2);
    for my $rid (@regexpids) {               # trim qual line if seqline was
      if(exists($trimmedlen->{$rid})) { 
        $line2= substr($line2,0, $trimmedlen->{$rid});
      }
    }
    print  "$qual$line2\n";
    $i = 0;
    $trimmedlen={};
  }
}                                       # LINE

close $IN1 || die "$fastq[0]: $!";
close $IN2 || die "$fastq[1]: $!";

for my $rid (@regexpids) { 
  warn "trimmed $ntrimmed->{$rid} poly${rid}'s from the reads (totalling $ntrimmedtotal->{$rid} nucleotides)\n"
      if exists($ntrimmed->{$rid}) && $ntrimmed->{$rid} > 0;
}
