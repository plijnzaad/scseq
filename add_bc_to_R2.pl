#!/usr/bin/perl -w -s
## In CELSeq2, read1 contains (in that order) CBC, UMI, polyT, and read2 contains
## the mRNA (starts of fragments). 
## This script puts the CBC and UMI from read1 in front of read2 and prints it
## to stdout.
##
## The current protocols have an artefact that tends to produces long
## ranges of polyA (and to a lesser extent polyT) Specifying
## e.g. -trim=A=12,T=18 will delete any occurrence of AAAAAAAAAAAA.*$
## TTTTTTTTTTTTTTTTTT.*$ from the read (the quality lines are trimmed in
## the same way). (The numbers suggested in the usage message correspond
## to roughly 0.1% of the actual occurrences in the human transcriptome)
## 

use strict;

our($fastq, $rb_len, $cbc_len, $trim);

if (!($fastq)){
  die "Usage: $0 -fastq=s_R1.fastq[.gz],s_R2.fastq[.gz] -rb_len=6  -cbc_len=8 [ -trim=A18,T18 ] | gzip >  s_cbc.fastq.gz ";
}

## dirty hack to avoid having to use epxlicit vars
my $regexps ={};

if (defined($trim)) { 
  my @nucs=split(',', $trim);
  for my $nt (@nucs) { 
    my($nuc, $num)= ($nt =~ /^([ACGT])(\d+)/);
    die "expected string like -trim=A18,T18" unless $nuc && $num;
    $regexps->{$nuc}=$num;
  }
}

my $ntrimmed={};
my $ntrimmedtotal={};

my @nucs = keys  %$regexps;
for my $nuc (@nucs) { 
  my $re = '(' . $nuc x $regexps->{$nuc} . ".*)";
  $regexps->{$nuc}= qr/$re/;
  $ntrimmed->{$nuc}=0;
  $ntrimmedtotal->{$nuc}=0;
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
    for my $nuc (@nucs) { 
      if( $line2 =~ $regexps->{$nuc} ) { 
        my $newlen=length($line2) - length($1);
        $trimmedlen->{$nuc}=$newlen;    # remember for the qual line
        $line2= substr($line2,0, $newlen);
        $ntrimmed->{$nuc}++;
        $ntrimmedtotal->{$nuc} += length($1)
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
    for my $nuc (@nucs) {               # trim qual line if seqline was
      if(exists($trimmedlen->{$nuc})) { 
        $line2= substr($line2,0, $trimmedlen->{$nuc});
      }
    }
    print  "$qual$line2\n";
    $i = 0;
    $trimmedlen={};
  }
}                                       # LINE

close $IN1 || die "$fastq[0]: $!";
close $IN2 || die "$fastq[1]: $!";

for my $nuc (@nucs) { 
  warn "trimmed $ntrimmed->{$nuc} poly${nuc}'s from the reads (totalling $ntrimmedtotal->{$nuc} nucleotides)\n"
      if exists($ntrimmed->{$nuc}) && $ntrimmed->{$nuc} > 0;
}
