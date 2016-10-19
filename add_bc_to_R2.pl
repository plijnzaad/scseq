#!/usr/bin/perl -w -s
## In CELSeq2, read1 contains (in order) CBC, UMI, polyT, and read2 contains
## the mRNA (starts of fragments). 
## This script puts the CBC and UMI from read1 in front of read2 and prints it
## to stdout.

## use strict;

our($fastq, $rb_len, $cbc_len);

if (!($fastq)){
  die "usage: -fastq=s_R1.fastq,s_R2.fastq -rb_len=6  -cbc_len=8 [ -A=10 -C=12 -T=18 -G=14 ] > s_cbc.fastq ";
}

## dirty hack to avoid having to use epxlicit vars
my $polynucs ={};

$regexps->{'A'}=$A if defined($A);
$regexps->{'C'}=$C if defined($C);
$regexps->{'G'}=$G if defined($G);
$regexps->{'T'}=$T if defined($T);
$ncleanedup={};

my @nucs = keys  %$regexps;
for my $nuc (@nucs) { 
  $re = '(' . $nuc x $regexps->{$nuc} . ".*)";
  $regexps->{$nuc}= qr/$re/;
  $ncleanedup->{$nuc}=0;
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
    for my $nuc (@nucs) { 
      if( $line2 =~ $regexps->{$nuc} ) { 
        $newlen=length($line2) - length($1);
        $trimmedlen->{$nuc}=$newlen;    # remember for the qual line
        $line2= substr($line2,0, $newlen);
        $ncleanedup->{$nuc}++;
      }
    }
    print  "$bar$line2";
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
    for my $nuc (@nucs) {               # trim qual line if seqline was
      if(exists($trimmedlen->{$nuc})) { 
        $line2= substr($line2,0, $trimmedlen->{$nuc});
      }
    }
    print  "$qual$line2";
    $i = 0;
    $trimmedlen={};
  }
}                                       # LINE

close $IN1 || die "$fastq[0]: $!";
close $IN2 || die "$fastq[1]: $!";


