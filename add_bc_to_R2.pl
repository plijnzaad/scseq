#!/usr/bin/perl -w -s
## In CELSeq2, read1 contains (in order) CBC, UMI, polyT, and read2 contains
## the mRNA (starts of fragments). 
## This script puts the CBC and UMI from read1 in front of read2 and prints it
## to stdout.

if (!($fastq)){
  die "usage: -fastq=s_R1.fastq,s_R2.fastq -rb_len=6  -cbc_len=8 [ -A=10 -C=12 -T=18 -G=14 ] > s_cbc.fastq ";
}

die "no -rb_len specified" unless $rb_len > 0; # length of the UMI
die "no -cbc_len specified" unless $cbc_len > 0; # length of the cell bar code

my $prefix_len = $cbc_len + $rb_len;
my $barcode_quality='F';                # i.e. 37
## if the quality is too low, bwa makes the BC:Z:<barcodesequence> all lowercase,
## and it is not mapped anyway due to -B N flag.

my $A_regexp= 'A' x $opt_A;
my $C_regexp= 'C' x $opt_C;
my $G_regexp= 'G' x $opt_G;
my $T_regexp= 'T' x $opt_T;


@fastq = split(/\,/,$fastq);

# open fastq file
my $cat = "cat ";
$cat = "zcat " if $fastq =~ /\.gz/;

open($IN1, "$cat $fastq[0] |") || die "$fastq[0]: $!";
open($IN2, "$cat $fastq[1] |") || die "$fastq[1]: $!";

$i = 0; 

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
    my $qual= $barcode_quality  x $prefix_len;
    print  "$qual$line2";
    $i = 0;
  }
}                                       # LINE

close $IN1 || die "$fastq[0]: $!";
close $IN2 || die "$fastq[1]: $!";


