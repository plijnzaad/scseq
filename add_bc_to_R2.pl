#!/usr/bin/perl -w -s
## In CELSeq2, read1 contains (in order) CBC, UMI, polyT, and read2 contains
## the mRNA (starts of fragments). 
## This script puts the CBC and UMI from read1 in front of read2, and gives it
## 

if (!($fastq)){
  die "usage: -fastq=s_R1.fastq,s_R2.fastq -rb_len=6  [ -cbc_len=8 ] > s_cbc.fastq ";
}

die "no -rb_len specified" unless $rb_len > 0; # length of the UMI
$cbc_len = 8 if !$cbc_len;

my $prefix_len = $cbc_len + $rb_len;
my $barcode_quality='!';                # i.e. 0 
# (This used to be 'F', i.e. 37, but we don't want this to map ...)

@fastq = split(/\,/,$fastq);

# open fastq file
open($IN1, "<", $fastq[0]) || die "$fastq[0]: $!";
open($IN2, "<", $fastq[1]) || die "$fastq[1]: $!";

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

close $IN1;
close $IN2;


