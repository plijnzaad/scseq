#!/usr/bin/perl -w -s
# add barcode
if (!($fastq)){
    die "usage: -fastq=s_X_1.fastq,s_X_2.fastq -rb_len=6"
}


@fastq = split(/\,/,$fastq);


# create out file name
$out = $fastq[1];
$out =~ s/(\.)\w+$/\_cbc.fastq/;
# set random barcodelength
die "no -rb_len specified" unless $rb_len > 0;
# open fastq file
open($IN1, "<", $fastq[0]) || die "$fastq[0]: $!";
open($IN2, "<", $fastq[1]) || die "$fastq[1]: $!";
# set line counter to 0
$i = 0;

# read through fastq file
LINE:
while( not eof $IN1 and not eof $IN2){
	$line1 = <$IN1>;
	$line2 = <$IN2>; 
	chomp $line1;
	chomp $line2;	
# extract index from first line of fastq file	
	if ($i == 0){
		print  $line2."\n";
		$i++;
		next;
	}
# print index in second line of fastq file
	if ($i == 1){
		$bar = substr($line1,0,8+$rb_len);
		print  $bar.$line2."\n";
		$i++;
		next;
	}
	if ($i == 2){
		print  $line2."\n";
		$i++;
		next;
	}
# print extra quality score characters 
	if ($i == 3){
		print  "F" x (8+$rb_len).$line2."\n";
		$i = 0;
	}
}                                       # LINE

close $IN1;
close $IN2;
close OUT;


