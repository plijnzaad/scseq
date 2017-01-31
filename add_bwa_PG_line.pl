#!/usr/bin/perl -w
#

use strict;

my $usage = "Usage (e.g.) :  bwa samse -n N reference.fa output.sai input.fastq.gz  |  add_bwa_PG_line.pl  bwa aln -B 0 -t 4  -q 40 -n 0.02 -k 2 -l 200 reference.fa input.fastq.gz |  samtools view -h -b -   > file.bam";

die $usage unless @ARGV >= 3 && $ARGV[0] eq 'bwa';
my @args = @ARGV; 
@ARGV=();

while(<>) { 
  if( /^\@PG.*bwa.*sam[ps]e/) { 
    my($tag, $id, $pn, $vn, @rest)=split("\t");
    print join("\t", ($tag, $id, $pn, $vn, "CL:".join(" ", @args)))."\n";
    print;
  } else { 
    print;
  }
}


