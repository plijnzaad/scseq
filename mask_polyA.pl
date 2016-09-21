#!/usr/bin/perl -w

## replace any run of A's longer than n (default 10) with as many N's

use strict;

use Getopt::Long;
use Bio::SeqIO;

my($in, $out, $n, $help);
$n=10;

my $usage="Usage: mask_polyA.pl --in input.fa --out output.fa  -n number of As (default: 10)\n" ;

GetOptions('i|in=s' => \$in,
           'o|out=s'=> \$out,
           'n=i'    => \$n,
           'h|help' => \$help
    ) || die $usage;

die $usage if $help || !$in  || !$out ;

my $regexp = "([Aa]{$n,})";

my $input  = Bio::SeqIO->new(-file => $in , '-format' => 'fasta');
my $output = Bio::SeqIO->new(-file=> "> $out", '-format' => 'fasta');

while ( my $seq=$input->next_seq()  ) { 
  my $s=$seq->seq();

  ## $s= join("", map { $_ = 'N' x length($_) if /$regexp/; $_}  split(/$regexp/, $s)); # slightly slower

  while ( $s =~ /$regexp/g ){
    my $x = $&;
    $x =~ s/[Aa]/N/g;
    $s = substr($s,0,pos($s) - length($x)).$x.substr($s,pos($s),length($s) - pos($s));
  }

  $seq->seq($s);
  $output->write_seq($seq);
}

