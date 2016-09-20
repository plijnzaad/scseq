#!/usr/bin/perl -w -s

## replace any run of A's longer than n (default 10) with as many N's

use tools;

if (scalar @ARGV == 1){
    die "usage: -in=input.fa -n=number of As (default: 10) -out=output.fa\n" if $ARGV[0] eq "help";
}

$n = 10 if !$n;
$seqs = fasta2hash($in);
$m = "(A|a){".$n.",}";

open(OUT,">",$out) or die "$out:$!";
foreach $k (keys %$seqs){
  $s = $seqs->{$k};
  while ( $s =~ /$m/g ){
    $x = $&;
    $x =~ s/[Aa]/N/g;
    $s = substr($s,0,pos($s) - length($x)).$x.substr($s,pos($s),length($s) - pos($s));
  }
  print OUT ">".$k."\n".$s."\n";
}
close(OUT);

