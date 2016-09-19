#!/usr/bin/perl -w -s

### this is very very slow, and very badly coded (meaningless names,
### globals all over the place). Better download the CCDS transcriptome

use tools;

if (scalar @ARGV == 1){
    die "usage: -in=INPUT.gtf -ref=GENOME.fa -f=ID field name   > transcripts.fa 2> log.txt \n" if $ARGV[0] eq "help";
}
$f = "transcript_id" if !$f;

$genome = {};
warn "reading genome from $ref ...\n";
$genome=fasta2hash($ref);
$j = 0;
$flag  = 0;
%first = ();

## open(TERM, ">> /dev/tty") or die "/dev/tty: $!";

print STDERR "reading $in ...";
open(IN,"<",$in) || die "$in: $!";
while(<IN>){
  chomp;
  @F = split(/\t/);
  next if $F[2] =~ /transcript/;
  if ( !exists($genome->{$F[0]}) ){
    $ns{$F[0]} = 1;
    next;
  }
  @G = split(/\s/,$F[$#F]);
  $id = get_id($_,$f);
  if (!exists($first{$id})){
    if ($flag){
      $seq = get_sequence(\%exons, $genome, $region);
      $seq = revcompl($seq) if $str eq "-";
      print ">".$ID."\n".$seq."\n";
    }
    $j ++;
    print STDERR $j."\r" if $j % 100 ==0; 
    $flag       = 1;
    $first{$id} = 1;
    $ID         = $id;
    $region     = $F[0];
    $str        = $F[6];
    %exons      = ();
  }
  $exons{$F[3]} = $F[4];
}
print STDERR "\ndone\n";

### also finish the very last sequence:
$seq = get_sequence(\%exons,$genome,$region);
$seq = revcompl($seq) if $str eq "-";
print ">".$id."\n".$seq."\n" if exists($genome->{$region});
close(IN);

foreach (sort keys %ns ){
  print STDERR $_."\n";
}

print STDERR "done\n";
exit 0;

### ------------------------------------------------------------------------


sub get_sequence {
  $qhash = shift; # exons
  $rhash = shift; # genome
  $key   = shift; # region
  $flag  = 1;
  my $seq="";
  foreach $k (sort {$a<=>$b} keys %$qhash){
    if ($flag){
      $seq = substr($rhash->{$key}, $k-1, $qhash->{$k} - $k + 1);
    }else{
      $seq = $seq.substr($rhash->{$key}, $k-1, $qhash->{$k} - $k + 1);
    }
    $flag = 0;
  }
  return $seq;
}

sub get_id {
  $x = shift;
  $f = shift;
  chomp($x);
  %h = ();
  @H = split(/\t/,$_);
  $H[$#H] =~ s/[\";]//g;
  @G = split(/\s/,$H[$#H]);

  $i = 0;
  while ( $i < $#G ){
    $key = $G[$i];
    $i++;
    $value = $G[$i];
    $h{$key} = $value;
  }
  return($h{$f});
}
