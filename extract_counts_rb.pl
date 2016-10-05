#!/usr/bin/perl -w -s

use tools;

if (! ( $in && $outc && $outb && $outt )) { 
    die "
Usage: -bl=4  # UMI-length              \
       -in=aggr_counts.tsv                  \
       -outc=READ_COUNTS.csv            \
       -outb=BARCODE_COUNTS.cvs         \
       -outt=TRANSCRIPT_COUNTS.csv\n";  
}

$bl = 4 if !$bl;
$bn = 4 ** $bl;
$flag = 0;
@ina = split(/\,/,$in);
foreach $k (@ina){
  open(IN,"<",$k) || die "$k: $!";
  while(<IN>){
    chomp;
    if ( $_ =~ /GENEID/ ){
      @title = split(/\t/);
      if ($_ =~ /CLASS/){
	$flag = 1;
      }
      next;
    }
    if ( $flag ){
      @F = split(/\t/);
      $anno{$F[1]} = $F[0];
      @F = @F[1..$#F];
    }else{
      @F = split(/\t/);
    }    
    next if $F[1] =~ /N/;
    if (!exists($rc{$F[0]})){
      for $i (2..$#F){
	${$rc{$F[0]}}[$i - 2] = 0;
	${$bc{$F[0]}}[$i - 2] = 0;
      }  
    }
    for $i (2..$#F){
	${$rc{$F[0]}}[$i - 2] += $F[$i];
        ${$bc{$F[0]}}[$i - 2] += min(1,$F[$i]) if $F[$i] > 0 && ! exists($seen{$F[0]}{$F[1]}{$i});
	$seen{$F[0]}{$F[1]}{$i} = 1 if $F[$i] > 0;   
      }
   
  }
  close(IN);
}

open(OUTC,">",$outc) || die "$outc:$!";
open(OUTB,">",$outb) || die "$outb:$!";;
open(OUTT,">",$outt) || die "$outt:$!";;
if ( $flag ){
  print OUTC join("\t",(@title[0..1],@title[3..$#title]))."\n";
  print OUTB join("\t",(@title[0..1],@title[3..$#title]))."\n";
  print OUTT join("\t",(@title[0..1],@title[3..$#title]))."\n";
}else{
  print OUTC join("\t",($title[0],@title[2..$#title]))."\n";
  print OUTB join("\t",($title[0],@title[2..$#title]))."\n";
  print OUTT join("\t",($title[0],@title[2..$#title]))."\n";
}
foreach (sort keys %rc){
  print OUTC $anno{$_}."\t" if ( $flag );
  print OUTC $_;
  foreach $n ( @{$rc{$_}} ){
    print OUTC "\t".$n;
  }
  print OUTC "\n";
}
foreach (sort keys %bc){
  print OUTB $anno{$_}."\t" if ( $flag );
  print OUTT $anno{$_}."\t" if ( $flag );
  print OUTB $_;
  print OUTT $_;
  foreach $n ( @{$bc{$_}} ){
    print OUTB "\t".$n;
    if ($n == $bn){ $n = $n - .5; }
    print OUTT "\t".(-log(1 - $n/$bn)*$bn);
  }
  print OUTB "\n";
  print OUTT "\n";
}
close(OUTC);
close(OUTB);
close(OUTT);
