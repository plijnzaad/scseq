#!/usr/bin/perl -w -s

use tools;

### converts UCSC table file (mysql table with headers
### '#bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds')
### to a GTF file (-out argument), and mapping from gene -> isoform transcript ids
### (in format: gene.'__'. chromosome \t refseqid '__' chromosome '__' 1
### The addition of chromosome name is needed for alternate locus groups (i.e. alternate partial assemblies; typically
### for part of the other half of a diploid genome), as they often have the same genes mapping there.
### Looks like the --utr=3,5 only outputs the 3 or 5' UTRs, in weird format: 
### ---->>>>------>>>>------====  (--- : untransl exon, >>>: intron, ===: translated exon)
### |<- - - - - - - - - - ->| is transcript coordinates (should have been called 5UTR or 3UTR)
###                   |<- ->| is the untranslated part of the first exon (unclear what to call this)

### See also genePredToGtf and 
### http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format

warn "### NOTE: the GRCh38 RefSeq table has duplicate gene names for e.g. RNU1-1; basically this table is corrupt, 
### leading to some huge transcripts (with exons from both strands !?!?) ... UCSC has been notified
" ;

if (@ARGV == 0 && (!$in || !$out || !$m)) {
  die "usage: -in=INPUT.ucsc_format -out=OUTPUT.gtf -m=gene2isoforms.tsv -utr=0,3,5\n";
}

$utr = 0 if !$utr;
open(OUT,">",$out);
open(OUT2,">",$m);
open(IN,"<",$in) or die "$in: $!";
while(<IN>){
  chomp;
  next if /^\#/;
  ($dum,$name,$chr,$str,$start,$end,$cdsStart,$cdsEnd,$dum,$exstart,$exend,$dum,$name2)=split(/\t/);
  if($start <= 0 ) {
    print STDERR "transcript start <= 0, adjusted to 1. Line $.: $_\n";
    $start=1;
    $exstart =~ s/0,/1,/;
  }
  @S = split(/\,/,$exstart);            # exon starts
  @E = split(/\,/,$exend);
  $sid = $name2."__".$chr;              # gene_id
  $t   = $name."__".$chr;               # refseqid, i.e. transcript name
  $count{$t}++;
  $tid = $name."__".$chr."__".$count{$t}; # transcript_id
  push(@{$n{$sid}},$tid);                 # mappings of gene_id
  if ( $utr == 0 ){
    print OUT join("\t",($chr,"UCSC","transcript",$start,$end,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\";"))."\n";
    for $i (0..$#S){
      print OUT join("\t",($chr,"UCSC","exon",$S[$i],$E[$i],0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\"; exon \"".($i + 1)."\";"))."\n";
    }
  }else{
    if ( ( ( $utr == 5 && $str eq "+" ) || ( $utr == 3 && $str eq "-" ) )  && $start < $cdsStart ){
      print OUT join("\t",($chr,"UCSC","transcript",$start,$cdsStart - 1,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\";"))."\n";
      for $i (0..$#S){
	if ( $S[$i] < $cdsStart ){
	  if ( $E[$i] < $cdsStart ){ $e = $E[$i] }else{ $e = $cdsStart - 1};
	  print OUT join("\t",($chr,"UCSC","exon",$S[$i],$e,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\"; exon \"".($i + 1)."\";"))."\n";
	}
      }
    }
    if ( ( ( $utr == 5 && $str eq "-" ) || ( $utr == 3 && $str eq "+" ) ) && $end > $cdsEnd ){
      print OUT join("\t",($chr,"UCSC","transcript",$cdsEnd + 1,$end,0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\";"))."\n";
      for $i (0..$#S){
	if ( $E[$i] > $cdsEnd ){
	  if ( $S[$i] > $cdsEnd ){ $s = $S[$i] }else{ $s = $cdsEnd + 1};
	  print OUT join("\t",($chr,"UCSC","exon",$s,$E[$i],0,$str,".", "gene_id \"".$name2."\"; transcript_id \"".$tid."\"; exon \"".($i + 1)."\";"))."\n";
	}
      }
    }
  }
}
close(IN);
close(OUT);

foreach ( sort keys %n ){
  print OUT2 $_."\t".join("|",@{$n{$_}})."\n";
}
close(OUT2);
