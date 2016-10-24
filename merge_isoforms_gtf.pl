#!/usr/bin/perl -s -w
## This script merges all isoforms into one (often non-existing)
## transcript: all exons are strung together into one gene model (which
## may have non-existing exon juctions, and will have existing junctions missing)
## The results is a  gtf file with names and coordinates; to create a fasta file
## out of this, use  gffread.

## Note: in the Hg38 RefSeq there are genes with several copies that all have the
## same name (like RNU1-1) (EnsEMBL does this properly, e.g. RN1-1 and
## RNU1-3). The result is that all copies are strung together into one
## supertranscript, *even if they lie on different strands* !!

use tools;

if (!@ARGV && (!$in || !$cl || !$out)) {
  die "usage: -in=INPUT.gtf -cl=clusters_list.tsv -out=OUTPUT.gtf\n";
}

## clusters_list.tsv comes from the -m option of ucsc2gtf.pl, is a table (reverse)mapping gene_id to
## transcript_ids

%stop  = ();
%strand = ();
%chr=();
open(IN,"<",$cl) || die "$cl: $!";
while(<IN>){
    chomp;
    ($cl_name,$e) = split(/\t/);
    @me = split(/\|/,$e);
    foreach $m (@me){
	$clust{$m} = $cl_name;
    }
}
close(IN);

die "$cl: no clusters read " unless int(keys %clust);

open(IN,"<",$in) || die "$in: $!";
while(<IN>){
    chomp;
    next if $_ =~ /^\#/;
    @F = split(/\t/);
    next if $F[2] eq "transcript";
    $info = $F[$#F];
    $info =~ s/(\"|\;)//g;
    ($dum,$dum,$dum,$tid) = split(/\s+/,$info);
    $id = $clust{$tid};
    if (!exists($stop{$id}{$F[3]})){
	$stop{$id}{$F[3]}  = $F[4];
	$strand{$id}       = $F[6];
	$chr{$id}          = $F[0];
    }else{
	$stop{$id}{$F[3]} = $F[4] if $F[4] > $stop{$id}{$F[3]};
    }
}
close(IN);

open(OUT,">",$out);
foreach $id (sort keys %stop){
    @start    = sort {$a <=> $b} keys %{$stop{$id}};
    $tmp_stop = -1;
    @exon     = ();

    foreach $i (0..$#start){
	$s = $start[$i];
	if ($i == 0){
	    $tmp_min = $s;
	    $tmp_max = $stop{$id}{$s};
	    $tmp_start = $s;
	    $tmp_stop = $stop{$id}{$s};
	}else{
	    $tmp_min = min($tmp_min, $s);
	    $tmp_max = max($tmp_max,$stop{$id}{$s});
	    if ($s > $tmp_stop ){
		push(@exon,join("\t",($chr{$id},"merge","exon",$tmp_start,$tmp_stop,0,$strand{$id},".","gene_id \"".$id."\"; transcript_id \"".$id."\"; exon \"".$i."\";"))."\n");
		$tmp_start = $s;
	    }
	    $tmp_stop  = max( $tmp_stop, $stop{$id}{$s});
	}
    }
    push(@exon,join("\t",($chr{$id},"merge","exon",$tmp_start,$tmp_stop,0,$strand{$id},".","gene_id \"".$id."\"; transcript_id \"".$id."\"; exon \"".($#start + 1)."\";"))."\n");
    print OUT join("\t",($chr{$id},"merge","transcript",$tmp_min,$tmp_max,0,$strand{$id},".","gene_id \"".$id."\"; transcript_id \"".$id."\""))."\n";
    foreach $e (@exon){
	print OUT $e;
    }
}
close(OUT);
