#!/usr/bin/perl -w -s
## script to demultiplex CELSeq2 single-cell RNA seq data, do the bookkeeping and convert reads to txtp counts
## original writtten by Dominic Grün and Lennart Kester
use tools;

use strict;

use Carp;
use File::Basename;

use mismatch;

my $version=getversion($0);
warn "Running $0, version $version\n";

our ($sam, $barfile, $umi_len, $cbc_len, $allow_mm, $protocol, $rescue_umi_Ns);

if ( !($sam && $barfile && $umi_len && $cbc_len) ) { 
  die "Usage: $0 -sam=sam.sam -barfile=barfile.csv [-allow_mm=1] -umi_len=UMILENGTH -cbc_len=CBCLENGTH [ -protocol=1 ] [ -rescue_umi_Ns=1 ]";
}

$protocol =2 unless $protocol;

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($barfile); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'

sub bycode {                            # sort the barcodes by their ids (which may contain prefixes)
  my ($aa,$bb) = ($a,$b);
  $aa=$barcodes->{$aa}; $aa =~ s/[A-Za-z_]//g;
  $bb=$barcodes->{$bb}; $bb =~ s/[A-Za-z_]//g; 
  $aa <=> $bb;
}
my @cells = sort bycode (keys %$barcodes); # @cells bar codes sorted by their id's (e.g. c1, c2, ... )

my $mismatch_REs=undef;

if ($allow_mm) { 
  $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, allowed_mismatches =>$allow_mm);
}
$barcodes_mixedcase=undef;              # not used in remainder, delete to avoid confusion


my $nreads = 0;
my $nreverse=0;
my $ninvalidUMI=0;
my $nrescued_invalidUMI=0;

my $nignored=0;

my $ninvalidCBC=0;
my $nmapped_invalidCBC=0;
my $nrescued_invalidCBC=0;

my $nmapped=0;
my $nunimapped=0;

my $tc = {};

# read through sam file create a hash with all genes and cells and extract mapped reads into the hash
my $cat = "cat ";
my $samtools = "samtools";                    # alternative: sambamba, might be faster
$cat = "$samtools view " if $sam =~ /\.bam$/;

open(IN,"$cat $sam |") || die "$sam: $!";

SAMLINE:
while( <IN> ) {
  chomp $_;

  if (substr($_,1,2) eq "SQ" ){
    next SAMLINE;
  }

  if (substr($_,1,2) eq "PG" ){
    warn "$0: found: $_";                 # PL:should check if it contains bwa
    next SAMLINE;
  }

  my @r1 = split("\t",$_);
  my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,@rest)=@r1;

  my ($cbc,$umi);
  my(@parts)=split(':', $QNAME);
  for my $tag ( @parts ) { 
    $cbc= $1 if $tag =~ /cbc=([A-Z]+)/i;
    $umi= $1 if $tag =~ /umi=([A-Z]+)/i;
  }
  die "$0: could not find cbc= or umi= in id $QNAME of file $sam " unless $cbc &&  $umi;

  my $X0 = 0;
  my $dum = 'NA';

  foreach my $el (@rest){
    # ($dum,$dum,$NM) = split(":",$el) if ($el =~ /^NM\:/); # NM: number of mismatches
    # ($dum,$dum,$XA) = split(":",$el) if ($el =~ /^XA\:/); # XA: number of alternative hits (chr,pos,CIGAR,NM;)+
    ($dum,$dum,$X0) = split(":",$el) if ($el =~ /^X0\:/); # X0: number of best hits
  }
  
  $nmapped += ($X0 > 0); # $X0 is number of locations to which the read maps
  $nunimapped += ($X0 == 1);
  $nreverse += ($FLAG == 16);

  if (! exists $barcodes->{$cbc} && $allow_mm) { 
    $cbc=mismatch::rescue($cbc, $mismatch_REs);      # gives back the barcode without mismatches (if it can be found)
    $nrescued_invalidCBC += defined($cbc);
  } 

  ## count only reads with valid barcode, uniquely mapping in the sense orientation:
  if ($cbc && exists $barcodes->{$cbc}){
    if ($X0 == 1 && $FLAG != 16){ # flag==16: read is reverse strand
      $tc->{$RNAME}{$cbc}{$umi}++; 
      # note: invalid umi's are filtered out later!
    } else {
      $nignored++;
      $tc->{'#IGNORED'}{$cbc}{$umi} ++;
      ## keep track of some subsets of this
      $tc->{'#unmapped'}{$cbc}{$umi} += ($X0 == 0 );
      $tc->{'#multimapped'}{$cbc}{$umi} += ($X0 > 1 );
      $tc->{'#reverse'}{$cbc}{$umi} += ($FLAG != 16); # (may overlap with multimappers)
    }
  } else { 
    $ninvalidCBC++;
    $nmapped_invalidCBC += ($X0 > 0);
  } 
  $nreads++;
  warn int($nreads/1000000) . " million reads processed\n" if ($nreads % 1000000 == 0 );
}                                       # SAMLINE
close(IN) || die "$cat $sam: $!";


my $coutt = $sam;
my $coutb = $sam;
my $coutc = $sam;
my $sout = $sam;

$coutt   =~ s/(\.)\w+$/\.coutt\.csv/;
$coutb   =~ s/(\.)\w+$/\.coutb\.csv/;
$coutc   =~ s/(\.)\w+$/\.coutc\.csv/;
$sout   =~ s/(\.)\w+$/\.sout/;

open(OUTT, "> $coutt") || die "$coutt: $!";
open(OUTB, "> $coutb") || die "$coutb: $!";
open(OUTC, "> $coutc") || die "$coutc: $!";

print OUTB "GENEID";
print OUTC "GENEID";
print OUTT "GENEID";

foreach my $cbc (@cells){
  my $id=$barcodes->{$cbc};
  print OUTB "\t$id";
  print OUTC "\t$id";
  print OUTT "\t$id";
}	
print OUTB "\n";
print OUTC "\n";
print OUTT "\n";	

my $trc = 0;

## print read counts, umi counts and transcript counts
my $bn = 4 ** $umi_len;                 # max #(different umis)

GENE:
foreach my $gene (sort keys %$tc) {
  print OUTB $gene;
  print OUTT $gene;
  print OUTC $gene;
CELL:
  foreach my $cbc (@cells) {
    my $n = 0;                             # distinct UMIs for this gene+cell
    my $rc = 0;                            # total reads for this gene+cell
    my $umihash=$tc->{$gene}{$cbc};
    my @umis = keys %{$umihash};

    if ( $rescue_umi_Ns ) {             # preprocessing to rescue UMI's containing N's
      my @Ns=grep(/N/i, @umis);
      if (@Ns) { 
        warn "*** this code has not been tested yet";
        warn "*** gene $gene cell $cbc has several N-containing UMIs, or with a UMI containing several N's:\n"  
            . join('\n***', @umis) if ( @Ns > 1 || grep(/N.*N/i, @Ns) ); # prolly rare
        my @noNs=grep(! /N/i, @umis);

        my $h=cleanup_umis(\@noNs, \@Ns, $umihash );
        $umihash = $h->{umis};
        $tc->{$gene}{$cbc} = $umihash;
        @umis = keys %{$umihash};

        $ninvalidUMI += $h->{discarded};
        $nrescued_invalidUMI += $h->{rescued};
      }
    }

  UMI:
    foreach my $umi (@umis) {
      if ($umi =~ /N/i  && $gene !~ /#/) { 
        confess "gene $gene, cbc $cbc, umi $umi contains N, should not happen when rescueing UMIs" 
            if $rescue_umi_Ns;    # should have become 'X' or disappeared altogether
        $ninvalidUMI ++ ;
        next UMI;
      }
      my $reads=$tc->{$gene}{$cbc}{$umi};
      $n += ($reads > 0);
      $rc += $reads; # total valid (=uniquely sense-mapped) reads for this gene+cell
    }                                   # UMI
    $trc += $rc unless $gene =~ /^#/;
    $n = $n - 0.5 if ($n == $bn); # saturation correction PL: @@@ keep count of this?
    my $txpts = $n;                      # used only for '#IGNORED' etc. @@@fix this
    $txpts = -log(1 - ($n/$bn)) * $bn unless ($gene =~ /^#/ ); # binomial/Poisson correction

    print OUTB "\t$n";
    print OUTC "\t$rc";
    print OUTT "\t$txpts";
  }                                     # CBC/CELL
  print OUTT "\n";
  print OUTB "\n";
  print OUTC "\n";
}                                       # GENE

sub cleanup_umis { 
  ## costly ...
  my ($noNs,$Ns, $umihash)=@_;          # umihash contains read counts per umi
  my @newNs=();
  my($nrescued, $ndiscarded)=(0,0);

 UMI:
  for my $N ( @$Ns ) { 
    if ($umihash->{$N} >1) { 
      # only keep singles, as 2 x ACTN could have come from ACTT and ACTA
      delete $umihash->{$N};
      $ndiscarded++;
      next UMI;
    }
    my $re=$N;
    $re =~ s/[Nn]/./g;
    my @hits=grep( qr/$re/, @$noNs) > 1;

    if ( @hits ) {                 
      # e.g. /ACT./ ~ ACTG but could represent ACTA => 2 umis
      delete $umihash->{$N};
      $ndiscarded++;
      next UMI;
    }
    ## In case you can't tolerate any N's: 
    ## my $new=$N;
    ## $new =~ s/[Nn]/X/g;   
    ## $umihash->{$new} = $umihash->{$N};
    ## delete $umihash->{$N};
    $nrescued++;
  }                                     # UMI

  { umihash => $umihash,
    discarded=>$ndiscarded, 
    rescued=>$nrescued,
  };
}                                       # cleanup_umis

sub stat_format { 
  my($part, $total)=@_;
  sprintf("%s / %s   = %.1f %%\n", commafy($part), commafy($total), 100*$part/$total);
}

open (SOUT, "> $sout") || die "$sout: $!";

print SOUT "number of mapped reads: " , stat_format($nmapped, $nreads);
print SOUT "uniquely mapping reads: ", stat_format($nunimapped, $nreads);
print SOUT "uniquely with valid cbc and umi: " , stat_format($trc, $nreads);
print SOUT "valid barcode, invalid UMI: " , stat_format($ninvalidUMI, $nreads);
print SOUT "rescued invalid CBC: " , stat_format($nrescued_invalidCBC, $nreads);
print SOUT "invalid CBC: " , stat_format($ninvalidCBC, $nreads);
print SOUT "mapped read, but invalid CBC: " , stat_format($nmapped_invalidCBC, $nunimapped);
print SOUT "total reads = unique&valid + ignored + invalidCBC + invalidUMI:\n" 
    .     sprintf("%d = %d + %d + %d + %d\n", $nreads,$trc, $nignored, $ninvalidCBC, $ninvalidUMI);

$nreads /= 100;
print SOUT "%% unique&valid + ignored + invalidCBC + invalidUMI:\n" 
    .     sprintf("100%% = %.1f + %.1f + %.1f + %.1f\n", 
                  $trc/$nreads, $nignored/$nreads, $ninvalidCBC/$nreads, $ninvalidUMI/$nreads);

close SOUT;
