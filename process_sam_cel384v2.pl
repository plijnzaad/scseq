#!/usr/bin/perl -w
## script to demultiplex CELSeq2 single-cell RNA seq data, do the bookkeeping and convert reads to txpt counts
## original writtten by Dominic Grün and Lennart Kester
use strict;

use Carp;
use File::Basename;
use Getopt::Long;

use tools;

use mismatch;

my $version=getversion($0);
warn "Running $0, version $version\nwith arguments:\n  @ARGV\n";

our ($barfile, $umi_len, $cbc_len, $allow_mm, $prefix, $ref, $rescue_umis, $help);

my $usage = "
Usage: $0 --barcodefile barcodes.csv --umi_len UMILENGTH --cbc_len CBCLENGTH   [--prefix name ] [ -allow_mm=1 ] [ --rescue_umis ] file.bam [ file2.bam ...] 

Arguments: 

--barcodefile FILE  File with cell bar codes (format: id \\t sequence)
--umi_len N         Length of the UMIs

file.bam ...        Name(s) of the bam file(s), typically from several lanes. If
                    more than one they are assumed to have the same header.
Options:
--help              This message
--allow_mm N        How many mismatches to allow in the cell barcodes (default: 0)

--prefix name       Prefix for the four output files: NAME.coutt.csv, NAME.coutb.csv, NAME.coutc.csv and NAME.sout
                    Default is the name of the first bam file without extension and lane number.
--ref name          Name of the reference genome (only for logging)
";

die $usage unless GetOptions('barcodefile=s'=> \$barfile,
                             'umi_len=i'=> \$umi_len,
                             'cbc_len=i'=> \$cbc_len,
                             'allow_mm=i'=> \$allow_mm,
                             'prefix=s' => \$prefix,
                             'ref=s' => \$ref,
                             'rescue_umis'=> \$rescue_umis,
                             'help|h' => \$help);
my @bams=@ARGV;

if ( $help  || !(@bams && $barfile && $umi_len && $cbc_len) ) { 
  die $usage;
}

my $maxumis = 4 ** $umi_len;

my $barcodes_mixedcase = mismatch::readbarcodes($barfile); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'


sub byletterandnumber { # sort the barcodes by their ids (which may contain prefixes)
  my ($aa,$bb) = ($a,$b);
  $aa=$barcodes->{$aa}; 
  $bb=$barcodes->{$bb}; 
  my($re)= qr/([A-Za-z_]*)([0-9]+)/;
  my ($Sa, $Na) = ($aa =~ $re);
  my ($Sb, $Nb) = ($bb =~ $re);
  ($Sa cmp $Sb) || ($Na <=>  $Nb);
}

my @wells = sort byletterandnumber (keys %$barcodes); 

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

## mismatched cell barcodes:
my $nmmCBC=0;
my $nmapped_mmCBC=0;
my $nrescued_mmCBC=0;

my $nmapped=0;
my $nunimapped=0;

my $tc = {};

# read through sam file create a hash with all genes and cells and extract mapped reads into the hash
my $samtools = "samtools";                    # alternative: sambamba, might be faster

if (!$prefix) { 
  $prefix=$bams[0];
  $prefix =~ s/\.bam//;
  $prefix =~ s/[-_,.]L0*[1-4][-_,.]//g if @bams > 1;
  warn "using prefix $prefix\n";
}

for my $f (@bams) {
  die "$f: $!" unless $f; 
}
my $cmd = "($samtools view -H $bams[0]; " . join("; ", map { "samtools view $_" } @bams) . ")";

open(IN,"$cmd |") || die "$cmd: $!";
my $saturation    = "$prefix-saturation.txt";
my $wellsaturation    = "$prefix-wellsaturation.txt";

open(SATURATION, "> $saturation") || die "$saturation: $!";
my @headers=qw(reads nmapped genes umis txpts);
print SATURATION "#" . join("\t", @headers) . "\n";
print SATURATION join("\t", (1) x int(@headers) )."\n";

open(WELLSATURATION, "> $wellsaturation") || die "$wellsaturation: $!";

my $sample_every = 10_000;
my $genes_seen={};                      # cumulative
my $umis_seen={};                       # cumulative
my $nrefgenes=0;
my $nERCCs=0;

sub umi_correction { 
  my($n, $maxumis)=@_;
  $n= $n - 0.5 if $n >= $maxumis;
  sprintf("%.2f", -$maxumis*log(1 - ($n/$maxumis)));
}

HEADERLINE:
while(<IN> ) {
  ## loop is shifted one readline so final printing of SATURATION is OK
  if ( /^\@SQ/ ){
    my ($tag, $name, @rest)=split("\t", $_);
    if ($name =~ /^ERCC-0/) { $nERCCs++; } else { $nrefgenes++; }
    next HEADERLINE;
  }
  if ( /^\@PG/ ){
    warn "$0: found: $_";                 # PL:should check if it contains bwa
    next HEADERLINE;
  }
  last HEADERLINE unless /^@/;
}

READ:
while(1) { 
  chomp $_;
  my @r1 = split("\t",$_);
  my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,@rest)=@r1;

  my ($cbc,$umi);
  my(@parts)=split(':', $QNAME);
  for my $tag ( @parts ) { 
    $cbc= $1 if $tag =~ /cbc=([A-Z]+)/i;
    $umi= $1 if $tag =~ /umi=([A-Z]+)/i;
  }
  die "$0: could not find cbc= or umi= in id $QNAME of the files @bams " unless $cbc &&  $umi;

  my $X0 = 0;
  my $dum = 'NA';

  foreach my $el (@rest){
    # ($dum,$dum,$NM) = split(":",$el) if ($el =~ /^NM\:/); # NM: number of mismatches
    # ($dum,$dum,$XA) = split(":",$el) if ($el =~ /^XA\:/); # XA: number of alternative hits (chr,pos,CIGAR,NM;)+
    ($dum,$dum,$X0) = split(":",$el) if ($el =~ /^X0\:/); # X0: number of best hits (bwa-specific!)
  }
  
  $nmapped += ($X0 > 0); # $X0 is number of (equally) optimal locations to which the read maps (one was chosen at random)
  $nunimapped += ($X0 == 1);
  my $reverse= !!($FLAG & 16);
  $nreverse +=  $reverse;

  if (! exists $barcodes->{$cbc} && $allow_mm) { 
    $cbc=mismatch::rescue($cbc, $mismatch_REs);      # gives back the barcode without mismatches (if it can be found)
    $nrescued_mmCBC += defined($cbc);
  } 

  ## count only reads with valid barcode, uniquely mapping in the sense orientation:
  if ($cbc && exists $barcodes->{$cbc}){
    if ($X0 == 1 && ! $reverse){ 
      $tc->{$RNAME}{$cbc}{$umi}++; 
      $genes_seen->{$RNAME}++;          ## unless $RNAME =~ /^ERCC-/ (slowish)
      $umis_seen->{$RNAME.$umi}++;
      # note: invalid umi's are filtered out later
    } else {
      $nignored++;
      $tc->{'#IGNORED'}{$cbc}{$umi} ++;
      ## keep track of some (non-exclusive!) subsets of this
      $tc->{'#unmapped'}{$cbc}{$umi} += ($X0 == 0 );
      $tc->{'#multimapped'}{$cbc}{$umi} += ($X0 > 1 );
      $tc->{'#reverse'}{$cbc}{$umi} += $reverse; # (may overlap with multimappers)
    }
  } else { 
    $nmmCBC++;
    $nmapped_mmCBC += ($X0 > 0);
  } 
  $nreads++;
} continue { 
  if ($nreads % $sample_every == 0 || eof(IN) ) { 
    my $g=int(keys(%$genes_seen));
    my $u=int(keys(%$umis_seen));
    print SATURATION  join("\t", 
                     ($nreads, 
                      $nmapped,         # includes non-unique, reverse and invalid CBC
                      $g, 
                      $u,
                      umi_correction($u,$maxumis*$g))) . "\n";
  }
  warn int($nreads/1000_0000) . " million reads processed\n" if $nreads % 1000_000 == 0;
  last READ if eof(IN);
  $_ = <IN>;
}                                        # READ

close(IN) || die "$cmd: $!";
close(SATURATION) || die "$saturation: $!";

my $coutt   = "$prefix.coutt.csv";
my $coutb   = "$prefix.coutb.csv";
my $coutc   = "$prefix.coutc.csv";
my $sout    = "$prefix.sout";

open(OUTT, "> $coutt") || die "$coutt: $!";
open(OUTB, "> $coutb") || die "$coutb: $!";
open(OUTC, "> $coutc") || die "$coutc: $!";
open (SOUT, "> $sout") || die "$sout: $!";

print OUTB "GENEID";
print OUTC "GENEID";
print OUTT "GENEID";

foreach my $cbc (@wells){
  my $id=$barcodes->{$cbc};
  print OUTB "\t$id";
  print OUTC "\t$id";
  print OUTT "\t$id";
}	
print OUTB "\n";
print OUTC "\n";
print OUTT "\n";	

## gather read counts, umi counts and transcript counts
my $trc = 0;

GENE:
foreach my $gene (sort keys %$tc) {
  print OUTB $gene;
  print OUTT $gene;
  print OUTC $gene;
WELL:
  foreach my $cbc (@wells) {
    my $n = 0;                          # distinct UMIs for this gene+cell
    my $rc = 0;                         # total reads for this gene+cell
    my $umihash=$tc->{$gene}{$cbc};
    my @umis = keys %{$umihash};

    if ( $rescue_umis && $gene !~ /#/ ) { # preprocessing to rescue UMI's containing N's
      my @Ns=grep(/N/i, @umis);
      if (@Ns) { 
        my @noNs=grep(! /N/i, @umis);

        my $oldus=join(',', keys %$umihash);
        my $h=cleanup_umis(\@noNs, \@Ns, $umihash );
        my $newus=join(',', keys %$umihash);
        warn "UMIn: $gene\t$cbc\t$oldus=>$newus\tdisc:$h->{discarded} resc:$h->{rescued}\n";
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
            if $rescue_umis;    # should have become 'X' or disappeared altogether
        $ninvalidUMI ++ ;
        next UMI;
      }
      my $reads=$tc->{$gene}{$cbc}{$umi};
      $n += ($reads > 0);
      $rc += $reads; # total valid (=uniquely sense-mapped) reads for this gene+cell
    }                                   # UMI
    $trc += $rc unless $gene =~ /^#/;
    $n = $n - 0.5 if ($n == $maxumis); # saturation correction PL: @@@ keep count of this?
    my $txpts = $n;                      # used only for '#IGNORED' etc. @@@fix this
    $txpts = -log(1 - ($n/$maxumis)) * $maxumis unless ($gene =~ /^#/ ); # binomial correction

    print OUTB "\t$n";
    print OUTC "\t$rc";
    print OUTT "\t$txpts";
  }                                     # WELL
  print OUTT "\n";
  print OUTB "\n";
  print OUTC "\n";
}                                       # GENE

sub cleanup_umis { 
  ## costly. Deletes invalid UMIs from the hash, converts useable ones and returns number of deleted and rescued UMIs
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
    $re="^$re\$";
    my @hits=grep( /$re/, @$noNs);

    if ( @hits ) {                 
      # e.g. /ACT./ ~ ACTG but could represent ACTA => 2 umis
      delete $umihash->{$N};
      $ndiscarded++;
      next UMI;
    }
    ## In case you can't tolerate any N's: 
    my $new=$N;
    $new =~ s/[Nn]/X/g;   
    $umihash->{$new} = $umihash->{$N};
    delete $umihash->{$N};
    $nrescued++;
  }                                     # UMI

  { discarded=>$ndiscarded, 
    rescued=>$nrescued,
  };
}                                       # cleanup_umis

sub stat_format { 
  my($part, $total)=@_;
  sprintf("%s / %s   = %.1f %%\n", commafy($part), commafy($total), 100*$part/$total);
}

$ref = "unknown" unless $ref;

print SOUT "reference transcriptome: $ref\n";
print SOUT "        contains $nERCCs ERCCs and " . commafy($nrefgenes) . " genes\n";
print SOUT "number of mapped reads: " , stat_format($nmapped, $nreads);
print SOUT "uniquely mapping reads: ", stat_format($nunimapped, $nreads);
print SOUT "uniquely with valid cbc and umi: " , stat_format($trc, $nreads);
print SOUT "valid barcode, invalid UMI: " , stat_format($ninvalidUMI, $nreads);
print SOUT "rescued mismatched CBC: " , stat_format($nrescued_mmCBC, $nreads);
print SOUT "unknown CBC: " , stat_format($nmmCBC, $nreads );
print SOUT "mapped read, but  CBC: " , stat_format($nmapped_mmCBC, $nunimapped);
print SOUT "total reads = unique&valid + ignored + mismatched CBC + invalidUMI:\n" 
    .     sprintf("%d = %d + %d + %d + %d\n", $nreads,$trc, $nignored, $nmmCBC, $ninvalidUMI);

$nreads /= 100;
print SOUT "%% unique&valid + ignored + mmCBC + invalidUMI:\n" 
    .     sprintf("100%% = %.1f + %.1f + %.1f + %.1f\n", 
                  $trc/$nreads, $nignored/$nreads, $nmmCBC/$nreads, $ninvalidUMI/$nreads);

close SOUT;
