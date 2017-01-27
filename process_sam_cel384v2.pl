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

our ($barfile, $umi_len, $cbc_len, $allow_mm, $prefix, $ref, $help);

my $usage = "
Usage: $0 --barcodefile barcodes.csv --umi_len UMILENGTH --cbc_len CBCLENGTH  [ -allow_mm=1 ] [--prefix name ] file.bam [ file2.bam ...]

Arguments: 

--barcodefile=FILE  File with cell bar codes (format: id \\t sequence)
--umi_len=N         Length of the UMIs

file.bam ...        Name(s) of the bam file(s), typically from several lanes. If
                    more than one they are assumed to have the same header.
Options:
--help              This message
--allow_mm=N        How many mismatches to allow in the cell barcodes (default: 0)

--prefix=name       Prefix for the four output files: NAME.coutt.csv, NAME.coutb.csv, NAME.coutc.csv and NAME.sout
                    Default is the name of the first bam file without extension and lane number.
--ref=name          Name of the reference genome (only for logging)
";

die $usage unless GetOptions('barcodefile=s'=> \$barfile,
                             'umi_len=i'=> \$umi_len,
                             'cbc_len=i'=> \$cbc_len,
                             'allow_mm=i'=> \$allow_mm,
                             'prefix=s' => \$prefix,
                             'ref=s' => \$ref,
                             'help|h' => \$help);
my @bams=@ARGV;

if ( $help  || !(@bams && $barfile && $umi_len && $cbc_len) ) { 
  die $usage;
}

my $barcodes_mixedcase = mismatch::readbarcodes($barfile); ## eg. $h->{'AGCGtT') => 'M3'
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
my $nignored=0;

my $ninvalidCBC=0;
my $nmapped_invalidCBC=0;
my $nrescued_invalidCBC=0;

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

my $nrefgenes=0;
my $nERCCs=0;

SAMLINE:
while( <IN> ) {
  chomp $_;

  if (substr($_,1,2) eq "SQ" ){
    my (@tag, $name)=split("\t", $_);
    if ($name =~ /ERCC/) { $nERCCs++; } else { $nrefgenes ++; }
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
  die "$0: could not find cbc= or umi= in id $QNAME of the files @bams " unless $cbc &&  $umi;

  my $X0 = 0;
  my $dum = 'NA';

  foreach my $el (@rest){
    # ($dum,$dum,$NM) = split(":",$el) if ($el =~ /^NM\:/); # NM: number of mismatches
    # ($dum,$dum,$XA) = split(":",$el) if ($el =~ /^XA\:/); # XA: number of alternative hits (chr,pos,CIGAR,NM;)+
    ($dum,$dum,$X0) = split(":",$el) if ($el =~ /^X0\:/); # X0: number of best hits
  }
  
  $nmapped += ($X0 > 0); # $X0 is number of locations to which the read maps
  $nunimapped += ($X0 == 1);
  my $reverse=!!($FLAG & 16);
  $nreverse +=  $reverse;

  if (! exists $barcodes->{$cbc} && $allow_mm) { 
    $cbc=mismatch::rescue($cbc, $mismatch_REs);      # gives back the barcode without mismatches (if it can be found)
    $nrescued_invalidCBC += defined($cbc);
  } 

  ## count only reads with valid barcode, uniquely mapping in the sense orientation:
  if ($cbc && exists $barcodes->{$cbc}){
    if ($X0 == 1 && ! $reverse){ 
      $tc->{$RNAME}{$cbc}{$umi}++; 
      # note: invalid umi's are filtered out later!
    } else {
      $nignored++;
      $tc->{'#IGNORED'}{$cbc}{$umi} ++;
      ## keep track of some subsets of this
      $tc->{'#unmapped'}{$cbc}{$umi} += ($X0 == 0 );
      $tc->{'#multimapped'}{$cbc}{$umi} += ($X0 > 1 );
      $tc->{'#reverse'}{$cbc}{$umi} += $reverse; # (may overlap with multimappers)
    }
  } else { 
    $ninvalidCBC++;
    $nmapped_invalidCBC += ($X0 > 0);
  } 
  $nreads++;
  warn int($nreads/1000000) . " million reads processed\n" if ($nreads % 1000000 == 0 );
}                                       # SAMLINE
close(IN) || die "$cmd: $!";

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
  UMI:
    foreach my $umi (keys %{$tc->{$gene}{$cbc}}) {
      if ($umi =~ /N/i  && $gene !~ /#/) { 
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
  }                                     # CELL
  print OUTT "\n";
  print OUTB "\n";
  print OUTC "\n";
}                                       # GENE

sub stat_format { 
  my($part, $total)=@_;
  sprintf("%s / %s   = %.1f %%\n", commafy($part), commafy($total), 100*$part/$total);
}

$ref = "unknown" unless $ref;

print SOUT "reference transcriptome: $ref\n";
print SOUT "        contains $nERCCs ERCCs and " . commafy($nERCCs) . "\n";
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
