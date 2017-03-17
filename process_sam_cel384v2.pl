#!/usr/bin/perl -w
## script to demultiplex CELSeq2 single-cell RNA seq data, do the bookkeeping and convert reads to txpt counts
## original writtten by Dominic Grün and Lennart Kester. Heavily modified by Philip Lijnzaad <plijnzaad@gmail.çom>.

use strict;

use Carp;
use File::Basename;
use FileHandle;
use Getopt::Long;

use tools;

use mismatch;

my $version=getversion($0);
warn "Running $0, version $version\nwith arguments:\n  @ARGV\n";

our ($barfile, $umi_len, $cbc_len, $allow_mm, $prefix, $ref, $help);

my $usage = "
Usage: $0 --barcodefile barcodes.csv --umi_len UMILENGTH --cbc_len CBCLENGTH   [--prefix name ] [ -allow_mm=1 ] file.bam [ file2.bam ...] 

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

### There used to be code to deal with UMIs containing N's, but they are extremely rare and not worth rescuing given
### the complexity of it. They were removed after version v15  (14 march 2017, commit 8fcc58d11a20b0f2cbdc978edfa3742c4e32a900)

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

my @cbcs = sort byletterandnumber (keys %$barcodes); # e.g. 'AGCGTT'
my @wells = map { $barcodes->{$_} } @cbcs;           # e.g. 'B6'

my $mismatch_REs=undef;

if ($allow_mm) { 
  $mismatch_REs = mismatch::convert2mismatchREs(barcodes=>$barcodes_mixedcase, allowed_mismatches =>$allow_mm);
}
$barcodes_mixedcase=undef;              # not used in remainder, delete to avoid confusion


my $nreads = 0;                         # everything
my $ninvalidUMI=0;
my $nantisense=0;
my $nunisense=0;

my $nvalid=0;
my $nignored=0;

## mismatched cell barcodes:
my $nmmCBC=0;
my $nmapped_mmCBC=0;
my $nrescued_mmCBC=0;

my $nmapped=0;
my $nunimapped=0;
my $nrefgenes=0;
my $nERCCs=0;

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

open(SATURATION, "> $saturation") || die "$saturation: $!";
my @headers=qw(reads nmapped nvalid genes umis txpts);
print SATURATION "#" . join("\t", @headers) . "\n";
print SATURATION join("\t", (1) x int(@headers) )."\n";

## well saturation: open three files
my $wellsat_prefix    = "$prefix-wellsat";
my $wellsat_files={};

for my $type ( qw(genes umis)  ) {
  my $file="$wellsat_prefix-$type.txt";
  my $fh = FileHandle->new("> $file") || die "$file: $!";
  $wellsat_files->{$type}=$fh;
  print $fh "reads\t" . join("\t", @wells) . "\n";
}

my $sample_every = 10_000;
my $genes_seen={};                      # cumulative
my $umis_seen={};                       # cumulative
my $wellwise_seen={};                   # keys are those of %$wellsat_files

sub umi_correction { 
  my($n, $maxumis)=@_;
  return 'NA' if $maxumis <= 1;
  $n= $n - 0.5 if $n >= $maxumis;
  sprintf("%.2f", -$maxumis*log(1 - ($n/$maxumis)));
}

HEADERLINE:
while(<IN> ) {
  ## loop is shifted one readline so final printing of SATURATION is not skipped
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

  if ($umi =~ /N/i) { # very rare (< 0.1%), not worth rescuing. See v15 for code that did.
    $ninvalidUMI++;
    $nreads++;
    next READ;
  }

  my $X0 = 0;
  my $dum = 'NA';
  foreach my $el (@rest){
    # ($dum,$dum,$NM) = split(":",$el) if ($el =~ /^NM\:/); # NM: number of mismatches
    # ($dum,$dum,$XA) = split(":",$el) if ($el =~ /^XA\:/); # XA: number of alternative hits (chr,pos,CIGAR,NM;)+
    ($dum,$dum,$X0) = split(":",$el) if ($el =~ /^X0\:/); # X0: number of best hits (bwa-specific!)
  }
  
  $nmapped += ($X0 > 0); # $X0 is number of (equally) optimal locations to which the read maps (one was chosen at random)
  $nunimapped += ($X0 == 1);
  my $antisense= !!($FLAG & 16);
  $nantisense +=  $antisense;
  $nunisense += ( ($X0 == 1) && !$antisense);

  if (! exists $barcodes->{$cbc} && $allow_mm) { 
    $cbc=mismatch::rescue($cbc, $mismatch_REs);      # gives back the barcode without mismatches (if it can be found)
    $nrescued_mmCBC += defined($cbc);
  } 
  
  ## count only reads with valid barcode, uniquely mapping in the sense orientation:
  if ($cbc && exists $barcodes->{$cbc}) {
    if ($X0 == 1 && ! $antisense){ 
      $nvalid++;
      $tc->{$RNAME}{$cbc}{$umi}++; 

      $genes_seen->{$RNAME}++;          ## unless $RNAME =~ /^ERCC-/ (slower)
      $umis_seen->{$RNAME.$umi}++;

      $wellwise_seen->{'genes'}{$cbc}{$RNAME}++;
      $wellwise_seen->{'umis'}{$cbc}{$RNAME.$umi}++;
    } else {
      $nignored++;
      $tc->{'#IGNORED'}{$cbc}{$umi} ++;
      ## keep track of some (non-exclusive!) subsets of this
      $tc->{'#unmapped'}{$cbc}{$umi} += ($X0 == 0 );
      $tc->{'#multimapped'}{$cbc}{$umi} += ($X0 > 1 );
      $tc->{'#antisense'}{$cbc}{$umi} += $antisense; # (may overlap with multimappers)
    }
  } else { 
    $nmmCBC++;
    $nmapped_mmCBC += ($X0 > 0);
  } 
  $nreads++;
} continue { 
  if ($nreads % $sample_every == 0 || eof(IN) ) { 
    my $g=int(keys(%$genes_seen));      # includes those 
    my $u=int(keys(%$umis_seen));       #                 with unknown CBC
    print SATURATION  join("\t", 
                           ($nreads, 
                            $nmapped, # includes non-unique and/or antisense and/or invalid CBC
                            $nvalid, # only the valid reads (those used in downstream analysis)
                            $g, 
                            $u,
                            umi_correction($u,$maxumis*$g))) . "\n";
    
    ## note: for the wellwise counts, only print $nreads; get the rest from the 
    ## overall saturation counts
    my $fh;

    my @genecounts = map { int keys %{$wellwise_seen->{'genes'}{$_}} } @cbcs;
    $fh=$wellsat_files->{'genes'};
    print $fh "$nreads\t" .  join("\t", @genecounts) .  "\n";

    my @umicounts = map { int keys %{$wellwise_seen->{'umis'}{$_}} } @cbcs;
    $fh=$wellsat_files->{'umis'};
    print $fh "$nreads\t" .  join("\t", @umicounts) .  "\n";

  }
  warn int($nreads/1_000_000) . " million reads processed\n" if $nreads % 1_000_000 == 0;
  last READ if eof(IN);
  $_ = <IN>;
}                                        # READ

close(IN) || die "$cmd: $!";
close(SATURATION) || die "$saturation: $!";
foreach my $type (keys %$wellsat_files ) {
  my $fh= $wellsat_files->{$type};
  $fh->close() || die "Well saturation file for type $type :$!";
}

my $coutt   = "$prefix.coutt.csv";
my $coutb   = "$prefix.coutb.csv";
my $coutc   = "$prefix.coutc.csv";
my $sout    = "$prefix.sout";

open(OUTT, "> $coutt") || die "$coutt: $!";
open(OUTB, "> $coutb") || die "$coutb: $!";
open(OUTC, "> $coutc") || die "$coutc: $!";
open (SOUT, "> $sout") || die "$sout: $!";

print OUTT "GENEID\t".join("\t", @wells)."\n";
print OUTB "GENEID\t".join("\t", @wells)."\n";
print OUTC "GENEID\t".join("\t", @wells)."\n";

## gather read counts, umi counts and transcript counts
my $trc = 0;
my $nsaturated_umis=0;

GENE:
foreach my $gene (sort keys %$tc) {
  print OUTB $gene;
  print OUTT $gene;
  print OUTC $gene;
WELL:
  foreach my $cbc (@cbcs) {
    my $n = 0;                          # distinct UMIs for this gene+well
    my $rc = 0;                         # total reads for this gene+well
    my $umihash=$tc->{$gene}{$cbc};
    my @umis = keys %{$umihash};

  UMI:
    foreach my $umi (@umis) {
      my $reads=$tc->{$gene}{$cbc}{$umi};
      $n += ($reads > 0);
      $rc += $reads; # total valid (=uniquely sense-mapped) reads for this gene+cell
    }                                   # UMI
    $trc += $rc unless $gene =~ /^#/;
    $nsaturated_umis += ($n == $maxumis);
    $n = $n - 0.5 if ($n == $maxumis);  # saturation correction
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

sub stat_format { 
  my($part, $total)=@_;
  sprintf("%s / %s   = %.1f %%\n", commafy($part), commafy($total), 100*$part/$total);
}

$ref = "unknown" unless $ref;

print SOUT "reference transcriptome: $ref\n";
print SOUT "        contains $nERCCs ERCCs and " . commafy($nrefgenes) . " genes\n";
print SOUT "invalid UMI (all other numbers refer to valid UMIs): " , stat_format($ninvalidUMI, $nreads);
print SOUT "number of mapped reads: " , stat_format($nmapped, $nreads);
print SOUT "uniquely mapping reads: ", stat_format($nunimapped, $nreads);
print SOUT "antisense mapping reads: ", stat_format($nantisense, $nreads);
print SOUT "uniquely sense mapping reads: ", stat_format($nunisense, $nreads);
print SOUT "uniquely sense mapping with valid CBC: " , stat_format($trc, $nreads);
warn "*** trc ($trc) should be == nvalid ($nvalid) *** \n" if $trc != $nvalid;
print SOUT "rescued mismatched CBC: " , stat_format($nrescued_mmCBC, $nreads);
print SOUT "unknown CBC: " , stat_format($nmmCBC, $nreads );
print SOUT "mapped read, but unknown CBC: " , stat_format($nmapped_mmCBC, $nunisense);
print SOUT "total reads = unique&valid + ignored + unknown CBC + invalidUMI:\n" 
    .     sprintf("%d = %d + %d + %d + %d\n", $nreads,$trc, $nignored, $nmmCBC, $ninvalidUMI);
$nreads /= 100;
print SOUT "%% unique&valid + ignored + mmCBC + invalidUMI:\n" 
    .     sprintf("100%% = %.1f + %.1f + %.1f + %.1f\n", 
                  $trc/$nreads, $nignored/$nreads, $nmmCBC/$nreads, $ninvalidUMI/$nreads);
## note: should be slightly less than 100%: we're threw out the invalid UMIs
print SOUT "gene+well combinations that used up all umis: $nsaturated_umis\n";
close SOUT;
