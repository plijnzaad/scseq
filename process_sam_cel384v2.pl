#!/usr/bin/perl -w -s
use tools;

use strict;

use Carp;
use File::Basename;

use mismatch;

my ($script,$path) = fileparse($0);
warn "Running $0\n";

our ($sam, $barfile, $rb_len, $cbc_len, $allow_mm);

if ( !($sam && $barfile && $rb_len && $cbc_len) ) { 
  die "Usage: $script -sam=sam.sam -barfile=barfile.csv -rb_len=UMILENGTH -cbc_len=CBCLENGTH [ -allow_mm=1] ";
}

my $barcodes_mixedcase = mismatch::readbarcodes_mixedcase($barfile); ## eg. $h->{'AGCGtT') => 'M3'
my $barcodes = mismatch::mixedcase2upper($barcodes_mixedcase);     ## e.g. $h->{'AGCGTT') => 'M3'

sub bycode {                            # sort the barcodes by their ids
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
$barcodes_mixedcase=undef;

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
    warn "$script: found: $_";                 # PL:should check if it contains bwa
    next SAMLINE;
  }

  my @r1 = split("\t",$_);
  my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,@rest)=@r1;

  my $RBC=shift(@rest);   # first is BC:Z: tag, contains barcode (== UMI + CBC)

  if (!defined($RBC) || $RBC !~ /BC:Z/) { 
    die "expected a BC:Z tag with UMI+cell barcode ($sam line $.)
Is this a sam file from bwa with input from add_bc_to_R2.pl output?";
  }
  my $UMI = substr($RBC,5,$rb_len);
  my $cbc = substr($RBC,(5+$rb_len), $cbc_len); # cell barcode

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
    $cbc=rescue($cbc, $barcodes, $mismatch_REs);      # gives back the barcode without mismatches (if it can be found)
    $nrescued_invalidCBC += defined($cbc);
  }
  if (exists $barcodes->{$cbc}){
    if ($X0 == 1 && $FLAG != 16){ # flag==16: read is reverse strand, i.e. doesn't map properly
      $tc->{$RNAME}{$cbc}{$UMI}++; # only uniquely mapping reads are counted
      # note: invalid umi's are filtered out later!
    } else {
      $nignored++;
      $tc->{'#IGNORED'}{$cbc}{$UMI} ++;
      $tc->{'#unmapped'}{$cbc}{$UMI} += ($X0 == 0 );
      $tc->{'#multimapped'}{$cbc}{$UMI} += ($X0 > 1 );
      $tc->{'#reverse'}{$cbc}{$UMI} += ($FLAG != 16);
    }
  } else { 
    $ninvalidCBC++;
    $nmapped_invalidCBC += ($X0 > 0);
  } 
  $nreads++;
  warn int($nreads/1000000) . " million reads processed\n" if ($nreads % 1000000 == 0 );
}                                       # SAMLINE
close(IN) || die "$cat $sam: $!";

my $bn = 4 ** $rb_len;

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
      $rc += $reads; # total valid (=uniquely mapped) reads for this gene+cell
    }                                   # UMI
    $trc += $rc unless $gene =~ /^#/;      # virtual genes such as '#unmapped';
    $n = $n - 0.5 if ($n == $bn); # saturation correction PL: @@@ keep count of this?
    my $txpts=(-log(1 - ($n/$bn))*$bn); # binomial correction

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
