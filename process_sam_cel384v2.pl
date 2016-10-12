#!/usr/bin/perl -w -s
use Carp;
use strict;
use tools;

our ($sam, $barfile, $rb_len, $cbc_len);

if ( !($sam && $barfile && $rb_len && $cbc_len) ) { 
  die "usage: $0 -sam=sam.sam -barfile=barfile.csv -rb_len=UMILENGTH -cbc_len=CBCLENGTH";
}

my @cells = ();
my %bar = ();
my %tc = ();

# open barcode file, read all 384 barcodes
my $log = $sam;
$log   =~ s/(\.)\w+$/\.log/;
open(LOG, ">", $log) || die "$log: $!";
print LOG "Process sam started\n";

open(IN,"<",$barfile) || die "$barfile:$!";
while(<IN>){   # lines look like ^10 \t GTCGTTCC$ Better use strings for barcode ids!!
  chomp;
  my($id,$barcode) = split("\t",$_);
  $bar{$barcode} = $id;                 # e.g. $bar{'GTCGTTCC'}=> 10 . 
  push(@cells, $barcode);
}
close(IN);

my $nreads = 0;
my $tot_map_reads=0;
my $tot_map_reads_u=0;
my $noBCtag=0;

# read through sam file create a hash with all genes and cells and extract mapped reads into the hash

open(IN,"<",$sam) || die "$sam: $!";

SAMLINE:
while( <IN> ) {
  chomp $_;

  if (substr($_,1,2) eq "SQ" ){     # PL: keep count
    next SAMLINE;
  }

  if (substr($_,1,2) eq "PG" ){
    print LOG "$_";                 # PL:should check if it contains bwa
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
  
  $tot_map_reads += ($X0 > 0); # $X0 is number of locations to which the read maps
  $tot_map_reads_u += ($X0 == 1);

  if (exists $bar{$cbc}){
    if ($X0 == 1 && $FLAG != 16){       # flag==16: read is reverse strand, i.e. doesn't map properly
      $tc{$RNAME}{$cbc}{$UMI}++;
      ## @@: only uniquely mapping reads are counted !
    }
  }  
  $nreads++;
  print LOG int($nreads/1000000) . " million reads processed\n" if ($nreads % 1000000 == 0 );
}                                       # SAMLINE
close(IN);
close(LOG);

my $bn = 4 ** $rb_len;

# create a gene array to loop through while writing the data frame;

my @genes = keys %tc;

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
  print OUTB "\t".$bar{$cbc};
  print OUTC "\t".$bar{$cbc};
  print OUTT "\t".$bar{$cbc};
}	
print OUTB "\n";
print OUTC "\n";
print OUTT "\n";	

my $trc = 0;

GENE:
foreach my $gene (sort @genes){
  print OUTB $gene;
  print OUTT $gene;
  print OUTC $gene;
CELL:
  foreach my $cbc (@cells) {
    my $n = 0;                             # distinct UMIs for this gene+cell
    my $rc = 0;                            # total reads for this gene+cell
  UMI:
    foreach my $umi (keys %{$tc{$gene}{$cbc}}) {
      next UMI if ($umi =~ /N/i);
      my $reads=$tc{$gene}{$cbc}{$umi};
      $n += ($reads > 0);
      $rc = $rc + $reads; # total reads for this gene+cell
    }                                   # UMI
    $trc += $rc; 
    $n = $n - 0.5 if ($n == $bn); # saturation correction @@@ keep count of this?
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

print SOUT "number of reads: " , stat_format($tot_map_reads, $nreads);
print SOUT "uniquely mapped reads: ", stat_format($tot_map_reads_u, $nreads);
print SOUT "uniquely mapped reads with valid barcode: " , stat_format($trc, $tot_map_reads_u);
## (original version used wrong ratio!)

close SOUT;
