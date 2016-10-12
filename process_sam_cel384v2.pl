#!/usr/bin/perl -w -s
use List::Util 'sum';
use Carp;

if ( !($sam && $barfile && $rb_len && $cbc_len) ) { 
  die "usage: $0 -sam=sam.sam -barfile=barfile.csv -rb_len=UMILENGTH -cbc_len=CBCLENGTH";
}

@cells = ();
%bar = ();
%tc = ();

# open barcode file, read all 384 barcodes
$log = $sam;
$log   =~ s/(\.)\w+$/\.log/;
open(LOG, ">", $log) || die "$log: $!";
print LOG "Process sam started\n";

open(IN,"<",$barfile) || die "$barfile:$!";
while(<IN>){   # lines look like ^10 \t GTCGTTCC$ Better use strings for barcode ids!!
  chomp;
  my($id,$barcode) = split(/\t/,$_);
  $bar{$barcode} = $id;                 # e.g. $bar{'GTCGTTCC'}=> 10 . 
  push(@cells, $barcode);
}
close(IN);

$nreads = 0;

# read through sam file create a hash with all genes and cells and extract mapped reads into the hash

open(IN,"<",$sam) || die "$sam: $!";

SAMLINE:
while( <IN> ) {
  chomp $_;
  $line1 = $_;

  if (substr($line1,1,2) eq "SQ" ){     # PL: keep count
    next SAMLINE;
  }

  if (substr($line1,1,2) eq "PG" ){
    print LOG "$line1";                 # PL:should check if it contains bwa
    next SAMLINE;
  }

  @r1 = split(/\t/,$line1);
  my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,@rest)=@r1;
  my $RBC=shift(@rest);   # first is BC:Z: tag, contains barcode (== UMI + CBC)
  my $UMI = substr($RBC,5,$rb_len);
  my $cbc = substr($RBC,(5+$rb_len), $cbc_len); # cell barcode

  $NM = "NA";
  $XA = "NA";
  $X0 = 0;

  foreach $el (@rest){
    ($dum,$dum,$NM) = split(/\:/,$el) if ($el =~ /^NM\:/); # NM: number of mismatches
    ($dum,$dum,$XA) = split(/\:/,$el) if ($el =~ /^XA\:/); # XA: number of alternative hits (chr,pos,CIGAR,NM;)+
    ($dum,$dum,$X0) = split(/\:/,$el) if ($el =~ /^X0\:/); # X0: number of best hits
  }
  ## $X0 = number of locations to which the read maps
  if ($X0 > 0){
    $tot_map_reads++;                   # PL: wrong ...
    ## $tot_map_reads += $X0; # PL
  }
  if ($X0 == 1){
    $tot_map_reads_u++;
  }
  if (exists $bar{$cbc}){
    if ($X0 == 1 && $FLAG != 16){       # flag==16: read is reverse strand, i.e. doesn't map properly
      $tc{$RNAME}{$cbc}{$UMI}++;
    }
  }
  
  $nreads++;
  print LOG int($nreads/1000000) . " million reads processed\n" if ($nreads % 1000000 == 0 );
}                                       # SAMLINE
close(IN);
close(LOG);

my $bn = 4 ** $rb_len;

# create a gene array to loop through while writing the data frame;

@genes = keys %tc;

$coutt = $sam;
$coutb = $sam;
$coutc = $sam;
$sout = $sam;
$trc = 0;
$coutt   =~ s/(\.)\w+$/\.coutt\.csv/;
$coutb   =~ s/(\.)\w+$/\.coutb\.csv/;
$coutc   =~ s/(\.)\w+$/\.coutc\.csv/;

open(OUTT, "> $coutt") || die "$coutt: $!";
open(OUTB, "> $coutb") || die "$coutb: $!";
open(OUTC, "> $coutc") || die "$coutc: $!";

print OUTB "GENEID";
print OUTC "GENEID";
print OUTT "GENEID";

foreach $cbc (@cells){
  print OUTB "\t".$bar{$cbc};
  print OUTC "\t".$bar{$cbc};
  print OUTT "\t".$bar{$cbc};
}	
print OUTB "\n";
print OUTC "\n";
print OUTT "\n";	

GENE:
foreach $gene (sort @genes){
  print OUTB $gene;
  print OUTT $gene;
  print OUTC $gene;
CELL:
  foreach $cbc (@cells){
    my $n = 0;                             # distinct UMIs for this gene+cell
    my $rc = 0;                            # total reads for this gene+cell
  UMI:
    foreach $umi (keys %{$tc{$gene}{$cbc}}){
      next UMI if ($umi =~ /N/i);        # @@@ keep count of these, for stat purposes
      my $reads=$tc{$gene}{$cbc}{$umi};
      $n += ($reads > 0);
      $rc = $rc + $reads; # total reads for this gene+cell
    }                                   # UMI
    $trc += $rc; 
    $n = $n - 0.5 if ($n == $bn); # saturation correction @@@ again: keep count of this
    my $txpts=(-log(1 - ($n/$bn))*$bn);

    print OUTB "\t$n";
    print OUTC "\t$rc";
    print OUTT "\t$txpts";
  }                                     # CELL
  print OUTT "\n";
  print OUTB "\n";
  print OUTC "\n";
}                                       # GENE

$sout   =~ s/(\.)\w+$/\.sout/;
open (SOUT, "> $sout") || die "$sout: $!";
print SOUT "number of reads: $nreads\n";
print SOUT "number of mapped reads: $tot_map_reads\n";
print SOUT "fraction of reads mapped: ".$tot_map_reads/$nreads."\n";
print SOUT "number of uniquely mapped reads: $tot_map_reads_u\n";
print SOUT "fraction of reads mapped uniquely: ".$tot_map_reads_u/$nreads."\n";
print SOUT "number of mapped reads with valid barcode: $trc\n";
print SOUT "fraction of reads mapped with valid barcode: ".$trc/$nreads."\n";

close SOUT;
