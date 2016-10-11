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
while(<IN>){	
  chomp;
  @line = split(/\t/,$_);
  $bar{$line[1]} = $line[0];
  push(@cells, $line[1]);
}
close(IN);

$i = 0;
$x = 1;
$tot_map_reads = 0;
$tot_map_reads_u = 0;

#create variables to store the different components of the sam file for R1;
$QNAME = ();
$FLAG = ();
$RNAME = ();
$POS = ();
$MAPQ = ();
$CIGAR = ();
$MRNM = ();
$MPOS = ();
$ISIZE = ();
$SEQ = ();
$QUAL = ();
$NM = ();
$XA = ();
$X0 = ();
$RBC = ();
$UMI = "NA";
@bar_map = ();
%bar_map_h = ();

# read through sam file create a hash with all genes and cells and extract mapped reads into the hash

open(IN,"<",$sam) || "$sam:$!";
SAMLINE:
while( <IN> ) {
  chomp $_;
  $line1 = $_;

  if (substr($line1,1,2) eq "SQ" ){
    next SAMLINE;
  }
  if (substr($line1,1,2) eq "PG" ){
    print LOG "$line1";
    next SAMLINE;
  }
# here all the information is extracted from the sam-file line;
  @r1 = split(/\t/,$line1);
  ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,$RBC)=@r1;
  $cellbc = substr($RBC,($rb_len)+5,8);
  $UMI = substr($RBC,5,$rb_len);
  $NM = "NA";
  $XA = "NA";
  $X0 = 0;

  foreach $el (@r1){
    ($dum,$dum,$NM) = split(/\:/,$el) if ($el =~ /^NM\:/);
    ($dum,$dum,$XA) = split(/\:/,$el) if ($el =~ /^XA\:/);
    ($dum,$dum,$X0) = split(/\:/,$el) if ($el =~ /^X0\:/);
    
  }
  ## $X0 = number of locations to which the read maps
  if ($X0 > 0){
    $tot_map_reads++;
  }
  if ($X0 == 1){
    $tot_map_reads_u++;
  }
  if (exists $bar{$cellbc}){
    if ($X0 == 1 && $FLAG != 16){
      $tc{$RNAME}{$cellbc}{$UMI}++;
    }
  }
  
  $i++;
  if ($i == $x*1000000){
    print LOG $x." million reads processed\n";
    $x++;
  }
}                                       # SAMLINE
close(IN);
close(LOG);


$bn = 4 ** $rb_len;

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

open(OUTT, ">", $coutt) || die "$coutt: $!";
open(OUTB, ">", $coutb) || die "$coutb: $!";
open(OUTC, ">", $coutc) || die "$coutc: $!";

print OUTB "GENEID";
print OUTC "GENEID";
print OUTT "GENEID";

foreach $cell (@cells){
  print OUTB "\t".$bar{$cell};
  print OUTC "\t".$bar{$cell};
  print OUTT "\t".$bar{$cell};
}	
print OUTB "\n";
print OUTC "\n";
print OUTT "\n";	
foreach $gene (sort @genes){
  print OUTB $gene;
  print OUTT $gene;
  print OUTC $gene;
  foreach $cell (@cells){
    $n = 0;
    $rc = 0;
    foreach $umi (keys %{$tc{$gene}{$cell}}){
      if ($umi =~ /N/i){
        next;
      }
      if ($tc{$gene}{$cell}{$umi} > 0.1){
        $n++;
      }	
      $rc = $rc + $tc{$gene}{$cell}{$umi};
    }
    print OUTB "\t".$n;
    if ($n == $bn) {
      $n = $n - 0.5;
    }
    print OUTC "\t".$rc;
    $trc = $trc + $rc; 
    print OUTT "\t".(-log(1 - ($n/$bn))*$bn);
  }
  print OUTT "\n";
  print OUTB "\n";
  print OUTC "\n";
}

$sout   =~ s/(\.)\w+$/\.sout/;
open (SOUT, ">", $sout);
print SOUT "number of reads: ".$i."\n";
print SOUT "number of mapped reads: ".$tot_map_reads."\n";
print SOUT "fraction of reads mapped: ".$tot_map_reads/$i."\n";
print SOUT "number of uniquely mapped reads: ".$tot_map_reads_u."\n";
print SOUT "fraction of reads mapped uniquely: ".$tot_map_reads_u/$i."\n";
print SOUT "number of mapped reads with valid barcode: ".$trc."\n";
print SOUT "fraction of reads mapped with valid barcode: ".$trc/$i."\n";
close SOUT;
