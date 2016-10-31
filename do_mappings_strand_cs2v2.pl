#!/usr/bin/perl -w -s
## driver script for single-cell RNAseq scripts, originally written by Dominic Grün
use tools;
use Carp;
use File::Basename;
use strict;

$,=" ";         # for interpolating arrays inside strings (default anyway)

our($r, $f1, $f2, $out, $bwaparams, $outdir, $t, $ind, $q, $aln_n, $aln_k, $l,
    $BR, $i, $npr, $nsam, $bar, $umi_len, $cbc_len, $trim, $xytrim, $allow_mm, $test, $protocol);

if (!($r && $f1 && $out && $t)){
  confess "Usage:  $0 -r=REFERENCE     \
                 -f1=READ1    \
                 -f2=READ2 (optional)    \
                 -out=OUTPUT_PREFIX    \
                 -bwaparams='...' (optional, overrides the -q -n -k -l options) \
                 -outdir=OUTPUT_DIRECTORY (optional)    \
                 -t=THREADS  (0)  \
                 -ind=is or bwtsw (default: is (<2GB reference)    \
                 -q=PHRED_QUALITY_FOR_TRIMMING (0..62, default 0)   \
                 -aln_n=edit distance  (bwa aln option -n, default 0.04)  \
                 -aln_k=edit distance in seed  (bwa aln option -k, default 2)   \
                 -l=SEED_LENGTH  (bwa -l option)  \
                 -i= 1 or 0 (1 if indexing is required, runs bwa index )    \
                 -npr=0,1,2: 0: map and process; 1: only map ; 2: only process
                 -nsam= 0 or 1 (1: do *not* produce new sam file (calls bwa samse/sampe)    \
                 -bar=cel-seq_barcodes.csv    \
                 -umi_len= length of UMI (default = 6)  \
                 -cbc_len= length of cell barcode (default: 8) \
                 -trim='A12,T=14' (optional, passed to preprocess_fastq.pl for trimming) \
                 -xytrim=10 (optional, passed to preprocess_fastq.pl for trimming) \
                 -allow_mm=N (optional, passed to process_sam_cel384v2, allows N mismatches in cell bar codes) \
                 -test=0 or 1 (latter runs in test mode, doesn't call external programs) 
                 -protocol=2 (2 (celseq2): umi=6, cbc=8; 1 (celseq1): umi=4, cbc=8 \
";
}

# hard coded:
my $n = 100; # maximal number of unpaired reads to output in XA tag
my $N = 100; # maximal number of paired reads to output in XA tag

$q = 0 if !$q; # base quality cutoff for trimming
$l = 200 if !$l; # seed length


$i = 0 if !$i; # create index
$npr  = 0 if !$npr;
$nsam = 0 if !$nsam;
$ind = "is" if !$ind;
$protocol = 2 if !$protocol;

$umi_len = ({1=>4,2=>6})->{$protocol} if !$umi_len;
$cbc_len = 8 if !$cbc_len;

$aln_n = 0.04 if !$aln_n; # edit distance
$aln_k = 2 if !$aln_k; # edit distance in seed

$bwaparams=" -q $q -n $aln_n -k $aln_k -l $l " unless $bwaparams;

$trim="" unless $trim;
$xytrim="" unless $xytrim;

$trim = "-trim=$trim" if $trim;
$xytrim = "-xytrim=$xytrim" if $xytrim;

if ($allow_mm) { 
  $allow_mm="-allow_mm=$allow_mm";
} else { 
  $allow_mm="";
}

$test = 0 if !$test;
if ($outdir){
  makedir($outdir);
  $out = "$outdir/$out";
}

confess "reference genome $r not found" unless -f $r;

if ($i){
  my $str = "bwa index -a $ind $r";
  print $str."\n";
  execute(cmd=>$str, merge=>1) if ($test == 0);
}

my ($name, $path, $ext)=fileparse($f2, ('.fastq', '.fastq.gz'));

my $cbc = "$path/${name}_cbc.fastq";       # output from add_bc_to_R2, input to bwa
my $gzip = " cat ";

if ( $ext =~ /\.gz/ ) { 
  $cbc .= ".gz";
  $gzip = " gzip ";
}

if ( $npr != 2 ) {                       # npr is 0 or 1: do mapping
  if ( -f $cbc  ) { 
    print "*** Seeing file $cbc, not running preprocess_fastq.pl to re-create it\n";
  } else { 
    my $str = "preprocess_fastq.pl -fastq=$f1,$f2 -umi_len=$umi_len -cbc_len=$cbc_len $trim $xytrim | $gzip > $cbc ";
    print $str."\n";
    execute(cmd=>$str, merge=>0) if ($test == 0);
    check_filesize(file=>$cbc, minsize=>1000);
  }

  my $sai = "$outdir/$name.sai"; 

  my $str = "bwa aln -B 0 -t $t $bwaparams $r $cbc > $sai";
  print $str."\n";
  execute(cmd=>$str, merge=>0) if ($test == 0);
  check_filesize(file=>$sai, minsize=>1000);
  
  if ( $nsam == 0 ){
    my $compress = "samtools view -h -b - ";
    $str = "bwa samse -n $n $r $sai $cbc | $compress > $out.bam";
  }
  print $str."\n";
  execute(cmd=>$str) if ($test == 0);
  check_filesize(file=>"$out.bam",minsize=>1000);
}                                       # npr!=2

if ( $npr == 0 || $npr == 2){
  my $str = "process_sam_cel384v2.pl -sam=$out.bam -barfile=$bar -umi_len=$umi_len -cbc_len=$cbc_len $allow_mm";
  print $str."\n";
  execute(cmd=>$str, merge=>1) if ($test == 0);
}                                       # if ( $npr == 0 || $npr == 2)
