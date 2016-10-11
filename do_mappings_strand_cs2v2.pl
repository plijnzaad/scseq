#!/usr/bin/perl -w -s
use tools;
use Carp;
use File::Basename;

$LIST_SEPARATOR=" ";                    # for interpolating arrays inside strings (default anyway)

if (!($r && $f1 && $out && $t)){
  confess "Usage:  $0 -r=REFERENCE     \
                 -f1=READ1    \
                 -f2=READ2 (optional)    \
                 -out=OUTPUT_PREFIX    \
                 -outdir=OUTPUT_DIRECTORY (optional)    \
                 -t=THREADS  (0)  \
                 -ind=is or bwtsw (default: is (<2GB reference)    \
                 -q=PHRED_QUALITY_FOR_TRIMMING (0..62, default 0)   \
                 -aln_n=edit distance  (bwa aln option -n, default 0.04)  \
                 -aln_k=edit distance in seed  (bwa aln option -k, default 2)   \
                 -l=SEED_LENGTH  (bwa -l option)  \
                 -BL=barcode length left read (bwa aln option -B, default is 0)    \
                 -BR=barcode length right read (bwa aln option -B, default 0)   \
                 -i= 1 or 0 (1 if indexing is required, runs bwa index )    \
                 -npr=0,1,2: 0: map and process; 1: only map ; 2: only process
                 -nsam= 0 or 1 (1: do *not* produce new sam file (calls bwa samse/sampe)    \
                 -bar=cel-seq_barcodes.csv    \
                 -rb_len= length of UMI (default = 6)  \
                 -cbc_len= length of cellseq2 barcode (default: 8) \
                 -test=0 or 1 (latter runs in test mode, doesn't call external programs) \
";
}

# hard coded:
$n = 100; # maximal number of unpaired reads to output in XA tag
$N = 100; # maximal number of paired reads to output in XA tag
$q = 0 if !$q; # base quality cutoff for trimming
$l = 200 if !$l; # seed length
$i = 0 if !$i; # create index
$BL = 0 if !$BL;
$BR = 0 if !$BR;
$npr  = 0 if !$npr;
$nsam = 0 if !$nsam;
$ind = "is" if !$ind;
$rb_len = 6 if !$rb_len;
$cbc_len = 8 if !$cbc_len;

$BR = $cbc_len+$rb_len;

$aln_n = 0.04 if !$aln_n; # edit distance
$aln_k = 2 if !$aln_k; # edit distance in seed

$test = 0 if !$test;
if ($outdir){
  makedir($outdir);
  $out = "$outdir/$out";
}

confess "reference genome $r not found" unless -f $r;

if ($i){
  $str = "bwa index -a $ind $r";
  print $str."\n";
  execute(cmd=>$str, merge=>1) if ($test == 0);
}

my ($name, $path, $ext)=fileparse($f2, ('.fastq'));
$sai = "$outdir/$name.sai"; 
$cbc=$f2; $cbc =~ s/(\.)\w+$/\_cbc.fastq/;

if ( -f $cbc  ) { 
  print "*** Seeing file $cbc, not running add_bc_to_R2.pl to re-create it\n";
} else { 
  $str = "add_bc_to_R2.pl -fastq=$f1,$f2 -rb_len=$rb_len -cbc_len=$cbc_len > $cbc ";
  print $str."\n";
  execute(cmd=>$str, merge=>0) if ($test == 0);
  check_filesize(file=>$cbc, minsize=>1000);
}

if ( $npr != 2 ){                       # npr is 0 or 1: do mapping
  die "obsolete, see before commit 61a2fce50246ce47 (2016-10-11 15:10:00)" if ($F[$i] =~ /txt/);
  $B = $BR;
  $str = "bwa aln -B $B -q $q -n $aln_n -k $aln_k -l $l -t $t $r $cbc > $sai";
  print $str."\n";
  execute(cmd=>$str, merge=>0) if ($test == 0);
  check_filesize(file=>$sai, minsize=>1000);
  
  if ( $nsam == 0 ){
    $str = "bwa samse -n $n $r $sai $cbc > $out.sam";
  }
  print $str."\n";
  execute(cmd=>$str) if ($test == 0);
  check_filesize(file=>"$out.sam",minsize=>1000);
}                                       # npr!=2

if ( $npr == 0 || $npr == 2){
  $s = 1;
  $s = 0 if $pflag;
  ## if ( $STRT ) # unknown, see before commit 61a2fce50246ce47 (2016-10-11 15:10:00)
  $str = "process_sam_cel384v2.pl -sam=$out.sam -barfile=$bar -rb_len=$rb_len";
  print $str."\n";
  execute(cmd=>$str, merge=>1) if ($test == 0);
}                                       # if ( $npr == 0 || $npr == 2)
