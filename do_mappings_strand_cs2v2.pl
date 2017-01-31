#!/usr/bin/perl -w
## driver script for the single-cell RNAseq scripts, originally written by Dominic Grün, 
## heavily modified by Lennert Kestler and Philip Lijnzaad
use strict;
use Carp;
use File::Basename;
use Getopt::Long;

use tools;

my $version=getversion($0);
warn "Running $0, version $version\nwith arguments:\n @ARGV\n";

$,=" ";         # for interpolating arrays inside strings (default anyway)

our($r, $f1, $f2, $out, $bwaparams, $outdir, $t, $ind, $q, $aln_n, $aln_k, $l,
    $i, $npr, $bar, $umi_len, $cbc_len, $trim, $allow_mm, $test, $protocol, $help, $cbc);

my $script=basename($0);

my $usage = "Usage:  $script ARGUMENTS  [ 2> logfile ]

Mandatory arguments:

  --r file.fa            Fasta file containing the reference transcriptome
  --f1 FILE_R1.fastq.gz  Input file with R1 reads (if read1 and read2 not yet preprocessed into *_cbc.fastq.gz)   
  --f2 FILE_R2.fastq.gz  Input file with R2 reads (if read1 and read2 not yet preprocessed into *_cbc.fastq.gz)   
  --cbc FILE_cbc.fastq.gz Input file with reads + barcodes  (if --f2 and --f2 are already preprocessed into _cbc.fastq.gz)
  --out OUTPUT_PREFIX    Prefix for the resulting bam file

Optional arguments:

  --outdir OUTPUT_DIRECTORY (optional) Directory into which to put results
  --t THREADS  Number of threads bwa can use  (default 1)
  --ind 'is' or 'bwtsw' (default: 'is' (<2GB reference)    
  --q PHRED_QUALITY_FOR_TRIMMING (0..62, default 0)   
  --aln_n edit distance  (bwa aln option -n, default 0.04)  
  --aln_k edit distance in seed  (bwa aln option -k, default 2)   
  --l SEED_LENGTH  (bwa -l option)  
  --i 1 or 0 (1 if indexing is required, runs bwa index on reference transcriptome. Default 0)
  --bwaparams '...' (optional, overrides the -q -n -k -l options) 
  --npr 0,1,2: 0: map and process; 1: only map ; 2: only process 
  --bar file  Name of file with CEL-seq2 barcodes (format: id \\t sequence)
  --umi_len N length of UMI (default = 6)  
  --cbc_len N length of cell barcode (default: 8) 
  --trim string Passed to preprocess_fastq.pl for trimming. Example: 'A12,T=14'
  --allow_mm N Passed to process_sam_cel384v2, allows N mismatches in cell bar codes
  --test N   0 or 1 (latter runs in test mode, doesn't call external programs) 
  --protocol 1 or 2  Default 2 (celseq2): umi=6, cbc=8; 1 (celseq1): umi=4, cbc=8 
";

die $usage unless GetOptions(
  "r|ref=s"	=> \$r,
  "1|f1|read1=s"=> \$f1,
  "2|f2|read2=s"=> \$f2,
  "cbc=s"=>     => \$cbc,
  "out|prefix=s"=> \$out,
  "bwaparams=s"	=> \$bwaparams,
  "outdir=s"	=> \$outdir,
  "t|threads=i"	=> \$t,
  "ind|indextype=s"=> \$ind,
  "q|phred_quality=i" => \$q,
  "aln_n=f"	=> \$aln_n,
  "aln_k=i"	=> \$aln_k,
  "l|seedlen=i"	=> \$l,
  "i|doindexing=i"  => \$i,
  "npr|processing=i"=> \$npr,
  "bar|barcodes=s"  => \$bar,
  "umi_len=i"	=> \$umi_len,
  "cbc_len=i"	=> \$cbc_len,
  "trim=s"	=> \$trim,
  "allow_mm=i"	=> \$allow_mm,
  "test=i"	=> \$test,
  "protocol=i"	=> \$protocol,
  "h|help"	=> \$help       
    );

die $usage if ($help || !$r || !$out);
die "Specify either --cbc, or --f1 and --f2; \nin the latter case\n" . $usage if (!!($f1 && $f2)  == !!$cbc);
die "Don't mix any of -[ql] or aln_{k,n} with --bwaparams! \n" . $usage if ( $bwaparams && ($q || $aln_n || $aln_k || $l)  );

# hard coded:
my $n = 100; # maximal number of unpaired reads to output in XA tag
my $N = 100; # maximal number of paired reads to output in XA tag

$q = 0 if !$q; # base quality cutoff for trimming
$l = 200 if !$l; # seed length

$t=1 if !$t;                            # threads
$i = 0 if !$i;                          # create index
$npr  = 0 if !$npr;
$ind = "is" if !$ind;
$protocol = 2 if !$protocol;

$umi_len = ({1=>4,2=>6})->{$protocol} if !$umi_len;
$cbc_len = 8 if !$cbc_len;

$aln_n = 0.04 if !$aln_n; # edit distance
$aln_k = 2 if !$aln_k; # edit distance in seed

$bwaparams=" -q $q -n $aln_n -k $aln_k -l $l " unless $bwaparams;

$trim="" unless $trim;

$trim = "--trim $trim" if $trim;

if ($allow_mm) { 
  $allow_mm="--allow_mm $allow_mm";
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
  my $str = "bwa index -a $ind $r 2>&1 ";
  print $str."\n";
  execute($str) if ($test == 0);
}

if ( $npr != 2 ) {                      # npr is 0 or 1: do mapping
  my ($name, $path, $ext);
  if ( $f1 && $f2  ) {                  # create cbc file now
    ($name, $path, $ext)=fileparse($f2, '.fastq.gz');
    die "$f2: extension must be .fastq.gz " unless $ext eq '.fastq.gz';
    $cbc = "$path/${name}_cbc.fastq.gz";
    die "$f1: $!" unless -f $f1;
    die "$f2: $!" unless -f $f2;
    print "Creating file $cbc from $f1 and $f2 ...\n";
    my($log1, $log2)=(openlog("preprocessLOG-$version"), openlog("preprocessZIP-$version"));
    my $str = "preprocess_fastq.pl --fastq $f1,$f2 --umi_len $umi_len --cbc_len $cbc_len $trim  2>$log1 | gzip > $cbc 2> $log2 ";
    print $str."\n";
    execute($str) if ($test == 0);
    check_filesize(file=>$cbc, minsize=>1000);
    dumplog($log1);
    dumplog($log2);
  } else {
    ($name, $path, $ext)=fileparse($cbc, '_cbc.fastq.gz');
    die "$cbc: extension must be _cbc.fastq.gz" unless $ext eq '_cbc.fastq.gz';
    die "$cbc: $!" unless -f $cbc;
  } 

  my $sai = "$outdir/$name.sai"; 
  my($log)=openlog("bwa_alnLOG-$version");
  my $mapping = "bwa aln -B 0 -t $t $bwaparams $r $cbc ";
  my $cmdline= "$mapping > $sai 2>$log";
  print "$cmdline\n";
  dumplog($log);
  my $status=execute($cmdline) if ($test == 0);
  confess "non-zero exit status: $status" unless $status==0;
  check_filesize(file=>$sai, minsize=>1000);
  
  my $sam2bam = "samtools view -h -b - ";
  my($samse_log, $PGline_log, $sam2bam_log)=(openlog("bwasamse-LOG-$version"),
                                             openlog("addbwaPGline-LOG-$version"),
                                             openlog("sam2bam-LOG-$version"));
  $cmdline = "bwa samse -n $n $r $sai $cbc 2>$samse_log | add_bwa_PG_line.pl $mapping 2> $PGline_log | $sam2bam > $out.bam 2> $sam2bam_log";

  print $cmdline."\n";
  $status=execute($cmdline) if ($test == 0);
  dumplog($samse_log);
  dumplog($PGline_log);
  dumplog($sam2bam_log);
  check_filesize(file=>"$out.bam",minsize=>1000);
  confess "non-zero exit status: $status" unless $status==0;
  unlink($sai);
}                                       # npr!=2

if ( $npr == 0 || $npr == 2){
  my($log)=openlog("process_samBOTH-$version");
  my $cmdline = "process_sam_cel384v2.pl --barcodefile $bar --umi_len $umi_len --cbc_len $cbc_len $allow_mm $out.bam > $log 2>&1 ";
  print $cmdline."\n";
  execute($cmdline) if ($test == 0);
  dumplog($log);
}                                       # if ( $npr == 0 || $npr == 2)
