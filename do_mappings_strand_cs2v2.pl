#!/usr/bin/perl -w -s
use tools;
use Carp;

$LIST_SEPARATOR=" ";                    # for interpolating arrays inside strings (default anyway)

if (!($r && $f1 && $out && $t)){
    confess "usage:  -r=REFERENCE     \
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
                 -gff=DATA.gff (optional, produces wig ?? files. Passed to process_sam_cel_seq.pl  )    \
                 -s_flag=1 (optional: strand specific mapping. Passed to process_sam_cel_seq.pl)    \
                 -u=1 (optional: only reads mapping to one strand. Passed to process_sam_cel_seq.pl )    \
                 -uniq=1 (optional, keep only uniquely mapped reads. Passed to process_sam_cel_seq.pl )    \
                 -npr=0,1,2: 0: map and process; 1: only map ; 2: only process
                 -nsam= 0 or 1 (1: do *not* produce new sam file (calls bwa samse/sampe)    \
                 -cel=0 or 1 (1: call process_sam_cel_seq.pl; 0: call process_sam_strand.pl (obsolete) \
                 -bar=cel-seq_barcodes.csv    \
                 -fstr= 1 or 0 ( if 1 only mappings to the sense strand are allowed; passed to process_sam_cel_seq.pl )    \
                 -anno= anno.csv ( optional (CEL-seq), when mapping to the genome, run get_anno.pl on sam file first)    \
                 -rb= 0 or 1 (optional (CEL-seq), use random barcodes. Passed to process_sam_cel_seq.pl)    \
                 -s_flag= 0 or 1 (separate output for sense and antisense hits, passed to process_sam_cel_seq.pl)  \
                 -rb_len= length of random barcode (default = 4), passed to process_sam_cel_seq.pl    \
                 -dprm= 0 or 1 (for CEL-seq: 1: remove pcr duplicates, passed to process_sam_cel_seq.pl) \
                 -test=0 or 1 (latter runs in test mode, doesn't call external programs) \
                 -cel384=0 or 1 (if 1, run the 384 well version)   
";
}

# hard coded:
$n = 100; # maximal number of unpaired reads to output in XA tag
$N = 100; # maximal number of paired reads to output in XA tag
$q = 0 if !$q; # base quality cutoff for trimming
$l = 200 if !$l; # seed length
$i = 0 if !$i; # build reference database
$BL = 0 if !$BL;
$BR = 0 if !$BR;
$npr  = 0 if !$npr;
$nsam = 0 if !$nsam;
$cel = 0 if !$cel;
$fstr = 0 if !$fstr;
$ind = "is" if !$ind;
$anno = 0 if !$anno;
$rb = 0 if !$rb;
$rb_len = 4 if !$rb_len;
$dprm = 0 if !$dprm;
$cel384 = 0 if !$cel384;
$rb_len = 6 if $cel384;                 ## @@ hardcoded? PL
$BR = 8+$rb_len if $cel384; # @@ hardcoded? PL
$aln_n = 0.04 if !$aln_n; # edit distance
$aln_k = 2 if !$aln_k; # edit distance in seed

$test = 0 if !$test;
if ($outdir){
    makedir($outdir);
    $out = $outdir."/".$out;
}
$pflag = 0;
$pflag = 1 if ($f1 && $f2);             # paired

confess "reference genome $r not found" unless -f $r;

if ($i){
    $str = "bwa index -a $ind $r";
    print $str."\n";
    execute(cmd=>$str, merge=>1) if ($test == 0);
}

@F = ($f1);
@F = (@F, $f2) if $pflag;

@fastq = @F;
@sai = @F;
@cbc = @F;

for $i (0..$#F){
  $fastq[$i] =~ s/\.\w+$/\.fastq/; # PL: only meaningful if there are .txt files
  $sai[$i] =~ s/\.\w+$/\.sai/;
}

$cbc[1] =~ s/(\.)\w+$/\_cbc.fastq/;

if ($cel384){
  $str = "add_bc_to_R2.pl -fastq=$f1,$f2 -rb_len=$rb_len > $cbc[1] ";
  print $str."\n";
  execute(cmd=>$str, merge=>1) if ($test == 0);
  check_filesize(file=>$cbc[1], minsize=>1000);
}

if ( $npr != 2 ){                       # npr is 0 or 1: do mapping
  for $i (0..$#fastq){
    if ($F[$i] =~ /txt/){
      $str = "qseq2fastq.pl -clean=1 -in=$F[$i] > $fastq[$i]";
      print $str."\n";
      execute(cmd=>$str) if $test==0;
    }
    if ( $i == 0 ) {
      $B = $BL;
    } else {
      $B = $BR;
    }
    if ( $cel384 == 0){
      $str = "bwa aln -B $B -q $q -n $aln_n -k $aln_k -l $l -t $t $r $fastq[$i] > $sai[$i]";
      print $str."\n";
      execute(cmd=>$str, merge=>0) if ($test == 0);
    }
    if ( $cel384 == 1 && $i == 1){
      $str = "bwa aln -B $B -q $q -n $aln_n -k $aln_k -l $l -t $t $r $cbc[$i] > $sai[$i]";
      print $str."\n";
      execute(cmd=>$str, merge=>0) if ($test == 0);
    }
    check_filesize(file=>$sai[$i], minsize=>1000);
  }
  
  if ( $nsam == 0 ){
    if ($pflag && $cel384 == 0){
      $str = "bwa sampe -n $n -N $N $r @sai @fastq > $out.sam";
    }elsif ($pflag && $cel384 == 1){
      $str = "bwa samse -n $n $r $sai[1] $cbc[1] > $out.sam";
    }else{
      $str = "bwa samse -n $n $r @sai @fastq > $out.sam";
    }
  }
  print $str."\n";
  execute(cmd=>$str) if ($test == 0);
  check_filesize(file=>"$out.sam",minsize=>1000);
}                                       # npr!=2

if ( $npr == 0 || $npr == 2){
    $s = 1;
    $s = 0 if $pflag;
    $str = "process_sam_strand.pl -in=$out.sam -s=$s";
    if ( $cel ){
      $str = "process_sam_cel_seq.pl -in=$out.sam -bc=$bar -anno=$anno -rb=$rb -rb_len=$rb_len";
      if ( $fstr ){
	$str .= " -fstr=$fstr";
      }
      if ( $dprm ){
	$str .= " -dprm=$dprm";
      }
    }
    if ( $STRT ){
      $str = "process_sam_STRT.pl -sam=$out.sam -barfile=$bar -BCset=$BCset -rb_len=$rb_len -se=$s";
    }
    if ( $cel384 ){
      $str = "process_sam_cel384v2.pl -sam=$out.sam -barfile=$bar -rb_len=$rb_len";
    }
    if ($gff){
	$str .= "-gff=$gff";
    }
    if ($s_flag){
      $str .= " -s_flag=$s_flag";
    }
    if ($u){
	$str .= " -u=$u";
    }
    if ($uniq){
	$str .= " -uniq=$uniq";
    }
    print $str."\n";
    execute(cmd=>$str, merge=>1) if ($test == 0);
}                                       # if ( $npr == 0 || $npr == 2)
