#!/usr/bin/perl -s -w

use tools;

if (scalar @ARGV == 1) {
    die "usage:    \
      -in=INPUTFILE.sam    \
      -bc=cel-seq_barcodes.csv     \
      -s_flag= 1 or 0 ( if 1 then produce separate files for sense and antisense strand )    \
      -u= 1 or 0 ( if 1 then only map reads that map to only one strand, optional)    \
      -uniq=1 (optional, only unique reads)    \
      -fstr= 1 or 0 ( if 1 only mappings to the sense strand are allowed )    \
      -anno= anno.csv ( optional (CEL-seq), when mapping to the genome, run get_anno.pl on sam file first)    \
      -rb= 0 or 1 (optional (CEL-seq), use random UMIs)    \
      -rb_len= length of UMI (default = 4)  \
      -dprm = 0 or 1 (for CEL-seq: 1: remove duplicates) \
" if $ARGV[0] eq "help"; 
}

$u = 0 if !$u;
$uniq = 0 if !$uniq;
$s_flag = 0 if !$s_flag;
$fstr = 0 if !$fstr;                    # if 1: ignore mappings to reverse strand
$anno = 0 if !$anno;                    # not used 
$rb = 0 if !$rb;                        # use umis
$rb_len = 4 if !$rb_len;                # umi length
$dprm = 0 if !$dprm;                    # duplicate removal

%t = (); # hash for reads per reference sequence
%tc = (); # hash for reads per barcode and reference sequence
%bcount = (); # count of each barcode
$n = 0; # number of reference sequences
$sh_single = 0; # number of reads discarded because too short
$min_l = 15; # minimal length required for a read to be included

$p_r   = 0; # number of reads
$u_m   = 0; # number of reads mapped
$n_bar = 0; # number of reads not mapped in a proper pair


@F = split(/\//,$in);
$in_name = $F[$#F];
$dir = join("/",@F[0..$#F - 1]);
$in_name =~ s/(\.)\w+$//;

$sout = $in; # output file for summary stats
if ($u){
    $sout =~ s/(\.)\w+$/\.u\.sout/;
}else{
    $sout =~ s/(\.)\w+$/\.sout/;
}

if ( $s_flag ){ # never used
##     $cout_s = $in; # output file for read counts per reference sequence in input file
##     $cout_a = $in; # output file for read counts per reference sequence in input file
##     if ($u){
## 	$cout_s =~ s/(\.)\w+$/\.u\.sense\.cout\.csv/;
## 	$cout_a =~ s/(\.)\w+$/\.u\.antisense\.cout\.csv/;
##     }else{
##     	$cout_s =~ s/(\.)\w+$/\.sense\.cout\.csv/;
## 	$cout_a =~ s/(\.)\w+$/\.antisense\.cout\.csv/;
##     }
}else{
    $cout   = $in;
    if ($u){
	$cout   =~ s/(\.)\w+$/\.u\.cout\.csv/;
    }else{
    	$cout   =~ s/(\.)\w+$/\.cout\.csv/;
    }
}

$bout = $in; # output file for summary stats
$bout =~ s/(\.)\w+$/\.bout/;

my $bl;                                 # barcode length

open(BC,'<',$bc) || die "$bc: $!";
while(<BC>){
  chomp;
  ($name, $barcode)=split(/\t/);
  $bar{$barcode} = $name;                  # $bar{'GACACTCA'} => 23
  $bl = length($barcode); 
}
close(BC);

if ( $anno ne "0" ){ # we don't use this
  die "-anno=1 : we don't have the get_anno.pl script (PL)";
}
##   open(ANNO,"<",$anno);
##   while(<ANNO>){
##     chomp;
##     @F = split(/\t/);
##     $r = $F[3]."\t".$F[4];
##     #$an{$F[0]}{$F[1]}{($F[2] - 1)} = $r;
##     $an{$F[0]}{$F[1]}{$F[2]} = $r;
##     if ( !$rb ){
##       $t{"+"}{$r} = 0;
##       $t{"-"}{$r} = 0;
##       foreach $k (keys %bar){
## 	$tc{$bar{$k}}{"+"}{$r} = 0;
## 	$tc{$bar{$k}}{"-"}{$r} = 0;
##       }
##     }
##     $n ++;
##   }
##   close(ANNO);
## }

if ( $rb ){
  $t{"+"} = (); # %t: reads per reference sequence (?)
  $t{"-"} = (); # %tc: reads per reference sequence + cellcode
  foreach $barcode (keys %bar){
    $tc{$bar{$barcode}}{"+"} = ();            # keys of %tc are 1 .. 96
    $tc{$bar{$barcode}}{"-"} = ();
  }
}

open(IN,'<',$in) || die "$in: $!";
$l_flag = 0;
LINE:

while(<IN>){
  if (/^\@/){                           # reference sequences
##     if (/^\@SQ/ && $anno eq "0" && !$rb){ ## anno not used anymore PL
##       chomp;
##       @r = split(/\t/);
##       $r[1] =~ s/SN://g;
##       $t{"+"}{$r[1]} = 0;
##       $t{"-"}{$r[1]} = 0;
##       foreach $k (keys %bar){
## 	$tc{$bar{$k}}{"+"}{$r[1]} = 0;
## 	$tc{$bar{$k}}{"-"}{$r[1]} = 0;
##       }
##       $n ++;
##     }
  }
  if ($l_flag == 0){
    @l = ($_);
    $l_flag = 1;
    next LINE;
  }else{
    push(@l,$_);                      # read in pairs => must be name sorted, not checked!
    $l_flag = 0;
    ## @d_flag = ();                       # not used
  }
 READ:
  for $i (0..$#l){
    ($QNAME[$i],$FLAG[$i],$RNAME[$i],$POS[$i],$MAPQ[$i],$CIGAR[$i],
     $MRNM[$i],$MPOS[$i],$ISIZE[$i],$SEQ[$i],$QUAL[$i])=split(/\t/,$l[$i]);
    
    @LINE = split(/\t/,$l[$i]);     # ugly, we already splitted $l[$i]
    $NM[$i] = "NA";                 # edit distance
    $XA[$i] = "NA";                 # alternative alignments
    ## $d_flag[$i] = 1;                #  not used
    foreach $el (@LINE[ 10 .. $#@LINE] ){
      ($dum,$dum,$NM[$i]) = split(/\:/,$el) if ($el =~ /^NM\:/);
      ($dum,$dum,$XA[$i]) = split(/\:/,$el) if ($el =~ /^XA\:/);
    }
    ## $d_flag[$i] = 0 if $XA[$i] eq "NA"; # not used
  }                                 # READ
 READ12:
  for $i (0..$#l){
    @flag   = split(//,reverse(dec2bin($FLAG[$i])));

    if ( $i == 0 ){                 # read1
      $dpflag = 0;
      $BAR = substr($SEQ[$i],0,$bl) if  $i == 0  && !$flag[4]; # $bl is cell barcode length
      $BAR = substr(revcompl($SEQ[$i]),0,$bl) if  $i == 0  &&  $flag[4]; # flag[4]==1: read is on reversestrand
      if ( $rb ){
        $RBAR = substr($SEQ[$i],$bl,$rb_len) if  $i == 0  && !$flag[4]; # UMI
        $RBAR = substr(revcompl($SEQ[$i]),$bl,$rb_len) if  $i == 0  &&  $flag[4];
      }
      $dpflag = 1 if exists($dpseen{$BAR}{$RNAME[1]}{$POS[1]});
      $dpseen{$BAR}{$RNAME[1]}{$POS[1]} = 1;
    } else {                        # read2

      next READ12 if $dprm && $dpflag == 1;

      if ( $uniq ){
        @ALT  = split(/;/,$XA[$i]);
        $u_flag = 0;

      ALT:
        foreach $xa (@ALT){
          last if $xa eq "NA" || $xa =~ /\s+/;
          ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
          $u_flag = 1 if ($nm <= $NM[$i]);
        }
        next if $u_flag;            # @@@ should this be next READ12 or next ALT?
      }

      if ($flag[4]){$str = "-";}else{$str = "+";}
      
      $start  = $POS[$i] - 1;
      $length = length($SEQ[$i]);
      $p_r ++;                      # read counter
      
      if (length($SEQ[$i]) < $min_l){
        $sh_single ++;
        ## why is there no 'next' here ?!?! @@@ PL
      }

      if ($flag[2] == 0){           # read was mapped
        $bcount{$BAR} = 0 if !exists($bcount{$BAR});
        $bcount{$BAR} ++ if   exists($bcount{$BAR});
        if (exists($bar{$BAR})){
          $u_m ++;                  # barcode known
          update_t($RNAME[$i], $NM[$i], $XA[$i], \%t, \%{$tc{$bar{$BAR}}}, $str, $u, $start,
                   # following arguments are ignored:
                   $length, \%w, $s_flag);
        }else{
          $n_bar ++;                # barcode unknown
        }
      } else { 
        # read was not mapped
      }
    }                               # end of read2
  }                                 # READ12
}                                       # while <IN>
close(IN);

open(BOUT,'>',$bout) || die "$bout: $!";
foreach (keys %bcount){
  print BOUT $_."\t".$bcount{$_};
  print BOUT "\t".$bar{$_}."\n" if exists($bar{$_});
  print BOUT "\t\n" if !exists($bar{$_});
}
close(BOUT);

open(SOUT,'>',$sout) || die "$sout: $!";

print SOUT "number of reference sequences:\t$n\n";
print SOUT "number of reads:\t".$p_r ."\n";
print SOUT "applied constraint on read length:\t$min_l (NOT TRUE, PL)\n";
print SOUT "number of reads discarded due to length constraint:\t$sh_single\n";
print SOUT "number of reads mapped with valid barcode:\t$u_m\n";
print SOUT "number of reads mapped without valid barcode:\t$n_bar\n";
print SOUT "fraction of reads mapped with valid barcode:\t".($u_m/$p_r)."\n";
print SOUT "fraction of reads mapped without valid barcode:\t".($n_bar/$p_r)."\n";

$a_sum = 0;
$s_sum = 0;

if ($s_flag){ # never used
##    open(COUTA,'>',$cout_a) if !$fstr;
##    open(COUTS,'>',$cout_s);
}else{
  open(COUT,'>',$cout) || die "$cout: $!";
}

$name = "GENEID";
if ($anno ne "0" ){ $name = "CLASS\t".$name}
if ( $rb ){ $name = $name."\tRBAR"}
if ($s_flag){                           # never used
##  print COUTA join("\t",($name, sort {$a <=> $b} keys %tc))."\n" if !$fstr;
##  print COUTS join("\t",($name, sort {$a <=> $b} keys %tc))."\n";
}else{
  print COUT  join("\t",($name, sort {$a <=> $b} keys %tc))."\n";
}
foreach $k (keys %{$t{"+"}}){
  $seen{$k} = 1;
}
foreach $k (keys %{$t{"-"}}){
  $seen{$k} = 1;
}

foreach $k (keys %seen){
  if ($s_flag){
##     print COUTS $k;
##     print COUTA $k if !$fstr;
##     foreach $id ( sort {$a <=> $b} keys %tc ){
##       $pl = 0;
##       $mi = 0;
##       $pl = $tc{$id}{"+"}{$k} if exists($tc{$id}{"+"}{$k});
##       $mi = $tc{$id}{"-"}{$k} if exists($tc{$id}{"-"}{$k});
##       print COUTS "\t".$pl;
##       print COUTA "\t".$mi if !$fstr;
##     }
##     print COUTS "\n";
##     print COUTA "\n" if !$fstr;
  }else{
    print COUT $k;
    foreach $id ( sort {$a <=> $b} keys %tc ){
      $pl = 0;
      $mi = 0;
      $pl = $tc{$id}{"+"}{$k} if exists($tc{$id}{"+"}{$k});
      $mi = $tc{$id}{"-"}{$k} if exists($tc{$id}{"-"}{$k});
      if ($fstr){
	print COUT "\t".$pl;
      }else{
	print COUT "\t".($pl + $mi);
      }
    }
    print COUT "\n";
  }
  $t{"+"}{$k} = 0 if !exists($t{"+"}{$k});
  $s_sum += $t{"+"}{$k};
  $t{"-"}{$k} = 0 if !exists($t{"-"}{$k});
  $a_sum += $t{"-"}{$k} if !$fstr;
}
print SOUT "number of reads mapped to sense strand:\t".$s_sum."\n";
print SOUT "number of reads mapped to antisense strand:\t".$a_sum."\n";
if ( $s_sum + $a_sum > 0 ){
    $fr = $s_sum/($s_sum + $a_sum);
}else{
    $fr = 0;
}
print SOUT "fraction of reads mapped to sense strand:\t".$fr."\n";
if ( $s_sum + $a_sum > 0 ){
    $fr = $a_sum/($s_sum + $a_sum);
}else{
    $fr = 0;
}
print SOUT "fraction of reads mapped to antisense strand:\t".$fr."\n";

close(SOUT);
if ($s_flag){                           # never used
##    close(COUTA) if !$fstr;
##    close(COUTS);
}else{
    close(COUT);
}

sub update_t {
  my $RNAME = shift;
  my $NM    = shift;
  my $XA    = shift;
  my $t     = shift;
  my $tc    = shift;
  my $st    = shift;
  my $u     = shift;
  my $start = shift;
  
  my @ALT   = split(/;/,$XA);
  my %HITS   = ();
  my %COORDS = ();
  my $h_nb  = 0;
  if ( ( $fstr && $st eq "+" ) || !$fstr ){
    if ($anno ne "0"){
      if ( exists($an{$RNAME}{$st}{$start})){
	$h_nb  = 1;
	fill_histo(\%{$HITS{$st}},$an{$RNAME}{$st}{$start});
      }
    }else{
      $h_nb  = 1;
      fill_histo(\%{$HITS{$st}},$RNAME);
    }
  }
  foreach $xa (@ALT){
    last if $xa eq "NA" || $xa =~ /\s+/;
    my ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
    if ($pos < 0){$xa_st = "-"; $pos = -$pos;}else{$xa_st = "+";}
    $pos = $pos - 1;
    next if $fstr && $xa_st eq "-";
    if ($nm <= $NM){
      if ($anno ne "0"){
## 	if ( exists($an{$alt_rname}{$xa_st}{$pos})){
## 	  $h_nb ++;
## 	  fill_histo(\%{$HITS{$st}},$an{$alt_rname}{$xa_st}{$pos});
## 	}
      } else{
	fill_histo(\%{$HITS{$xa_st}},$alt_rname);
	$h_nb ++;
      } 
    }
  }                                     # foreach $xa

  if (!$u || !(exists($HITS{"+"}) && exists($HITS{"-"}))){
    foreach $s (keys %HITS){
      foreach $rname (keys %{$HITS{$s}}){
	if ($rb){
 	  $$t{$s}{$rname."\t".$RBAR}  = 0 if !exists($$t{$s}{$rname."\t".$RBAR});
 	  $$tc{$s}{$rname."\t".$RBAR} = 0 if !exists($$tc{$s}{$rname."\t".$RBAR});
 	  $$t{$s}{$rname."\t".$RBAR}  += $HITS{$s}{$rname}/$h_nb;
 	  $$tc{$s}{$rname."\t".$RBAR} += $HITS{$s}{$rname}/$h_nb;
	}else{
	  $$t{$s}{$rname}  += $HITS{$s}{$rname}/$h_nb;
	  $$tc{$s}{$rname} += $HITS{$s}{$rname}/$h_nb;				    
	}
      }                                 # foreach $rname
    }                                   # foreach $s
  }
}                                       # update_t


