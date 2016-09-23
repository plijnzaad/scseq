#!/usr/bin/perl -s -w

# use lib '/hpc/hub_oudenaarden/bin/';
# use lib '/Users/d.grun/data/bin/';
use tools;

if (scalar @ARGV == 1)
{
    die "usage: -in=INPUTFILE.sam -bc=cel-seq_barcodes.csv  -s_flag= 1 or 0 ( if 1 then produce separate files for sense and antisense strand ) -u= 1 or 0 ( if 1 then only map reads that map to only one strand, optional) -uniq=1 (optional, only unique reads) -fstr= 1 or 0 ( if 1 only mappings to the sense strand are allowed ) -anno= anno.csv ( optional (CEL-seq), when mapping to the genome, run get_anno.pl on sam file first) -rb= 0 or 1 (optional (CEL-seq), use random barcodes) -rb_len= length of random barcode (default = 4)\n" if $ARGV[0] eq "help";
}

$u = 0 if !$u;
$uniq = 0 if !$uniq;
$s_flag = 0 if !$s_flag;
$fstr = 0 if !$fstr;
$anno = 0 if !$anno;
$rb = 0 if !$rb;
$rb_len = 4 if !$rb_len;
$dprm = 0 if !$dprm;

%t = (); # hash for reads per reference sequence
%tc = (); # hash for reads per barcode and reference sequence
%bcount = (); # count of each barcode
$n = 0; # number of reference sequences
$sh_single = 0; # number of reads discarded due to length constraint
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

if ( $s_flag ){
    $cout_s = $in; # output file for read counts per reference sequence in input file
    $cout_a = $in; # output file for read counts per reference sequence in input file
    if ($u){
	$cout_s =~ s/(\.)\w+$/\.u\.sense\.cout\.csv/;
	$cout_a =~ s/(\.)\w+$/\.u\.antisense\.cout\.csv/;
    }else{
    	$cout_s =~ s/(\.)\w+$/\.sense\.cout\.csv/;
	$cout_a =~ s/(\.)\w+$/\.antisense\.cout\.csv/;
    }
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


open(BC,'<',$bc);
while(<BC>){
  chomp;
  @F=split(/\t/);
  $bar{$F[1]} = $F[0];
  $bl = length($F[1]);
}
close(BC);

if ( $anno ne "0" ){
  open(ANNO,"<",$anno);
  while(<ANNO>){
    chomp;
    @F = split(/\t/);
    $r = $F[3]."\t".$F[4];
    #$an{$F[0]}{$F[1]}{($F[2] - 1)} = $r;
    $an{$F[0]}{$F[1]}{$F[2]} = $r;
    if ( !$rb ){
      $t{"+"}{$r} = 0;
      $t{"-"}{$r} = 0;
      foreach $k (keys %bar){
	$tc{$bar{$k}}{"+"}{$r} = 0;
	$tc{$bar{$k}}{"-"}{$r} = 0;
      }
    }
    $n ++;
  }
  close(ANNO);
}

if ( $rb ){
  $t{"+"} = ();
  $t{"-"} = ();
  foreach $k (keys %bar){
    $tc{$bar{$k}}{"+"} = ();
    $tc{$bar{$k}}{"-"} = ();
  }
}

open(IN,'<',$in);
$l_flag = 0;
while(<IN>){
  if (/^\@/){
    if (/^\@SQ/ && $anno eq "0" && !$rb){
      chomp;
      @r = split(/\t/);
      $r[1] =~ s/SN://g;
      $t{"+"}{$r[1]} = 0;
      $t{"-"}{$r[1]} = 0;
      foreach $k (keys %bar){
	$tc{$bar{$k}}{"+"}{$r[1]} = 0;
	$tc{$bar{$k}}{"-"}{$r[1]} = 0;
      }
      $n ++;
    }
  }else{
    if ($l_flag == 0){
      @l = ($_);
      $l_flag = 1;
      next;
    }else{
      push(@l,$_);
      $l_flag = 0;
      @d_flag = ();
      for $i (0..$#l){
	($QNAME[$i],$FLAG[$i],$RNAME[$i],$POS[$i],$MAPQ[$i],$CIGAR[$i],
	 $MRNM[$i],$MPOS[$i],$ISIZE[$i],$SEQ[$i],$QUAL[$i])=split(/\t/,$l[$i]);
	
	@LINE = split(/\t/,$l[$i]);
	$NM[$i] = "NA";
	$XA[$i] = "NA";
	$d_flag[$i] = 1;
	foreach $el (@LINE){
	  ($dum,$dum,$NM[$i]) = split(/\:/,$el) if ($el =~ /^NM\:/);
	  ($dum,$dum,$XA[$i]) = split(/\:/,$el) if ($el =~ /^XA\:/);
	}
	$d_flag[$i] = 0 if $XA[$i] eq "NA";
      }
      for $i (0..$#l){
	@flag   = split(//,reverse(dec2bin($FLAG[$i])));
	if ( $i == 0 ){
	  $dpflag = 0;
	  $BAR = substr($SEQ[$i],0,$bl) if  $i == 0  && !$flag[4];
	  $BAR = substr(revcompl($SEQ[$i]),0,$bl) if  $i == 0  &&  $flag[4];
	  if ( $rb ){
	    $RBAR = substr($SEQ[$i],$bl,$rb_len) if  $i == 0  && !$flag[4];
	    $RBAR = substr(revcompl($SEQ[$i]),$bl,$rb_len) if  $i == 0  &&  $flag[4];
	  }
	  $dpflag = 1 if exists($dpseen{$BAR}{$RNAME[1]}{$POS[1]});
	  $dpseen{$BAR}{$RNAME[1]}{$POS[1]} = 1;

	  #print STDERR ">".$RNAME[$i]."\n".revcompl(substr($SEQ[$i],$bl+$rb_len,length($SEQ[$i]) - ($bl + $rb_len)))."\n" if  $i == 0  && !$flag[4];
	  #print STDERR ">".$RNAME[$i]."\n".revcompl(substr(revcompl($SEQ[$i]),$bl+$rb_len,length($SEQ[$i]) - ($bl + $rb_len)))."\n" if  $i == 0  && $flag[4];
	}else{
	  #next if $uniq && $d_flag[$i] == 1;
	  next if $dprm && $dpflag == 1;
	  if ( $uniq ){
	    @ALT  = split(/;/,$XA[$i]);
	    $u_flag = 0;
	    foreach $xa (@ALT){
	      last if $xa eq "NA" || $xa =~ /\s+/;
	      ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
	      $u_flag = 1 if ($nm <= $NM[$i]);
	    }
	    next if $u_flag;
	  }



	  if ($flag[4]){$str = "-";}else{$str = "+";}
	  
	  #print STDERR join("",@flag)."\t".$SEQ[$i]."\t".$BAR."\n";
	  
	  $start  = $POS[$i] - 1;
	  $length = length($SEQ[$i]);
	  #print STDERR  join("",@flag)."\t".$str."\t".$RNAME[$i]."\t".$CIGAR[$i]."\t".$SEQ[$i]."\t".$start."\t".$length."\n";
	  $p_r ++;
	  
	  if (length($SEQ[$i]) < $min_l){
	    $sh_single += 1;
	  }
	  if ($flag[2] == 0){
	    $bcount{$BAR} = 0 if !exists($bcount{$BAR});
	    $bcount{$BAR} ++ if   exists($bcount{$BAR});
	    if (exists($bar{$BAR})){
	      $u_m ++;
	      #print STDERR  join("",@flag)."\t".$str."\t".$RNAME[$i]."\t".$CIGAR[$i]."\t".$SEQ[$i]."\t".$start."\t".$length."\n";
	      update_t($RNAME[$i], $NM[$i], $XA[$i], \%t, \%{$tc{$bar{$BAR}}}, $str, $u, $start, $length, \%w, $s_flag);
	    }else{
	      $n_bar ++;
	    }
	  }
	}
      }
    }
  }
}
close(IN);
open(BOUT,'>',$bout);
foreach (keys %bcount){
  print BOUT $_."\t".$bcount{$_};
  print BOUT "\t".$bar{$_}."\n" if exists($bar{$_});
  print BOUT "\t\n" if !exists($bar{$_});
}
close(BOUT);

open(SOUT,'>',$sout);

print SOUT "number of reference sequences:\t".$n."\n";
print SOUT "number of reads:\t".$p_r ."\n";
print SOUT "applied constraint on read length:\t".$min_l."\n";
print SOUT "number of reads discarded due to length constraint:\t".$sh_single."\n";
print SOUT "number of reads mapped with valid barcode:\t".$u_m."\n";
print SOUT "number of reads mapped without valid barcode:\t".$n_bar."\n";
print SOUT "fraction of reads mapped with valid barcode:\t".($u_m/$p_r)."\n";
print SOUT "fraction of reads mapped without valid barcode:\t".($n_bar/$p_r)."\n";

$a_sum = 0;
$s_sum = 0;

if ($s_flag){
    open(COUTA,'>',$cout_a) if !$fstr;
    open(COUTS,'>',$cout_s);
}else{
    open(COUT,'>',$cout);
}

$name = "GENEID";
if ($anno ne "0" ){ $name = "CLASS\t".$name}
if ( $rb ){ $name = $name."\tRBAR"}
if ($s_flag){
  print COUTA join("\t",($name, sort {$a <=> $b} keys %tc))."\n" if !$fstr;
  print COUTS join("\t",($name, sort {$a <=> $b} keys %tc))."\n";
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
    print COUTS $k;
    print COUTA $k if !$fstr;
    foreach $id ( sort {$a <=> $b} keys %tc ){
      $pl = 0;
      $mi = 0;
      $pl = $tc{$id}{"+"}{$k} if exists($tc{$id}{"+"}{$k});
      $mi = $tc{$id}{"-"}{$k} if exists($tc{$id}{"-"}{$k});
      print COUTS "\t".$pl;
      print COUTA "\t".$mi if !$fstr;
    }
    print COUTS "\n";
    print COUTA "\n" if !$fstr;
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
if ($s_flag){
    close(COUTA) if !$fstr;
    close(COUTS);
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
      #print STDERR join("\t",($RNAME,$st,$start,$an{$RNAME}{$st}{$start}))."\n";
      if ( exists($an{$RNAME}{$st}{$start})){
	$h_nb  = 1;
	#print STDERR join("\t",($RNAME,$st,$start,$an{$RNAME}{$st}{$start}))."\n";
	fill_histo(\%{$HITS{$st}},$an{$RNAME}{$st}{$start});
      }
      #print STDERR join("\t",($RNAME,$st,$start,$l[$i]))."\n" if !exists($an{$RNAME}{$st}{$start});
    }else{
      $h_nb  = 1;
      fill_histo(\%{$HITS{$st}},$RNAME);
      #print STDERR join("\t",("here",$st,$RNAME,$HITS{$st}{$RNAME},$h_nb))."\n";
    }
  }
  foreach $xa (@ALT){
    #print STDERR join("\t",($XA,"here",$st,$RNAME,$HITS{$st}{$RNAME},$h_nb))."\n";
    last if $xa eq "NA" || $xa =~ /\s+/;
    my ($alt_rname, $pos, $CIGAR, $nm) = split(/,/,$xa);
    if ($pos < 0){$xa_st = "-"; $pos = -$pos;}else{$xa_st = "+";}
    $pos = $pos - 1;
    next if $fstr && $xa_st eq "-";
    if ($nm <= $NM){
      if ($anno ne "0"){
	if ( exists($an{$alt_rname}{$xa_st}{$pos})){
	  $h_nb ++;
	  fill_histo(\%{$HITS{$st}},$an{$alt_rname}{$xa_st}{$pos});
	}
	#print STDERR join("\t",($alt_rname,$xa_st,$pos,$l[$i],"XA"))."\n" if !exists($an{$alt_rname}{$xa_st}{$pos});
      }else{
	fill_histo(\%{$HITS{$xa_st}},$alt_rname);
	$h_nb ++;
      } 
    }
  }
  if (!$u || !(exists($HITS{"+"}) && exists($HITS{"-"}))){
    foreach $s (keys %HITS){
      foreach $rname (keys %{$HITS{$s}}){
	#print STDERR join("\t",($s,$rname,$HITS{$s}{$rname},$h_nb))."\n";
	if ($rb){
	  $$t{$s}{$rname."\t".$RBAR}  = 0 if !exists($$t{$s}{$rname."\t".$RBAR});
	  $$tc{$s}{$rname."\t".$RBAR} = 0 if !exists($$tc{$s}{$rname."\t".$RBAR});
	  $$t{$s}{$rname."\t".$RBAR}  += $HITS{$s}{$rname}/$h_nb;
	  $$tc{$s}{$rname."\t".$RBAR} += $HITS{$s}{$rname}/$h_nb;
	  #print STDERR join("\t",($s,$rname,$$t{$s}{$rname."\t".$RBAR}))."\n";
	}else{
	  $$t{$s}{$rname}  += $HITS{$s}{$rname}/$h_nb;
	  $$tc{$s}{$rname} += $HITS{$s}{$rname}/$h_nb;				    
	}
      }
    }
  }
}


