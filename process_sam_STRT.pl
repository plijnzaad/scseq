#!/usr/bin/perl -w -s
use List::Util 'sum';

if (scalar @ARGV == 1){
    die "usage: -fastq=s_X_1.fastq -sam=sam.csv -barfile=barfile.csv -rb_len=rb_len -se=1 (for single end sequencing)" if $ARGV[0] eq "help";
}


@cells = ();
@cbc = ();
$i = 0;
$x = 1;
%bar = ();
%tc = ();
@BCsets = split(/\,/,$BCset);
$rb_len = 6;	
$se = 0 if !$se;

# open barcode file, read all 384 barcodes, combine them with the setbarcodes and put them in a hash, 
# the values are cell-numbers, 1-384 for the first set, 385-768 for the second set etc.;
$log = $sam;
$log   =~ s/(\.)\w+$/\.log/;
open(LOG, ">", $log);
print LOG "Process sam started\n";

for $i (0..$#BCsets){
	open(IN,"<",$barfile);
	while(<IN>){	
		chomp;
		@line = split(/\t/,$_);
		$bar{$BCsets[$i].".".substr($line[1],0,7)} = ($i * 384) + $line[0];
		push(@cells, $BCsets[$i].".".substr($line[1],0,7))
	}
}
close(IN);


# create a barcode file that holds all the cell specific barcodes from the fastq file, this barcode file start with 
# the header of the sam file so when reading through the sam file the lines match up
$bar = $sam;
$bar   =~ s/(\.)\w+$/\.bar/;
open(BAR, ">", $bar);
open(IN,"<",$sam);
while(<IN>){
	chomp;

	if (substr($_,1,2) eq "SQ" ){
		print BAR $_."\n";
		next;
	}
	if (substr($_,1,2) eq "PG" ){
		next;
		print BAR $_."\n";
	}
	else {
		last;
	}
}
close(IN);

# extract index sequences from fastq file and put them in the bar file
$k = 0;

open(IN,"<",$fastq);
while(<IN>){
  	chomp;
  	if ($k == $x*1000000){
		print LOG $x." million barcodes processed\n";
		$x++;
	}
#	if ($k > 1000000){
#		last;
#	}
  	if ($i == 0){
  		print BAR substr($_,-7)."\n";
  		if ($se == 0){
  			print BAR substr($_,-7)."\n";
  		}
  		$k++;

	  }  
  	$i++;
	if ($i == 4){
	  	$i = 0;
	}
	
  
}
close(IN);
close(BAR);
print LOG "all barcodes processed\n";

$i = 0;
$j = 0;
$x = 1;
$tot_map_reads = 0;

#create arrays to store the different lines of the sam file;
@r1 = ();
@r2 = ();

#create arrays to store the different components of the sam file for R1 and R2;
@QNAME = ();
@FLAG = ();
@RNAME = ();
@POS = ();
@MAPQ = ();
@CIGAR = ();
@MRNM = ();
@MPOS = ();
@ISIZE = ();
@SEQ = ();
@QUAL = ();
@NM = ();
@XA = ();
@X0 = ();
$RBC = ();
$setbc = "NA";
$UMI = "NA";
$TSO = "NA";
@bar_map = ();
%bar_map_h = ();

#print $se."\n";
# read through sam file create a hash with all genes and cells and extract mapped reads into the hash

open(IN,"<",$sam);
open(BAR,"<",$bar); 
while(not eof IN and not eof BAR){
	$line1 = <IN>;
	$line2 = <BAR>;
	chomp $line1;
	chomp $line2;

#read through the header of the sam-file and barcode file;

	if (substr($line1,1,2) eq "SQ" ){

		next;
	}
	if (substr($line1,1,2) eq "PG" ){
		next;
	}

#	if ($i > 1000000){
#		last;
#	}
	
# j = 0 for read 1, so here all the information is extracted for read 1;	
	if ($se == 0){
		if ($j == 0){
			@r1 = split(/\t/,$line1);
			($QNAME[$j],$FLAG[$j],$RNAME[$j],$POS[$j],$MAPQ[$j],$CIGAR[$j],$MRNM[$j],$MPOS[$j],
			$ISIZE[$j],$SEQ[$j],$QUAL[$j],$RBC)=split(/\t/,$line1);
			$setbc = substr($RBC,5,2);
			$UMI = substr($RBC,7,$rb_len);
			$TSO = substr($RBC,7+$rb_len,3);
			


			$NM[$j] = "NA";
			$XA[$j] = "NA";
			foreach $el (@r1){
				($dum,$dum,$NM[$j]) = split(/\:/,$el) if ($el =~ /^NM\:/);
				($dum,$dum,$XA[$j]) = split(/\:/,$el) if ($el =~ /^XA\:/);	
			}
		
			$j++;
			next;
		}
	
# j = 1 for read 2, so here all information is extracted for read 2;
	
		if ($j == 1){
			@r2 = split(/\t/,$line1);
			($QNAME[$j],$FLAG[$j],$RNAME[$j],$POS[$j],$MAPQ[$j],$CIGAR[$j],$MRNM[$j],$MPOS[$j],
			$ISIZE[$j],$SEQ[$j],$QUAL[$j])=split(/\t/,$line1);
			$NM[$j] = "NA";
			$XA[$j] = "NA";
			$X0[$j] = 0;
			foreach $el (@r2){
				($dum,$dum,$NM[$j]) = split(/\:/,$el) if ($el =~ /^NM\:/);
				($dum,$dum,$XA[$j]) = split(/\:/,$el) if ($el =~ /^XA\:/);			
				($dum,$dum,$X0[$j]) = split(/\:/,$el) if ($el =~ /^X0\:/);			
			}
			$currbar = $setbc.".".$line2;
			
# check if the read fullfills the criteria in the if statements;
			if ($X0[$j] > 0){
				$tot_map_reads++;
			}		
			if (exists $bar{$currbar}){
				if ($X0[$j] == 1 && $RNAME[0] eq $RNAME[1] && $TSO eq "GGG"){
#					if ($FLAG[1] == 153 or $FLAG[1] == 145){
						$tc{$RNAME[$j]}{$currbar}{$UMI}++;
#					}
#					if ($FLAG[1] == 137 or $FLAG[1] == 129){
#						$tc{$RNAME[$j]}{$currbar}{$UMI}++;
#					}
				}
			}
			$j = 0;		
		}
		$i++;
		if ($i == $x*1000000){
			print LOG $x." million reads processed\n";
			$x++;
		}	
	}
	if ($se == 1){
		@r1 = split(/\t/,$line1);
		($QNAME[$j],$FLAG[$j],$RNAME[$j],$POS[$j],$MAPQ[$j],$CIGAR[$j],$MRNM[$j],$MPOS[$j],
		$ISIZE[$j],$SEQ[$j],$QUAL[$j],$RBC)=split(/\t/,$line1);
		$setbc = substr($RBC,5,2);
		$UMI = substr($RBC,7,$rb_len);
		$TSO = substr($RBC,7+$rb_len,3);	

		$NM[$j] = "NA";
		$XA[$j] = "NA";
		$X0[$j] = 0;
		foreach $el (@r1){
			($dum,$dum,$NM[$j]) = split(/\:/,$el) if ($el =~ /^NM\:/);
			($dum,$dum,$XA[$j]) = split(/\:/,$el) if ($el =~ /^XA\:/);	
			($dum,$dum,$X0[$j]) = split(/\:/,$el) if ($el =~ /^X0\:/);			
		}
		$currbar = $setbc.".".$line2;
		if (!($RNAME[$j] eq "*")){
			$tot_map_reads++;
			push(@bar_map,$currbar);
			$bar_map_h{$currbar}++;

		}
		if (exists $bar{$currbar}){
			if ($X0[$j] == 1 && $TSO eq "GGG"){
					$tc{$RNAME[$j]}{$currbar}{$UMI}++;
			}
		}
		$i++;
		if ($i == $x*1000000){
			print LOG $x." million reads processed\n";
			$x++;
		}	
	}
}

	
close(IN);
close(LOG);
close(BAR);


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
open(OUTT, ">", $coutt);
open(OUTB, ">", $coutb);
open(OUTC, ">", $coutc);

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
#	print $gene."\n";
	print OUTB $gene;
	print OUTT $gene;
	print OUTC $gene;
	foreach $cell (@cells){
		$n = 0;
		$rc = 0;
		foreach $umi (keys %{$tc{$gene}{$cell}}){
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
print SOUT "number of mapped reads with valid barcode: ".$trc."\n";
print SOUT "fraction of reads mapped with valid barcode: ".$trc/$i."\n";
close SOUT;



#open (BAROUT, ">", "mapped_bar.csv");
#foreach $barcode (keys %bar_map_h){
#	print BAROUT $barcode."\t".$bar_map_h{$barcode}."\n";
#}
#close BAROUT;




