#!/usr/bin/perl -w -s
use List::Util 'sum';

if (scalar @ARGV == 1){
    die "usage: -fastq=s_X_1.fastq -sam=sam.sam -barfile=barfile.csv -rb_len=rb_len -se=1 (for single end sequencing) -BCset=AA" if $ARGV[0] eq "help";
}


@cells = ();
@cbc = ();
$i = 0;
$x = 1;
%bar = ();
%tc = ();
$rb_len = 6;	
$se = 0 if !$se;

# open barcode file, read all 384 barcodes, combine them with the setbarcodes and put them in a hash, 
# the values are cell-numbers, 1-384 for the first set, 385-768 for the second set etc.;
$log = $sam;
$log   =~ s/(\.)\w+$/\.log/;
open(LOG, ">", $log);
print LOG "Process sam started\n";

open(IN,"<",$barfile);
while(<IN>){	
	chomp;
#	print $_."\n";
	@line = split(/\t/,$_);
	$bar{$line[1]} = $line[0];
#	print $line[1]."\t".$line[0]."\n";
	push(@cells, $line[1]);
}

close(IN);

$i = 0;
$j = 0;
$x = 1;
$tot_map_reads = 0;
$tot_map_reads_u = 0;

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
while( <IN> ){
	chomp $_;
	$line1 = $_;

#read through the header of the sam-file and barcode file;

	if (substr($line1,1,2) eq "SQ" ){

		next;
	}
	if (substr($line1,1,2) eq "PG" ){
		next;
	}

#	if ($i > 100000){
#		last;
#	}
	
# j = 0 for read 1, so here all the information is extracted for read 1;	
	if ($se == 0){
		if ($j == 0){
			@r1 = split(/\t/,$line1);
			($QNAME[$j],$FLAG[$j],$RNAME[$j],$POS[$j],$MAPQ[$j],$CIGAR[$j],$MRNM[$j],$MPOS[$j],
			$ISIZE[$j],$SEQ[$j],$QUAL[$j],$RBC)=split(/\t/,$line1);
			$cellbc = substr($RBC,($rb_len)+5,8);
#			print "R1\t".$SEQ[$j]."\n";
#			print "cellbc\t".$cellbc."\t";
			$UMI = substr($RBC,5,$rb_len);
#			print "UMI\t".$UMI."\n";
			$NM[$j] = "NA";
			$XA[$j] = "NA";
			$X0[$j] = 0;

			foreach $el (@r1){
				($dum,$dum,$NM[$j]) = split(/\:/,$el) if ($el =~ /^NM\:/);
				($dum,$dum,$XA[$j]) = split(/\:/,$el) if ($el =~ /^XA\:/);
				($dum,$dum,$X0[$j]) = split(/\:/,$el) if ($el =~ /^X0\:/);			
	
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
#			print "R2\t".$SEQ[$j]."\n";
#			print "GENEID\t".$RNAME[$j]."\n";
# check if the read fullfills the criteria in the if statements;
			if ($X0[$j] > 0){
				$tot_map_reads++;
			}
			if ($X0[$j] == 1){
				$tot_map_reads_u++;
			}		
			if (exists $bar{$cellbc}){
				if ($X0[$j] == 1 ){
#					if ($FLAG[1] == 153 or $FLAG[1] == 145){
						$tc{$RNAME[$j]}{$cellbc}{$UMI}++;
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
	
}

	
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
print SOUT "number of uniquely mapped reads: ".$tot_map_reads_u."\n";
print SOUT "fraction of reads mapped uniquely: ".$tot_map_reads_u/$i."\n";
print SOUT "number of mapped reads with valid barcode: ".$trc."\n";
print SOUT "fraction of reads mapped with valid barcode: ".$trc/$i."\n";
close SOUT;



#open (BAROUT, ">", "mapped_bar.csv");
#foreach $barcode (keys %bar_map_h){
#	print BAROUT $barcode."\t".$bar_map_h{$barcode}."\n";
#}
#close BAROUT;




