#!/usr/bin/perl -w
# Usage: see usage message
use Getopt::Long;
use tools;
use strict;

my $version=getversion($0);
warn "Running $0, version $version\n";
warn "Arguments: @ARGV\n";

my($fastq, $umi_len, $cbc_len, $polytrim, $CBCbeforeUMI, $read1, $five, $three, $help);

$CBCbeforeUMI=0;                        # CELSeq2

## -xytrim worked, but didn't help so was removed at after commit 26547a0374 (2016-11-02 16:44:17)

my $usage = "
Usage: $0 --fastq s_R1.fastq.gz,s_R2.fastq.gz --umi_len 6 --cbc_len 8  [ OPTIONS ]  | gzip -n >  s_cbc.fastq.gz 

In CELSeq2, read1 contains (in that order) CBC, UMI, polyT, whereas read2
contains the mRNA.  This script takes the CBC and UMI from read1, and
appends e.g.':cbc=TACTGCTG:umi=GTCTTT' onto the read-id. The resulting
FASTQ file is written to stdout (where it usually will be gzipped before
being written to disk)

The current protocols have an artefact that tends to produces long
stretches of polyA (and to a lesser extent polyT). Specifying
e.g. -trim=A=12,T=18 will delete any occurrence of AAAAAAAAAAAA.* and
TTTTTTTTTTTTTTTTTT.* from read2. The quality lines are trimmed in the same
way. (These numbers correspond to roughly 0.1% of the actual occurrences
in the human transcriptome).

By default CELSeq2 is used, i.e. UMI precedes the cell bar code. Use -protocol=1
to swap them.

If also read1 should get its read-id changed, use the -read1 option; this will
write the amended reads to the (gzipped) fastq file (they will *not* be trimmed
if the -trim option is specified, because the meaningful  information,
if any, is beyond the polyT stretch). Note that the transcript is part of read2,
so this option is prolly only needed when debugging script, lab protocol or both.

Arguments:

    --fastq s_R1.fastq.gz,s_R2.fastq.gz # input files
    --umi_len=6      # length of the UMI
    --cbc_len=8      # length of the cell barcode

Options: 

    --read1 s_R1_cbc.fastqc.gz     # also save the read1's with the new id
    --polytrim=G12,A14             # trim read2 of any stretches of 12 G's and 14 A's (in that order) and beyond
    --CBCbeforeUMI                 # CELseq2 has first the  UMI, then the CBC, this option inverts that
    --five=6                       # Trim 6 nt from the 5'-side of read2
    --three=8                      # Trim 8 nt from the 3'-side of read2 (only for reads that were not polytrimmed)

Heavily adapted by <plijnzaad\@gmail.com> from the original written by Lennart Kester.
";

die $usage unless GetOptions('fastq=s'=> \$fastq,
                             'umi_len=i'=> \$umi_len,
                             'cbc_len=i'=> \$cbc_len,
                             'read1=s'=> \$read1,
                             'polytrim=s' => \$polytrim,
                             'CBCbeforeUMI' => \$CBCbeforeUMI,
                             'five=i' => \$five,
                             'three=i' => \$three,
                             'help|h' => \$help);
die $usage if $help;
die $usage unless $fastq && defined($umi_len) && defined($cbc_len);

my $regexps ={};
my @regexpids = ();                     # to maintain the order

if (defined($polytrim)) { 
  my @oligos=split(',', $polytrim);
  for my $oli (@oligos) { 
    my($nuc, $num)= ($oli =~ /^([ACGT])(\d+)/);
    die "$0: expected string like --trim=A18,T18" unless $nuc && $num;
    my $re = '(' . $nuc x $num . ".*)";
    $regexps->{$nuc}= qr/$re/;
    push(@regexpids, $nuc);
  }
}

my $ntrimmed={};
my $ntrimmedtotal={};

for my $rid (@regexpids) {              # rid=regexp-id
  $ntrimmed->{$rid}=0;
  $ntrimmedtotal->{$rid}=0;
}

die "$0: no -umi_len specified" unless $umi_len > 0; # length of the UMI
die "$0: no -cbc_len specified" unless $cbc_len > 0; # length of the cell bar code

my $prefix_len = $cbc_len + $umi_len;

die "$fastq: input files must be gzipped " unless $fastq =~ /\.gz,.*\.gz$/;

my @fastq = split(/\,/,$fastq);

my($IN1, $IN2);
open($IN1, "zcat $fastq[0] |") || die "$0: $fastq[0]: $!";
open($IN2, "zcat $fastq[1] |") || die "$0: $fastq[1]: $!";

if($read1) { 
  $read1 =~ s/\.fastq.*$//i;
  open(READ1, " | gzip -n > $read1.fastq.gz ") || die "read1: $!";
}

my ($line1, $line2, $bar);
my $polytrimmedlen={};

my (@lines1, @lines2);

my $nreads=0;
my $nemptyreads=0;

READ:
while( not eof $IN1 and not eof $IN2) {
  $polytrimmedlen={};
  for(my $i=0; $i<4;$i++) {             # 4 lines at a time
    $lines1[$i] = <$IN1>;
    $lines2[$i] = <$IN2>; 
  }
  die "expected '+' lines in files @fastq, line $." 
      unless $lines1[2] eq "+\n" && $lines2[2] eq "+\n";

### id line:
  chomp($lines2[0]);
  my($id, $rest)=split(' ',$lines2[0]);
  my(undef, $rest1)=split(' ',$lines1[0]) if $read1;

### sequence line:
  $bar = substr($lines1[1], 0, $prefix_len);
  my $umi=substr($bar,0, $umi_len);
  my $cbc=substr($bar, $umi_len, $cbc_len);
  if ($CBCbeforeUMI) { 
    $cbc=substr($bar,0, $cbc_len);
    $umi=substr($bar, $cbc_len, $umi_len);
  }

  $lines1[0] = "$id:cbc=$cbc:umi=$umi $rest1\n" if $read1;
  $lines2[0] = "$id:cbc=$cbc:umi=$umi $rest\n";

  ## do polytrimming, if any:
  my $seq2=$lines2[1];
  chomp($seq2);
  for my $rid (@regexpids) { 
    if( $seq2 =~ $regexps->{$rid} ) { 
      my $trimmed=length($1);
      my $newlen=length($seq2) - $trimmed;
      $seq2= substr($seq2,0, $newlen);
      $polytrimmedlen->{$rid}=$newlen;    # remember for the qual line
      $ntrimmed->{$rid}++;
      $ntrimmedtotal->{$rid} += $trimmed;
    }
  }

  my $qual2=$lines2[3];
  chomp($qual2);
  # apply same trimming to phred qualities
  for my $rid (@regexpids) { 
    if(exists($polytrimmedlen->{$rid})) { 
      $qual2= substr($qual2,0, $polytrimmedlen->{$rid});
    }
  }

  ## ordinary trimming:
  if($five) {
    $seq2 = substr($seq2, $five);
    $qual2 = substr($qual2, $five);
  }

  if ($three && ! int(keys %$polytrimmedlen))  { 
    $seq2 = substr($seq2, 0, -$three);
    $qual2 = substr($qual2, 0, -$three);
  }

  if ( $seq2 eq "") {                   # trimmed to zilch
    $nemptyreads ++;
    next READ;
  }

  $lines2[1]=$seq2 ."\n";
  $lines2[3]=$qual2."\n";

  print  join("", @lines2);
  print  READ1 join("", @lines1) if $read1;

  $nreads++;
}                                       # READ

close $IN1 || die "$0: $fastq[0]: $!";
close $IN2 || die "$0: $fastq[1]: $!";

if ($read1) { 
  close (READ1) || die "$0: could not close (or open1) $read1: $!";
}

warn "preprocessed $nreads reads\n";
for my $rid (@regexpids) { 
  warn "trimmed $ntrimmed->{$rid} poly${rid}'s from the reads (totalling $ntrimmedtotal->{$rid} nucleotides)\n"
      if exists($ntrimmed->{$rid}) && $ntrimmed->{$rid} > 0;
}
warn "$nemptyreads of the trimmed reads were trimmed to length 0 and therefore discarded\n";
