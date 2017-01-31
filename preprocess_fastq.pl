#!/usr/bin/perl -w
# Usage: see usage message
use Getopt::Long;
use tools;
use strict;

my $version=getversion($0);
warn "Running $0, version $version\n";
warn "Arguments: @ARGV\n";

my($fastq, $umi_len, $cbc_len, $trim, $CBCbeforeUMI, $read1, $help);

$CBCbeforeUMI=0;                        # CELSeq2

## -xytrim worked, but didn't help so was removed at after commit 26547a0374 (2016-11-02 16:44:17)

my $usage = "
Usage: $0 --fastq s_R1.fastq.gz,s_R2.fastq.gz --umi_len 6 --cbc_len 8  [ OPTIONS ]  | gzip >  s_cbc.fastq.gz 

In CELSeq2, read1 contains (in that order) CBC, UMI, polyT, whereas read2
contains the mRNA.  This script takes the CBC and UMI from read1, and
appends e.g.':cbc=TACTGCTG:umi=GTCTTT' onto the read-id(previously they
were prepended to the read2 sequence). The resulting FASTQ file is written
to stdout.

The current protocols have an artefact that tends to produces long
stretches of polyA (and to a lesser extent polyT). Specifying
e.g. -trim=A=12,T=18 will delete any occurrence of AAAAAAAAAAAA.* and
TTTTTTTTTTTTTTTTTT.* from read2. The quality lines are trimmed in the same
way. (The numbers suggested in the usage message correspond to roughly
0.1% of the actual occurrences in the human transcriptome).

By default CELSeq2 is used, i.e. UMI precedes the cell bar code. Use -protocol=1
to swap them.

If also read1 should get its read-id changed, use the -read1 option; this will
write the amended reads to the (gzipped) fastq file (they will *not* be trimmed
if the -trim option is specified, because the meaningful  information,
if any, is beyond the polyT stretch). Note that the transcript is part of read2,
so this option is rarely needed.

Arguments:

    --fastq s_R1.fastq.gz,s_R2.fastq.gz # input files
    --umi_len=6      # length of the UMI
    --cbc_len=8      # length of the cell barcode

Options: 

    --read1 s_R1_cbc.fastqc.gz     # also save the read1's with the new id
    --trim=A18,T18                 # trim read2 of any stretches of 18 A's and 18 T's
    --CBCbeforeUMI                 # CELseq2 has first the  UMI, then the CBC, this option inverts that

Heavily adapted by <plijnzaad\@gmail.com> from the original written by Lennart Kester.
";

die $usage unless GetOptions('fastq=s'=> \$fastq,
                             'umi_len=i'=> \$umi_len,
                             'cbc_len=i'=> \$cbc_len,
                             'read1=s'=> \$read1,
                             'trim=s' => \$trim,
                             'CBCbeforeUMI' => \$CBCbeforeUMI,
                             'help|h' => \$help);
die $usage if $help;
die $usage unless $fastq && defined($umi_len) && defined($cbc_len);

my $regexps ={};

if (defined($trim)) { 
  my @nucs=split(',', $trim);
  for my $nt (@nucs) { 
    my($nuc, $num)= ($nt =~ /^([ACGT])(\d+)/);
    die "$0: expected string like -trim=A18,T18" unless $nuc && $num;
    $regexps->{$nuc}=$num;
  }

  for my $nuc ( keys  %$regexps ) { 
    my $re = '(' . $nuc x $regexps->{$nuc} . ".*)";
    $regexps->{$nuc}= qr/$re/;
  }
}

my $ntrimmed={};
my $ntrimmedtotal={};
my @all=sort keys %$regexps;
my @regexpids=( grep(/^[ATCG]$/, @all ), grep(/_/, @all )); # first -trim

for my $rid (keys @regexpids) {              # rid=regexp-id
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
  open(READ1, " | gzip > $read1.fastq.gz ") || die "read1: $!";
}

my ($line1, $line2, $bar);
my $trimmedlen={};

my (@lines1, @lines2);

my $nreads=0;
my $nemptyreads=0;

READ:
while( not eof $IN1 and not eof $IN2) {
  $trimmedlen={};
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

### do trimming, if any:
  my $line2=$lines2[1];
  for my $rid (@regexpids) { 
    if( $line2 =~ $regexps->{$rid} ) { 
      my $newlen=length($line2) - (length($1)+1); # +1 for the \n
      $trimmedlen->{$rid}=$newlen;    # remember for the qual line
      $line2= substr($line2,0, $newlen) . "\n";
      $ntrimmed->{$rid}++;
      $ntrimmedtotal->{$rid} += length($1);
    }
  }
  if ( $line2 eq "\n") { 
    $nemptyreads ++;
    next READ;
  }
  $lines2[1]=$line2;

### line with Phred qualities:
  $line2=$lines2[3];
  for my $rid (@regexpids) {               # trim qual line if seqline was
    if(exists($trimmedlen->{$rid})) { 
      $line2= substr($line2,0, $trimmedlen->{$rid}) . "\n";
    }
  }
  $lines2[3]=$line2;

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
