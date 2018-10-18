#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use IO::Uncompress::AnyUncompress;

#......................................................................................

our $VERSION = "0.4.2";
our $EXE = basename($0);
our $STDIN = '/dev/stdin'; # '-' will be replaced with this
our $URL = 'https://github.com/tseemann/any2fasta';

#......................................................................................

sub usage {
  my($errcode) = @_;
  $errcode ||= 0;
  my $ofh = $errcode ? \*STDERR : \*STDOUT;
  print $ofh 
    "NAME\n  $EXE $VERSION\n",
    "SYNOPSIS\n  Convert various sequence formats into FASTA\n",
    "USAGE\n  $EXE [options] file.{gb,fa,fq,gff,gfa,clw,sth}[.gz,bz2,zip] > output.fasta\n",
    "OPTIONS\n",
    "  -h       Print this help\n",
    "  -v       Print version and exit\n",
    "  -q       No output while running, only errors\n",
    "  -n       Replace ambiguous IUPAC letters with 'N'\n",
    "  -l       Lowercase the sequence\n",
    "  -u       Uppercase the sequence\n",
    "HOMEPAGE\n  $URL\n",
    "END\n";
  exit($errcode);
}

sub version {
  print "$EXE $VERSION\n";
  exit(0);
}

#......................................................................................

my %opt;
getopts('vhqnlu', \%opt) or exit(-1);
$opt{v} and version();
$opt{h} and usage(0);
@ARGV or usage(1);

#......................................................................................

sub msg {
  print STDERR "@_\n" unless $opt{q};
}

sub err {
  print STDERR "ERROR: @_\n";
  exit(-1);
}

#......................................................................................

msg("This is $EXE $VERSION");

# regexp to function mapping
my @formats = (
  [ 'GENBANK',   qr/^LOCUS\s/,       \&parse_genbank   ],
  [ 'EMBL',      qr/^ID\s/,          \&parse_embl      ],
  [ 'GFF',       qr/^##gff/,         \&parse_gff       ],
  [ 'CLUSTAL',   qr/^(CLUST|MUSCL)/, \&parse_clustal   ],
  [ 'STOCKHOLM', qr/^# STOCKHOLM\s/, \&parse_stockholm ],
  [ 'FASTA',     qr/^>/,             \&parse_fasta     ],
  [ 'FASTQ',     qr/^@/,             \&parse_fastq     ],
  [ 'GFA',       qr/^[A-Z]\t/,       \&parse_gfa       ],
);

# loop over all positional command line arguments
my $processed=0;

for my $fname (@ARGV) {
  msg("Opening '$fname'");
  $fname = $STDIN if $fname eq '-';
  -d $fname and err("'$fname' is a directory not a file");
  -r $fname or err("'$fname' is not readable");

  # read first line to see if we have any data
  my $unzip = IO::Uncompress::AnyUncompress->new($fname);
  my $header = scalar(<$unzip>); # read first line
  $header or err("The input appears to be empty");

  # detect format from first line
  my $ok=0;
  for my $fmt (@formats) {
    if ($header =~ $fmt->[1]) {
      msg("Detected", $fmt->[0], "format");
      # read in the rest of the file now
      my @line = ($header, <$unzip>);
      my $lines = scalar(@line);
      msg("Read $lines lines from '$fname'");
      # run the parser
      my $nseq = $fmt->[2]->( \@line );
      msg("Wrote $nseq sequences from", $fmt->[0], "file.");
      $nseq or err("No sequences found in '$fname'");
      $ok++;
      last;
    }
  }
  $ok or err("Unfamilar format with first line: $header");
  $processed++;
}

msg("Processed $processed files.");

msg("Done.");
exit(0);

#......................................................................................

sub purify {
  my($dna) = @_;
  $dna =~ s/[^ATGCN\n\r-]/N/gi if $opt{n};
  $dna = lc($dna) if $opt{l};
  $dna = uc($dna) if $opt{u};
  return $dna;
}

#......................................................................................

sub parse_fasta {
  my($lines) = @_;
  my $count=0;
  for my $line (@$lines) {
    next if $line =~ m/^\s*$/;
    if ($line =~ m/^>/) {
      $count++;
    }
    else {
      $line = purify($line);
    }
    print $line;
  }
  return $count;
}

#......................................................................................

sub parse_fastq {
  my($lines) = @_;
  my $count=0;
  # jump 4 lines at a time
  for ( my $i=0 ; $i < $#{$lines} ; $i+=4 ) {
    print ">", substr($lines->[$i], 1);
    print purify( $lines->[$i+1] );
    $count++;
  }
  return $count;
}

#......................................................................................

sub parse_gff {
  my($lines) = @_;
  my $at_seq = 0;
  for my $line (@$lines) {
    $at_seq++ if $line =~ m/^>/;
    print purify($line) if $at_seq;;
  }
  return $at_seq;
}

#......................................................................................
# https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md 
# https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md 

sub parse_gfa {
  my($lines) = @_;
  my %S;
  my %P;
  my $count=0;
  for my $line (@$lines) {
    my(@x) = split m/\t/, $line;
    # this is NOT the original contigs, rather the unitigs
    # need to parse L (link) and P (path) to reconstruct the contigs
    if ($x[0] eq 'S') {
      print ">", $x[1], "\n", purify($x[2]), "\n";
      $count++;
    }
  }
  return $count;
}

#......................................................................................

sub parse_genbank {
  my($lines) = @_;
  my $acc = '';
  my $dna = '';
  my $in_seq = 0;
  my $count = 0;

  foreach (@$lines) {
    chomp;
    if (m{^//}) {
      $dna = purify($dna);
      print ">", $acc, "\n", purify($dna);
      $count++;
      $in_seq = 0;
      $dna = $acc = '';
      next;
    }
    elsif (m/^ORIGIN/) {
      # ORIGIN
      $in_seq = 1;
      next;
    }
    
    if ($in_seq) {
      #       421 ctctcaaact aaagccgtct cactctccat gagtcgttcg acagatcgcg ttttaaattg
      my $s = substr $_, 10;  # trim the coordinate prefix
      $s =~ s/\s//g;
      $dna .= $s . "\n";
    }
    else {
      # LOCUS  NZ_AHMY02000075  683 bp   DNA   linear   CON 23-NOV-2017
      if (m/^LOCUS\s+(\S+)/) {
        $acc = $1;
      }
    }
  }
  return $count;
}

#......................................................................................

sub parse_embl {
  my($lines) = @_;
  my $acc = '';
  my $in_seq = 0;
  my $dna = '';
  my $count = 0;

  foreach (@$lines) {
    chomp;
    if (m{^//}) {
      # end of record
      print ">", $acc, "\n", purify($dna);
      $count++;
      $in_seq = 0;
      $dna = $acc = '';
      next;
    }
    elsif (m/^SQ\s/) {
      # SQ   Sequence 569 BP; 145 A; 133 C; 152 G; 139 T; 0 other;
      $in_seq = 1;
      next;
    }

    if ($in_seq) {
      #      agtcgctttt aattatgaat gttgtaacta cattatcatc gctgtcagtc ttctggctgg        60
      s/[\s\d]//g;
      $dna .= $_ . "\n";
    }
    else {
      # ID   K02675; SV 1; linear; genomic DNA; STD; UNC; 569 BP.
      if (m/^ID\s+([^;]+)/) {
        $acc = $1;
      }
    }
  }
  return $count;
}

#......................................................................................

sub parse_clustal {
  my($lines) = @_;
  my %seq;
  foreach (@$lines) {
    next unless m/^(\S+)\s+([A-Z-]+)$/i; # uses '-' for gap
    $seq{$1} .= $2."\n";
  }
  for my $id (sort keys %seq) {
    print ">", $id, "\n", purify($seq{$id});
  }
  return scalar(keys %seq);
}

#......................................................................................

sub parse_stockholm {
  my($lines) = @_;
  my %seq;
  foreach (@$lines) {
    next if m/^#/;
    last if m{^//};
    next unless m/^(\S+)\s+([A-Z.]+)$/i;  # uses '.' for gap
    my($id,$sq) = ($1,$2);
    $sq =~ s/\./-/g;  # switch to standard FASTA '-' gap char
    $seq{$id} .= $sq . "\n";
  }
  for my $id (sort keys %seq) {
    print ">", $id, "\n", purify($seq{$id});
  }
  return scalar(keys %seq);
}

#......................................................................................
