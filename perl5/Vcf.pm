package Vcf;

our $VERSION = 'r953';

# http://vcftools.sourceforge.net/specs.html
# http://samtools.github.io/hts-specs/
#
# Authors: petr.danecek@sanger
# for VCF v3.2, v3.3, v4.0, v4.1, v4.2
#

=head1 NAME

Vcf.pm.  Module for validation, parsing and creating VCF files. 
         Supported versions: 3.2, 3.3, 4.0, 4.1, 4.2

=head1 SYNOPSIS

From the command line:
    perl -MVcf -e validate example.vcf
    perl -I/path/to/the/module/ -MVcf -e validate_v32 example.vcf

From a script:
    use Vcf;

    my $vcf = Vcf->new(file=>'example.vcf.gz',region=>'1:1000-2000');
    $vcf->parse_header();

    # Do some simple parsing. Most thorough but slowest way how to get the data.
    while (my $x=$vcf->next_data_hash()) 
    { 
        for my $gt (keys %{$$x{gtypes}})
        {
            my ($al1,$sep,$al2) = $vcf->parse_alleles($x,$gt);
            print "\t$gt: $al1$sep$al2\n";
        }
        print "\n";
    }

    # This will split the fields and print a list of CHR:POS
    while (my $x=$vcf->next_data_array()) 
    {
        print "$$x[0]:$$x[1]\n";
    }

    # This will return the lines as they were read, including the newline at the end
    while (my $x=$vcf->next_line()) 
    { 
        print $x;
    }

    # Only the columns NA00001, NA00002 and NA00003 will be printed.
    my @columns = qw(NA00001 NA00002 NA00003);
    print $vcf->format_header(\@columns);
    while (my $x=$vcf->next_data_array())
    {
        # this will recalculate AC and AN counts, unless $vcf->recalc_ac_an was set to 0
        print $vcf->format_line($x,\@columns); 
    }

    $vcf->close();

=cut


use strict;
use warnings;
use Carp;
use Exporter;
use Data::Dumper;
use POSIX ":sys_wait_h";

use vars qw/@ISA @EXPORT/;
@ISA = qw/Exporter/;
@EXPORT = qw/validate validate_v32/;

=head2 validate

    About   : Validates the VCF file.
    Usage   : perl -MVcf -e validate example.vcf.gz     # (from the command line)
              validate('example.vcf.gz');               # (from a script)
              validate(\*STDIN);
    Args    : File name or file handle. When no argument given, the first command line
              argument is interpreted as the file name.

=cut

sub validate
{
    my ($fh) = @_;

    if ( !$fh && @ARGV ) { $fh = $ARGV[0]; }

    my $vcf;
    if ( $fh ) { $vcf = fileno($fh) ? Vcf->new(fh=>$fh) : Vcf->new(file=>$fh); }
    else { $vcf = Vcf->new(fh=>\*STDIN); }

    $vcf->run_validation();
}


=head2 validate_v32

    About   : Same as validate, but assumes v3.2 VCF version.
    Usage   : perl -MVcf -e validate_v32 example.vcf.gz     # (from the command line)
    Args    : File name or file handle. When no argument given, the first command line
              argument is interpreted as the file name.

=cut

sub validate_v32
{
    my ($fh) = @_;

    if ( !$fh && @ARGV && -e $ARGV[0] ) { $fh = $ARGV[0]; }

    my %params = ( version=>'3.2' );

    my $vcf;
    if ( $fh ) { $vcf = fileno($fh) ? Vcf->new(%params, fh=>$fh) : Vcf->new(%params, file=>$fh); }
    else { $vcf = Vcf->new(%params, fh=>\*STDIN); }

    $vcf->run_validation();
}


=head2 new

    About   : Creates new VCF reader/writer. 
    Usage   : my $vcf = Vcf->new(file=>'my.vcf', version=>'3.2');
    Args    : 
                fh      .. Open file handle. If neither file nor fh is given, open in write mode.
                file    .. The file name. If neither file nor fh is given, open in write mode.
                region  .. Optional region to parse (requires tabix indexed VCF file)
                silent  .. Unless set to 0, warning messages may be printed.
                strict  .. Unless set to 0, the reader will die when the file violates the specification.
                version .. If not given, '4.0' is assumed. The header information overrides this setting.

=cut

sub new
{
    my ($class,@args) = @_;
    my $self = {@args};
    bless $self, ref($class) || $class;

    $$self{silent}    = 0 unless exists($$self{silent});
    $$self{strict}    = 0 unless exists($$self{strict});
    $$self{buffer}    = [];       # buffer stores the lines in the reverse order
    $$self{columns}   = undef;    # column names 
    $$self{mandatory} = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'] unless exists($$self{mandatory}); 
    $$self{reserved}{cols} = {CHROM=>1,POS=>1,ID=>1,REF=>1,ALT=>1,QUAL=>1,FILTER=>1,INFO=>1,FORMAT=>1} unless exists($$self{reserved_cols});
    $$self{recalc_ac_an} = 1;
    $$self{has_header} = 0;
    $$self{default_version} = '4.2';
    $$self{versions} = [ qw(Vcf3_2 Vcf3_3 Vcf4_0 Vcf4_1 Vcf4_2) ];
    if ( !exists($$self{max_line_len}) && exists($ENV{MAX_VCF_LINE_LEN}) ) { $$self{max_line_len} = $ENV{MAX_VCF_LINE_LEN} }
    $$self{fix_v40_AGtags} = $ENV{DONT_FIX_VCF40_AG_TAGS} ? 0 : 1;
    my %open_args = ();
    if ( exists($$self{region}) ) 
    { 
        $open_args{region}=$$self{region}; 
        if ( !exists($$self{print_header}) ) { $$self{print_header}=1; }
    }
    if ( exists($$self{print_header}) ) { $open_args{print_header}=$$self{print_header}; }
    return $self->_open(%open_args);
}

sub throw
{
    my ($self,@msg) = @_;
    confess @msg,"\n";
}

sub warn
{
    my ($self,@msg) = @_;
    if ( $$self{silent} ) { return; }
    if ( $$self{strict} ) { $self->throw(@msg); }
    warn @msg;
}

sub _open
{
    my ($self,%args) = @_;

    if ( !exists($$self{fh}) && !exists($$self{file}) ) 
    {
        # Write mode, the version must be supplied by the user 
        return $self->_set_version(exists($$self{version}) ? $$self{version} : $$self{default_version});
    }

    # Open the file unless filehandle is provided
    if ( !exists($$self{fh}) ) 
    {
        if ( !defined $$self{file} ) { $self->throw("Undefined value passed to Vcf->new(file=>undef)."); }
        my $cmd = "<$$self{file}";

        my $tabix_args = '';
        if ( exists($args{print_header}) && $args{print_header} ) { $tabix_args .= ' -h '; }
        $tabix_args .= qq['$$self{file}'];
        if ( exists($args{region}) && defined($args{region}) ) { $tabix_args .= qq[ '$args{region}']; }

        if ( -e $$self{file} && $$self{file}=~/\.gz/i )
        {
            if ( exists($args{region}) && defined($args{region}) )
            {
                $cmd = "tabix $tabix_args |";
            }
            else { $cmd = "gunzip -c '$$self{file}' |"; } 
        }
        elsif ( $$self{file}=~m{^(?:http|ftp)://} )
        {
            if ( !exists($args{region}) ) { $tabix_args .= ' .'; }
            $cmd = "tabix $tabix_args |";
        }
        open($$self{fh},$cmd) or $self->throw("$cmd: $!");
    }

    # Set the correct VCF version, but only when called for the first time
    my $vcf = $self;
    if ( !$$self{_version_set} ) 
    { 
        my $first_line = $self->next_line();
        $vcf = $self->_set_version($first_line);
        $self->_unread_line($first_line);
    }
    return $vcf;
}



=head2 open

    About   : (Re)Open file. No need to call this explicitly unless reading from a different 
              region is requested.
    Usage   : $vcf->open(); # Read from the start
              $vcf->open(region=>'1:12345-92345');
    Args    : region       .. Supported only for tabix indexed files

=cut

sub open
{
    my ($self,%args) = @_;
    $self->close();
    $self->_open(%args);
}


=head2 close

    About   : Close the filehandle
    Usage   : $vcf->close();
    Args    : none
	Returns : close exit status

=cut

sub close
{
    my ($self) = @_;
    if ( !$$self{fh} ) { return; }
    my $ret = close($$self{fh});
    delete($$self{fh});
    $$self{buffer} = [];
	return $ret;
}


=head2 next_line

    About   : Reads next VCF line.
    Usage   : my $vcf = Vcf->new(); 
              my $x   = $vcf->next_line();
    Args    : none

=cut

sub next_line
{
    my ($self) = @_;
    if ( @{$$self{buffer}} ) { return shift(@{$$self{buffer}}); }

    my $line;
    if ( !exists($$self{max_line_len}) ) 
    {
        $line = readline($$self{fh});
    }
    else
    {
        while (1)
        {
            $line = readline($$self{fh});
            if ( !defined $line ) { last; }

            my $len = length($line);
            if ( $len>$$self{max_line_len} && !($line=~/^#/) ) 
            { 
                if ( !($line=~/^([^\t]+)\t([^\t]+)/) ) { $self->throw("Could not parse the line: $line"); }
                $self->warn("The VCF line too long, ignoring: $1 $2 .. len=$len\n"); 
                next;
            }
            last;
        }
    }
    return $line;
}

sub _unread_line
{
    my ($self,$line) = @_;
    unshift @{$$self{buffer}}, $line;
    return;
}


=head2 next_data_array

    About   : Reads next VCF line and splits it into an array. The last element is chomped.
    Usage   : my $vcf = Vcf->new(); 
              $vcf->parse_header(); 
              my $x = $vcf->next_data_array();
    Args    : Optional line to parse

=cut

sub next_data_array
{
    my ($self,$line) = @_;
    if ( !$line ) { $line = $self->next_line(); }
    if ( !$line ) { return undef; }
    if ( ref($line) eq 'ARRAY' ) { return $line; }
    my @items = split(/\t/,$line);
    if ( @items<8 ) { $line=~s/\n/\\n/g; $self->throw("Could not parse the line, wrong number of columns: [$line]"); }
    chomp($items[-1]);
    return \@items;
}


=head2 set_samples

    About   : Parsing big VCF files with many sample columns is slow, not parsing unwanted samples may speed things a bit.
    Usage   : my $vcf = Vcf->new(); 
              $vcf->set_samples(include=>['NA0001']);   # Exclude all but this sample. When the array is empty, all samples will be excluded.
              $vcf->set_samples(exclude=>['NA0003']);   # Include only this sample. When the array is empty, all samples will be included.
              my $x = $vcf->next_data_hash();
    Args    : Optional line to parse

=cut

sub set_samples
{
    my ($self,%args) = @_;

    if ( exists($args{include}) )
    {
        for (my $i=0; $i<@{$$self{columns}}; $i++) { $$self{samples_to_parse}[$i] = 0; }
        for my $sample (@{$args{include}})
        {
            if ( !exists($$self{has_column}{$sample}) ) { $self->throw("The sample not present in the VCF file: [$sample]\n"); }
            my $idx = $$self{has_column}{$sample} - 1;
            $$self{samples_to_parse}[$idx]  = 1;
        }
    }
    
    if ( exists($args{exclude}) )
    {
        for (my $i=0; $i<@{$$self{columns}}; $i++) { $$self{samples_to_parse}[$i] = 1; }
        for my $sample (@{$args{exclude}})
        {
            if ( !exists($$self{has_column}{$sample}) ) { $self->throw("The sample not present in the VCF file: [$sample]\n"); }
            my $idx = $$self{has_column}{$sample} - 1;
            $$self{samples_to_parse}[$idx]  = 0;
        }
    }
}


sub _set_version
{
    my ($self,$version_line) = @_;

    if ( $$self{_version_set} ) { return $self; }
    $$self{_version_set} = 1;

    $$self{version} = $$self{default_version};
    if ( $version_line )
    {
        if ( $version_line=~/^(\d+(?:\.\d+)?)$/ )
        {
            $$self{version} = $1;
            undef $version_line;
        }
        elsif ( !($version_line=~/^##fileformat=/i) or !($version_line=~/(\d+(?:\.\d+)?)\s*$/i) ) 
        { 
			chomp($version_line);
            $self->warn("Could not parse the fileformat version string [$version_line], assuming VCFv$$self{default_version}\n"); 
            undef $version_line;
        }
        else
        {
            $$self{version} = $1;
        }
    }

    my $reader;
    if ( $$self{version} eq '3.2' ) { $reader=Vcf3_2->new(%$self); }
    elsif ( $$self{version} eq '3.3' ) { $reader=Vcf3_3->new(%$self); } 
    elsif ( $$self{version} eq '4.0' ) { $reader=Vcf4_0->new(%$self); }
    elsif ( $$self{version} eq '4.1' ) { $reader=Vcf4_1->new(%$self); }
    elsif ( $$self{version} eq '4.2' ) { $reader=Vcf4_2->new(%$self); }
    else 
    { 
        $self->warn(qq[The version "$$self{version}" not supported, assuming VCFv$$self{default_version}\n]);
        $$self{version} = '4.2';
        $reader = Vcf4_2->new(%$self);
    }

    $self = $reader;
    # When changing version, change also the fileformat header line
    if ( exists($$self{header_lines}) && exists($$self{header_lines}[0]{key}) && $$self{header_lines}[0]{key} eq 'fileformat' )
    {
        shift(@{$$self{header_lines}});
    }

    return $self;
}


#---------------------------------------

package VcfReader;
use base qw(Vcf);
use strict;
use warnings;
use Carp;
use Data::Dumper;

sub new
{
    my ($class,@args) = @_;
    my $self = {@args};
    bless $self, ref($class) || $class;
    return $self;
}


=head2 next_data_hash

    About   : Reads next VCF line and splits it into a hash. This is the slowest way to obtain the data.
    Usage   : my $vcf = Vcf->new(); 
              $vcf->parse_header(); 
              my $x = $vcf->next_data_hash();

              # Or having a VCF data line $line
              my $x = $vcf->next_data_hash($line);

    Args    : Optional line to parse.

=cut

sub next_data_hash
{
    my ($self,$line) = @_;
    if ( !$line ) { $line = $self->next_line(); }
    if ( !$line ) { return undef; }
    my @items;
    if ( ref($line) eq 'ARRAY' ) { @items = @$line; }
    else { @items = split(/\t/,$line); }
    chomp($items[-1]);

    my $cols = $$self{columns};
    if ( !$cols ) 
    { 
        $self->_fake_column_names(scalar @items - 9); 
        $cols = $$self{columns};
    }

    # Check the number of columns
    if ( scalar @items != scalar @$cols )  
    { 
        if ( $line=~/^\s*$/ ) { $self->throw("Sorry, empty lines not allowed.\n"); }
        my $c = substr($line,0,1);
        if ( $c eq '#' ) 
        { 
            if ( !$$self{header_parsed} ) { $self->throw("FIXME: parse_header must be called before next_data_hash.\n"); }
            else { $self->throw("Multiple header blocks (^#) not allowed.\n"); }
        }

        if ( $items[-1] eq '' )
        {
            my $nremoved = 0;
            while ( $items[-1] eq '' ) { pop(@items); $nremoved++; }
            if ( $nremoved && !$$self{trailing_tabs_warned} )
            {
                $self->warn("Broken VCF: empty columns (trailing TABs) starting at $items[0]:$items[1].\n");
                $$self{trailing_tabs_warned} = 1;
            }
        }
        if ( scalar @items != scalar @$cols ) 
        {
            my @test = split(/\s+/,$line);
            if ( scalar @test == scalar @$cols ) { $self->warn("(Were spaces used instead of tabs?)\n\n"); }
            else { $self->throw(sprintf "Wrong number of fields%s; expected %d, got %d. The offending line was:\n[%s]\n\n", 
                exists($$self{file}) ? "in $$self{file}" : '', scalar @$cols, scalar @items, join("\t",@items)); }

            @items = @test;
        }
    }
    my %out;

    # Mandatory fields
    $out{CHROM}  = $items[0];
    $out{POS}    = $items[1];
    $out{ID}     = $items[2];
    $out{REF}    = $items[3];
    $out{ALT}    = [ split(/,/,$items[4]) ];
    $out{QUAL}   = $items[5];
    $out{FILTER} = [ split(/;/,$items[6]) ];

    # INFO, e.g. NS=58;DP=258;AF=0.786;DB;H2
    if ( defined $items[7] )
    {
        my %hash;
        for my $info (split(/;/,$items[7]))
        {
            my ($key,$val) = split(/=/,$info);
            if ( !defined $key ) 
            { 
                $self->warn("Broken VCF file, empty INFO field at $items[0]:$items[1]\n");
                next;
            }
            if ( defined $val )
            {
                $hash{$key} = $val;
            }
            elsif ( exists($$self{header}{INFO}{$key}) )
            {
                $hash{$key} = $$self{header}{INFO}{$key}{default};
            }
            else
            {
                $hash{$key} = undef;
            }
        }
        $out{INFO} = \%hash;
    }

    # The FORMAT field may not be present. GT:GQ:DP:HQ
    my $format;
    if ( $$cols[8] || $items[8] )
    {
        $format = $out{FORMAT} = [ split(/:/,$items[8]) ];
        if ( (!$$format[0] || $$format[0] ne 'GT') && !$$self{ignore_missing_GT} ) { $self->warn("Expected GT as the first genotype field at $items[0]:$items[1]\n"); } 
    }

    # Genotype fields
    my %gtypes;
    my $check_nformat = $$self{drop_trailings} ? 0 : 1;
    for (my $icol=9; $icol<@items; $icol++)
    {
        if ( $items[$icol] eq '' ) { $self->warn("Empty column $$cols[$icol] at $items[0]:$items[1]\n"); next; }
        if ( exists($$self{samples_to_parse}) && !$$self{samples_to_parse}[$icol] ) { next; }

        my @fields = split(/:/, $items[$icol]);
        if ( $check_nformat && @fields != @$format ) 
        {
            $self->warn("Different number of fields in the format and the column $$cols[$icol] at $items[0]:$items[1] ("
                .scalar @fields." vs ".scalar @$format.": [",join(',',@fields),"] vs [",join(',',@$format),"])\n");
        }
        my %hash;
        for (my $ifield=0; $ifield<@fields; $ifield++)
        {
            $hash{$$format[$ifield]} = $fields[$ifield];
        }
        $gtypes{$$cols[$icol]} = \%hash;
    }
    $out{gtypes} = \%gtypes;

    return \%out;
}


=head2 parse_header

    About   : Reads (and stores) the VCF header.
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header();
    Args    : silent .. do not warn about duplicate header lines

=cut

sub parse_header
{
    my ($self,%args) = @_;

    # First come the header lines prefixed by ##
    while ($self->_next_header_line(%args)) { ; }

    # Now comes the column names line prefixed by #
    $self->_read_column_names();

    $$self{header_parsed} = 1;
}


=head2 _next_header_line

    About   : Stores the header lines and meta information, such as fields types, etc.
    Args    : silent .. do not warn about duplicate column names

=cut

sub _next_header_line
{
    my ($self,%args) = @_;
    my $line = $self->next_line();
    if ( !defined $line ) { return undef; }
    if ( substr($line,0,2) ne '##' )
    {
        $self->_unread_line($line);
        return undef;
    }

    my $rec = $self->parse_header_line($line);
    if ( $rec ) { $self->add_header_line($rec,%args); }

    return $rec;
}

=head2 get_header_line

    Usage   : $vcf->get_header_line(key=>'INFO', ID=>'AC')
              $vcf->get_header_line(key=>'FILTER', ID=>'q10')
              $vcf->get_header_line(key=>'reference')
              $vcf->get_header_line(key=>'contig',ID=>'20')
    Args    : Header line filter as in the example above
    Returns : List ref of header line hashes matching the filter

=cut

sub get_header_line
{
    my ($self,%filter) = @_;

    my $key = $filter{key};
    delete($filter{key});

    my $id = $filter{ID};

    my @out;
    while (my ($hline_key,$hline_hash) = each %{$$self{header}})
    {
        if ( $key ne $hline_key ) { next; }

        if ( defined $id ) 
        { 
            if ( !exists($$hline_hash{$id}) ) { next; }
            $hline_hash = $$hline_hash{$id};
        }

        my $match = 1;
        while (my ($fkey,$fval) = each %filter)
        {
            if ( !exists($$hline_hash{$fkey}) or $$hline_hash{$fkey} ne $fval ) 
            { 
                $match=0; 
                last; 
            }
        }
        if ( $match ) { push @out,$hline_hash }
    }
    return \@out;
}


=head2 add_header_line

    Usage   : $vcf->add_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'})
              $vcf->add_header_line({key=>'reference',value=>'1000GenomesPilot-NCBI36'})
    Args    : Header line hash as in the example above
              Hash with additional parameters [optional]
                silent .. do not warn about existing header keys
                append .. append timestamp to the name of the new one
    Returns : 

=cut

sub add_header_line
{
    my ($self,$rec,%args) = @_;

    if ( !%args ) { $args{silent}=0; }

    my $key = $$rec{key};
    if ( !$key ) { $self->throw("Missing key: ",Dumper($rec)); }

    if ( exists($$rec{Type}) )
    {
        if ( !exists($$rec{default}) )
        {
            my $type = $$rec{Type};
            if ( exists($$self{defaults}{$type}) ) { $$rec{default}=$$self{defaults}{$type}; }
            else { $$rec{default}=$$self{defaults}{default}; }
        }
        if ( !exists($$rec{handler}) )
        {
            my $type = $$rec{Type};
            if ( !exists($$self{handlers}{$type}) ) 
            { 
                $self->warn("Unknown type [$type]\n"); 
                $type = 'String';
                $$rec{Type} = $type;
            }
            if ( exists($$self{handlers}{$type}) ) { $$rec{handler}=$$self{handlers}{$type}; }
            else { $self->throw("Unknown type [$type].\n"); }
        }
    }

    if ( exists($$rec{ID}) )
    {
        my $id = $$rec{ID};
        if ( exists($$self{header}{$key}{$id}) ) { $self->remove_header_line(%$rec); }
        $$self{header}{$key}{$id} = $rec;
        push @{$$self{header_lines}}, $rec;
        return;
    }

    if ( $args{append} )
    {
        my @tm = gmtime(time);
        $key = sprintf "%s_%d%.2d%.2d", $key,$tm[5]+1900,$tm[4]+1,$tm[3];
        my $i = 1;
        while ( exists($$self{header}{$key.'.'.$i}) ) { $i++; }
        $key = $key.'.'.$i;
        $$rec{key} = $key;
    }

    if ( $self->_header_line_exists($key,$rec) ) { $self->remove_header_line(%$rec); }

    push @{$$self{header}{$key}}, $rec;
    if ( $$rec{key} eq 'fileformat' ) 
    { 
        unshift @{$$self{header_lines}}, $rec; 
    }
    else
    {
        push @{$$self{header_lines}}, $rec;
    }
}

sub _header_line_exists
{
    my ($self,$key,$rec) = @_;
    if ( !exists($$self{header}{$key}) ) { return 0; }
    if ( $key eq 'fileformat' ) { return 1; }
    for my $hrec (@{$$self{header}{$key}})
    {
        my $differ = 0;
        for my $item (keys %$rec)
        {
            if ( !exists($$hrec{$item}) ) { $differ=1; last; }
            if ( $$hrec{$item} ne $$rec{$item} ) { $differ=1; last; }
        }
        if ( !$differ ) { return $hrec; }
    }
    return 0;
}

=head2 remove_header_line

    Usage   : $vcf->remove_header_line(key=>'INFO', ID=>'AC')
    Args    :
    Returns : 

=cut

sub remove_header_line
{
    my ($self,%args) = @_;
    my $key = $args{key};
    my %to_be_removed;
    for (my $i=0; $i<@{$$self{header_lines}}; $i++)
    {
        my $line = $$self{header_lines}[$i];
        if ( $$line{key} ne $key ) { next; }
        if ( exists($args{ID}) )
        {
            if ( $args{ID} ne $$line{ID} ) { next; }
            delete($$self{header}{$key}{$args{ID}});
            splice(@{$$self{header_lines}},$i--,1);
        }
        elsif ( scalar keys %args==1 && exists($$self{header}{$key}) )
        {
            splice(@{$$self{header_lines}},$i--,1);
            $to_be_removed{$key} = 1;
        }
        else
        {
            my $to_be_removed = $self->_header_line_exists($key,\%args);
            if ( !$to_be_removed ) { next; }
            for (my $j=0; $j<@{$$self{header}{$key}}; $j++)
            {
                if ( $$self{header}{$key}[$j] eq $to_be_removed ) { splice(@{$$self{header}{$key}},$j,1); last; }
            }
            splice(@{$$self{header_lines}},$i--,1);
        }
    }
    for my $key (keys %to_be_removed) { delete($$self{header}{$key}); }
}


=head2 parse_header_line

    Usage   : $vcf->parse_header_line(q[##reference=1000GenomesPilot-NCBI36])
              $vcf->parse_header_line(q[##INFO=NS,1,Integer,"Number of Samples With Data"])
    Args    : 
    Returns : 

=cut

sub parse_header_line
{
    my ($self,$line) = @_;

    chomp($line);
    $line =~ s/^##//;
    
    if ( !($line=~/^([^=]+)=/) ) { return { key=>$line, value=>'' }; }
    my $key   = $1;
    my $value = $';

    my $desc;
    if ( $value=~/,\s*\"([^\"]+)\"\s*$/ ) { $desc=$1; $value=$`; }

    if ( !$desc ) { return { key=>$key, value=>$value }; }

    if ( $key eq 'INFO' or $key eq 'FORMAT' )
    {
        my ($id,$number,$type,@rest) = split(/,\s*/,$value);
        if ( !$type or scalar @rest ) { $self->throw("Could not parse the header line: $line\n"); }
        return { key=>$key, ID=>$id, Number=>$number, Type=>$type, Description=>$desc };
    }
    if ( $key eq 'FILTER' )
    {
        my ($id,@rest) = split(/,\s*/,$value);
        if ( !$id or scalar @rest ) { $self->throw("Could not parse the header line: $line\n"); }
        return { key=>$key, ID=>$id, Description=>$desc };
    }
    $self->throw("Could not parse the header line: $line\n"); 
}

=head2 _read_column_names

    About   : Stores the column names as array $$self{columns} and hash $$self{has_column}{COL_NAME}=index.
              The indexes go from 1.
    Usage   : $vcf->_read_column_names();
    Args    : none

=cut

sub _read_column_names
{
    my ($self) = @_;
    my $line = $self->next_line();
    if ( !defined $line or substr($line,0,1) ne '#' ) { $self->throw("Broken VCF header, no column names?"); }
    $$self{column_line} = $line;

    my @cols  = split(/\t/, substr($line,1));
    chomp($cols[-1]);

    my $nremoved = 0;
    for (my $i=0; $i<@cols; $i++)
    {
        if ( !($cols[$i]=~/^\s*$/) ) { next; }
        $self->warn(sprintf "Empty fields in the header line, the column %d is empty, removing.\n",$i+1+$nremoved);
        $nremoved++;
        splice(@cols,$i,1);
    }

    my $ncols = scalar @cols;
    if ( $ncols == 1 )
    {
        # If there is only one name, it can be space-separated instead of tab separated
        @cols  = split(/\s+/, $cols[0]);
        $ncols = scalar @cols;
        chomp($line);
        if ( $ncols <= 1 ) { $self->warn("Could not parse the column names. [$line]\n"); return; }
        $self->warn("The column names not tab-separated? [$line]\n");
    }

    my $fields  = $$self{mandatory};
    my $nfields = scalar @$fields;

    # Check the names of the mandatory columns
    if ( $ncols < $nfields ) 
    { 
        chomp($line); 
        $self->warn("Missing some of the mandatory column names.\n\tGot:      $line\n\tExpected: #", join("\t",@{$$self{mandatory}}),"\n"); 
        return; 
    }

    for (my $i=0; $i<$ncols; $i++)
    {
        if ( $cols[$i]=~/^\s+/ or $cols[$i]=~/\s+$/ ) 
        {
            $self->warn("The column name contains leading/trailing spaces, removing: '$cols[$i]'\n");
            $cols[$i] =~ s/^\s+//;
            $cols[$i] =~ s/\s+$//;
        }
        if ( $i<$nfields && $cols[$i] ne $$fields[$i] ) 
        { 
            $self->warn("Expected mandatory column [$$fields[$i]], got [$cols[$i]]\n"); 
            $cols[$i] = $$fields[$i];
        }
        $$self{has_column}{$cols[$i]} = $i+1;
    }
    $$self{columns} = \@cols;
    return;
}


=head2 _fake_column_names

    About   : When no header is present, fake column names as the default mandatory ones + numbers
    Args    : The number of genotype columns; 0 if no genotypes but FORMAT present; <0 if FORMAT and genotypes not present

=cut

sub _fake_column_names
{
    my ($self,$ncols) = @_;

    $$self{columns} = [ @{$$self{mandatory}} ];
    if ( $ncols>=0 ) { push @{$$self{columns}}, 'FORMAT'; }
    for (my $i=1; $i<=$ncols; $i++) { push @{$$self{columns}}, $i; }
}


=head2 format_header

    About   : Returns the header.
    Usage   : print $vcf->format_header();
    Args    : The columns to include on output [optional]

=cut

sub format_header
{
    my ($self,$columns) = @_;

    my $out = '';
    for my $line (@{$$self{header_lines}}) { $out .= $self->format_header_line($line); }

    # This is required when using the API for writing new VCF files and the caller does not add the line explicitly
    if ( !exists($$self{header_lines}[0]{key}) or $$self{header_lines}[0]{key} ne 'fileformat' ) 
    { 
        $out = "##fileformat=VCFv$$self{version}\n" .$out; 
    }
    if ( !$$self{columns} ) { return $out; }

    my @out_cols;
    if ( $columns )
    {
        @out_cols = @{$$self{columns}}[0..8];
        for my $col (@$columns)
        {
            if ( exists($$self{has_column}{$col}) ) { push @out_cols, $col; }
        }
    }
    else 
    { 
        @out_cols = @{$$self{columns}}; 
    }
    $out .= "#". join("\t", @out_cols). "\n";
    
    return $out;
}


=head2 format_line

    About   : Returns the header.
    Usage   : $x = $vcf->next_data_hash(); print $vcf->format_line($x);
              $x = $vcf->next_data_array(); print $vcf->format_line($x);
    Args 1  : The columns or hash in the format returned by next_data_hash or next_data_array.
         2  : The columns to include [optional]

=cut

sub format_line
{
    my ($self,$record,$columns) = @_;

    if ( ref($record) eq 'HASH' ) { return $self->_format_line_hash($record,$columns); }
    if ( ref($record) eq 'ARRAY' ) { return join("\t",@$record)."\n"; }
    $self->throw("FIXME: todo .. " .ref($record). "\n");
}


=head2 recalc_ac_an

    About   : Control if the AC and AN values should be updated.
    Usage   : $vcf->recalc_ac_an(1); $x = $vcf->next_data_hash(); print $vcf->format_line($x);
    Args 1  : 0 .. never recalculate
              1 .. recalculate if present
              2 .. recalculate if present and add if missing

=cut

sub recalc_ac_an
{
    my ($self,$value) = @_;
    if ( $value eq '0' || $value eq '1' || $value eq '2' ) { $$self{recalc_ac_an} = $value; }
    return;
}

=head2 get_tag_index

    Usage   : my $idx = $vcf->get_tag_index('GT:PL:DP:SP:GQ','PL',':');
    Arg 1   : Field
        2   : The tag to find
        3   : Tag separator
    Returns : Index of the tag or -1 when not found

=cut

sub get_tag_index
{
    my ($self,$field,$tag,$sep) = @_;
    if ( !defined $field ) { return -1; }
    my $idx = 0;
    my $prev_isep = 0;
    my $isep = 0;
    while (1)
    {
        $isep = index($field,':',$prev_isep);
        if ( $isep==-1 ) 
        {
            if ( substr($field,$prev_isep) eq $tag ) { return $idx; }
            else { return -1; }
        }
        if ( substr($field,$prev_isep,$isep-$prev_isep) eq $tag ) { return $idx; }
        $prev_isep = $isep+1;
        $idx++;
    }
}

=head2 remove_field

    Usage   : my $field = $vcf->remove_field('GT:PL:DP:SP:GQ',1,':');    # returns 'GT:DP:SP:GQ'
    Arg 1   : Field
        2   : The index of the field to remove
        3   : Field separator
    Returns : Modified string

=cut

sub remove_field
{
    my ($self,$string,$idx,$sep) = @_;
    my $isep = -1;
    my $prev_isep = 0;
    my $itag = 0;
    while ($itag!=$idx)
    {
        $isep = index($string,$sep,$prev_isep);
        # The index may be out of range, VCFv4.1 allows omitting empty fields
        if ( $isep==-1 ) { return $string; }
        $prev_isep = $isep+1;
        $itag++;
    }
    my $out;
    if ( $isep>=0 ) { $out = substr($string,0,$isep); }
    my $ito=index($string,$sep,$isep+1);
    if ( $ito!=-1 ) 
    { 
        if ( $isep>=0 ) { $out .= ':' }
        $out .= substr($string,$ito+1); 
    }
    if ( !defined $out ) { return '.'; }
    return $out;
}

=head2 replace_field

    Usage   : my $col = $vcf->replace_field('GT:PL:DP:SP:GQ','XX',1,':');    # returns 'GT:XX:DP:SP:GQ'
    Arg 1   : Field
        2   : Replacement
        3   : 0-based index of the field to replace
        4   : Field separator
    Returns : Modified string

=cut

sub replace_field
{
    my ($self,$string,$repl,$idx,$sep) = @_;
    my $isep = -1;
    my $prev_isep = 0;
    my $itag = 0;
    while ($itag!=$idx)
    {
        $isep = index($string,$sep,$prev_isep);
        if ( $isep==-1 ) 
        { 
            # the out of range index may be OK, VCFv4.1 allows omitting empty fields
            if ( $$self{version}<4.1 ) 
            { 
                $self->throw("The index out of range ($string,$repl,$idx,$sep), missing fields not supported in VCFv$$self{version}."); 
            }
            while ( $itag<$idx ) { $string .= ':'; $itag++; }
            $string .= $repl;
            return $string;
        }
        $prev_isep = $isep+1;
        $itag++;
    }
    my $out;
    if ( $isep>=0 ) { $out = substr($string,0,$isep+1); }
    my $ito = index($string,$sep,$isep+1);
    if ( $ito==-1 )
    {
        $out .= $repl;
    }
    else
    { 
        $out .= $repl;
        $out .= ':';
        $out .= substr($string,$ito+1); 
    }
    if ( !defined $out ) { return '.'; }
    return $out;
}

=head2 get_info_field

    Usage   : my $line  = $vcf->next_line;
              my @items = split(/\t/,$line); 
              $af = $vcf->get_info_field('DP=14;AF=0.5;DB','AF');    # returns 0.5
              $af = $vcf->get_info_field('DP=14;AF=0.5;DB','DB');    # returns 1
              $af = $vcf->get_info_field('DP=14;AF=0.5;DB','XY');    # returns undef
    Arg 1   : The VCF line broken into an array
        2   : The tag to retrieve
    Returns : undef when tag is not present, the tag value if present, or 1 if flag is present

=cut

sub get_info_field
{
    my ($self,$info,$tag) = @_;

    my $info_len = length($info);
    my $tag_len = length($tag);
    my $idx = 0;
    while (1)
    {
        $idx = index($info,$tag,$idx);
        if ( $idx==-1 ) { return undef; }
        if ( $idx!=0 && substr($info,$idx-1,1) ne ';' ) { $idx += $tag_len; next; }
        if ( $tag_len+$idx >= $info_len ) { return 1; }

        my $follows = substr($info,$idx+$tag_len,1);
        if ( $follows eq ';' ) { return 1; }

        $idx += $tag_len;
        if ( $follows ne '=' ) { next; }

        $idx++;
        my $to = index($info,';',$idx);
        return $to==-1 ? substr($info,$idx) : substr($info,$idx,$to-$idx);
    }
}

=head2 get_field

    Usage   : my $line  = $vcf->next_line;
              my @items = split(/\t/,$line); 
              my $idx = $vcf->get_tag_index($$line[8],'PL',':'); 
              my $pl  = $vcf->get_field($$line[9],$idx) unless $idx==-1;
    Arg 1   : The VCF line broken into an array
        2   : The index of the field to retrieve
        3   : The delimiter [Default is ':']
    Returns : The tag value

=cut

sub get_field
{
    my ($self,$col,$idx,$delim) = @_;

    if ( !defined $delim ) { $delim=':'; }
    my $isep = 0;
    my $prev_isep = 0;
    my $itag = 0;
    while (1)
    {
        $isep = index($col,$delim,$prev_isep);
        if ( $itag==$idx ) { last; }
        if ( $isep==-1 ) { return '.'; }    # This is valid, missing fields can be ommited from genotype columns
        $prev_isep = $isep+1;
        $itag++;
    }
    return $isep<0 ? substr($col,$prev_isep) : substr($col,$prev_isep,$isep-$prev_isep);
}

=head2 get_sample_field

    Usage   : my $line  = $vcf->next_line;
              my @items = split(/\t/,$line); 
              my $idx = $vcf->get_tag_index($$line[8],'PL',':'); 
              my $pls = $vcf->get_sample_field(\@items,$idx) unless $idx==-1;
    Arg 1   : The VCF line broken into an array
        2   : The index of the field to retrieve
    Returns : Array of values

=cut

sub get_sample_field
{
    my ($self,$cols,$idx) = @_;
    my @out;
    my $n = @$cols;
    for (my $icol=9; $icol<$n; $icol++)
    {
        my $col = $$cols[$icol];
        my $isep = 0;
        my $prev_isep = 0;
        my $itag = 0;
        while (1)
        {
            $isep = index($col,':',$prev_isep);
            if ( $itag==$idx ) { last; }
            if ( $isep==-1 ) { return '.'; }    # This is valid, missing fields can be ommited from genotype columns
            $prev_isep = $isep+1;
            $itag++;
        }
        my $val = $isep<0 ? substr($col,$prev_isep) : substr($col,$prev_isep,$isep-$prev_isep);
        push @out,$val;
    }
    return \@out;
}


=head2 split_mandatory

    About   : Faster alternative to regexs, extract the mandatory columns
    Usage   : my $line=$vcf->next_line; my @cols = $vcf->split_mandatory($line);
    Arg     : 
    Returns : Pointer to the array of values

=cut

sub split_mandatory
{
    my ($self,$line) = @_;
    my @out;
    my $prev = 0;
    for (my $i=0; $i<7; $i++)
    {
        my $isep = index($line,"\t",$prev);
        if ( $isep==-1 ) { $self->throw("Could not parse the mandatory columns: $line"); }
        push @out, substr($line,$prev,$isep-$prev);
        $prev = $isep+1;
    }
    my $isep = index($line,"\t",$prev);
    if ( $isep!=-1 )
    {
        push @out, substr($line,$prev,$isep-$prev-1);
    }
    else
    {
        push @out, substr($line,$prev);
    }
    return \@out;
}


=head2 split_gt

    About   : Faster alternative to regexs
    Usage   : my ($a1,$a2,$a3) = $vcf->split_gt('0/0/1'); # returns (0,0,1)
    Arg     : Diploid genotype to split into alleles
    Returns : Array of values

=cut

sub split_gt
{
    my ($self,$gt) = @_;
    my @als;
    my $iprev = 0;
    while (1)
    {
        my $isep = index($gt,'/',$iprev);
        my $jsep = index($gt,'|',$iprev);
        if ( $isep<0 or ($jsep>=0 && $jsep<$isep) ) { $isep = $jsep; }
        push @als, $isep<0 ? substr($gt,$iprev) : substr($gt,$iprev,$isep-$iprev);
        if ( $isep<0 ) { return (@als); }
        $iprev = $isep+1;
    }
    return (@als);
}

=head2 split_by

    About   : Generalization of split_gt
    Usage   : my ($a1,$a2,$a3) = $vcf->split_gt('0/0|1',qw(| /)); # returns (0,0,1)
    Arg     : Diploid genotype to split into alleles
    Returns : Array of values

=cut

sub split_by
{
    my ($self,$str,@seps) = @_;
    my @out;
    my $iprev = 0;
    while (1)
    {
        my $min;
        for my $sep (@seps)
        {
            my $idx = index($str,$sep,$iprev);
            if ( $idx==-1 ) { next; }
            if ( !defined $min or $idx<$min ) { $min=$idx }
        }
        push @out, defined $min ? substr($str,$iprev,$min-$iprev) : substr($str,$iprev);
        if ( !defined $min ) { return @out; }
        $iprev = $min+1;
    }
    return (@out);
}

=head2 decode_genotype

    About   : Faster alternative to regexs
    Usage   : my $gt = $vcf->decode_genotype('G',['A','C'],'0/0'); # returns 'G/G'
    Arg   1 : Ref allele
          2 : Alt alleles
          3 : The genotype to decode
    Returns : Decoded GT string

=cut

sub decode_genotype
{
    my ($self,$ref,$alt,$gt) = @_;
    my $isep = 0;
    my $out;
    while (1)
    {
        my $i = index($gt,'/',$isep);
        my $j = index($gt,'|',$isep);
        if ( $i==-1 && $j==-1 )
        {
            my $idx = substr($gt,$isep);
            if ( $idx eq '.' )
            {
                $out .= $idx;
            }
            else
            {
                if ( $idx>@$alt ) { $self->throw("The genotype index $idx in $gt is out of bounds: ", join(',',@$alt)); }
                $out .= $idx==0 ? $ref : $$alt[$idx-1];
            }
            return $out;
        }
        if ( $i!=-1 && $j!=-1 && $i>$j ) { $i=$j; }
        elsif ( $i==-1 ) { $i=$j }

        my $idx = substr($gt,$isep,$i-$isep);
        if ( $idx eq '.' )
        {
            $out .= $idx;
        }
        else
        {
            if ( $idx>@$alt ) { $self->throw("The genotype index $idx in $gt out of bounds: ", join(',',@$alt)); }
            $out .= $idx==0 ? $ref : $$alt[$idx-1];
        }
        $out .= substr($gt,$i,1);
        $isep = $i+1;
    }
}


sub _format_line_hash
{
    my ($self,$record,$columns) = @_;

    if ( !$$self{columns} ) 
    { 
        my $ngtypes = scalar keys %{$$record{gtypes}};
        if ( !$ngtypes && !exists($$record{FORMAT}) ) { $ngtypes--; }
        $self->_fake_column_names($ngtypes); 
    }
    my $cols = $$self{columns};

    # CHROM  POS     ID      REF
    my $out;
    $out .= $$record{CHROM} . "\t";
    $out .= $$record{POS} . "\t";
    $out .= (defined $$record{ID} ? $$record{ID} : '.') . "\t";
    $out .= $$record{REF} . "\t";

    # ALT
    $out .= join(',',@{$$record{ALT}} ? @{$$record{ALT}} : '.');

    # QUAL
    $out .= "\t". $$record{QUAL};

    # FILTER
    $out .= "\t". join(';',$$record{FILTER} ? @{$$record{FILTER}} : '.');

    # Collect the gtypes of interest
    my $gtypes;
    if ( $columns )
    {
        # Select only those gtypes keys with a corresponding key in columns.
        for my $col (@$columns) { $$gtypes{$col} = $$record{gtypes}{$col}; }
    }
    else
    {
        $gtypes = $$record{gtypes};
    }

    # INFO
    # .. calculate NS, AN and AC, but only if recalc_ac_an is set
    my $needs_an_ac = $$self{recalc_ac_an}==2 ? 1 : 0;
    my @info;
    while (my ($key,$value) = each %{$$record{INFO}})
    {
        if ( $$self{recalc_ac_an}>0 )
        {
            if ( $key eq 'AN' ) { $needs_an_ac=1; next; }
            if ( $key eq 'AC' ) { $needs_an_ac=1; next; }
        }
        if ( defined $value )
        {
            push @info, "$key=$value";
        }
        elsif ( $key ne '.' )
        {
            push @info, $key;
        }
    }
    if ( $needs_an_ac )
    {
        my $nalt = scalar @{$$record{ALT}};
        if ( $nalt==1 && $$record{ALT}[0] eq '.' ) { $nalt=0; } 
        my ($an,$ac) = $self->calc_an_ac($gtypes,$nalt);
        push @info, "AN=$an","AC=$ac";
    }
    if ( !@info ) { push @info, '.'; }
    $out .= "\t". join(';', sort @info);

    # FORMAT, the column is not required, it may not be present when there are no genotypes
    if ( exists($$cols[8]) && defined $$record{FORMAT} )
    {
        $out .= "\t". join(':',@{$$record{FORMAT}});
    }

    # Genotypes: output all columns or only a selection?
    my @col_names = $columns ? @$columns : @$cols[9..@$cols-1];
    my $nformat = defined $$record{FORMAT} ? @{$$record{FORMAT}} : 0;
    for my $col (@col_names)
    {
        my $gt = $$gtypes{$col};
        my $can_drop = $$self{drop_trailings};
        my @gtype;
        for (my $i=$nformat-1; $i>=0; $i--)
        {
            my $field = $$record{FORMAT}[$i];
            if ( $i==0 ) { $can_drop=0; }

            if ( exists($$gt{$field}) ) 
            {
                $can_drop = 0;
                if ( ref($$gt{$field}) eq 'HASH' ) 
                {
                    # Special treatment for Number=[AG] tags
                    unshift @gtype, $self->format_AGtag($record,$col,$$gt{$field},$field);
                }
                else
                { 
                    unshift @gtype,$$gt{$field}; 
                }
            }
            elsif ( $can_drop ) { next; }
            elsif ( exists($$self{header}{FORMAT}{$field}{default}) ) { unshift @gtype,$$self{header}{FORMAT}{$field}{default}; $can_drop=0; }
            else { $self->throw(qq[No value for the field "$field" and no default available, column "$col" at $$record{CHROM}:$$record{POS}.\n]); }
        }
        $out .= "\t" . join(':',@gtype);
    }

    $out .= "\n";
    return $out;
}

sub calc_an_ac
{
    my ($self,$gtypes,$nalleles) = @_;
    my $sep_re = $$self{regex_gtsep};
    my ($an,%ac_counts);
    if ( defined $nalleles )
    {
        for (my $i=1; $i<=$nalleles; $i++) { $ac_counts{$i}=0; }
    }
    $an = 0;
    for my $gt (keys %$gtypes)
    {
        my $value = $$gtypes{$gt}{GT};
        if ( !defined $value ) { next; } # GT may not be present
        my ($al1,$al2) = split($sep_re,$value);
        if ( defined($al1) && $al1 ne '.' )
        {
            $an++;
            if ( $al1 ne '0' ) { $ac_counts{$al1}++; }
        }
        if ( defined($al2) && $al2 ne '.' )
        {
            $an++;
            if ( $al2 ne '0' ) { $ac_counts{$al2}++; }
        }
    }
    my @ac;
    for my $ac ( sort { $a <=> $b } keys %ac_counts) { push @ac, $ac_counts{$ac}; }
    if ( !@ac ) { @ac = ('0'); }
    return ($an,join(',',@ac),\@ac);
}

sub _validate_alt_field
{
    my ($self,$values,$ref) = @_;

    for (my $i=0; $i<@$values; $i++)
    {
        for (my $j=0; $j<$i; $j++)
        {
            if ( $$values[$i] eq $$values[$j] ) { return "The alleles not unique: $$values[$i]"; }
        }
        if ( $$values[$i] eq $ref ) { return "REF allele listed in the ALT field??"; }
    }
    return undef;
}

=head2 validate_alt_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_alt_field($$x{ALT});
    Args    : The ALT arrayref
    Returns : Error message in case of an error.

=cut

sub validate_alt_field
{
    my ($self,$values,$ref) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }

    my $ret = $self->_validate_alt_field($values,$ref);
    if ( $ret ) { return $ret; }
    
    my @err;
    for my $item (@$values)
    {
        if ( $item=~/^[ACTGN]$/ ) { next; }
        elsif ( $item=~/^I[ACTGN]+$/ ) { next; }
        elsif ( $item=~/^D\d+$/ ) { next; }

        push @err, $item;
    }
    if ( !@err ) { return undef; }
    return 'Could not parse the allele(s) [' .join(',',@err). ']';
}

=head2 event_type

    Usage   :   my $x = $vcf->next_data_hash(); 
                my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x,'NA00001');
                for my $allele (@$alleles)
                {
                    my ($type,$len,$ht) = $vcf->event_type($x,$allele);
                }
              or
                my ($type,$len,$ht) = $vcf->event_type($ref,$al);
    Args    : VCF data line parsed by next_data_hash or the reference allele
            : Allele
    Returns :   's' for SNP and number of SNPs in the record
                'i' for indel and a positive (resp. negative) number for the length of insertion (resp. deletion)
                'r' identical to the reference, length 0
                'o' for other (complex events) and the number of affected bases
                'b' breakend
                'u' unknown

=cut

sub event_type
{
    my ($self,$rec,$allele) = @_;

    my $ref = $rec;
    if ( ref($rec) eq 'HASH' ) 
    { 
        if ( exists($$rec{_cached_events}{$allele}) ) { return (@{$$rec{_cached_events}{$allele}}); }
        $ref = $$rec{REF};
    }

    my ($type,$len,$ht);
    if ( $allele eq $ref or $allele eq '.' ) { $len=0; $type='r'; $ht=$ref; }
    elsif ( $allele=~/^[ACGT]$/ ) { $len=1; $type='s'; $ht=$allele; }
    elsif ( $allele=~/^I/ ) { $len=length($allele)-1; $type='i'; $ht=$'; }
    elsif ( $allele=~/^D(\d+)/ ) { $len=-$1; $type='i'; $ht=''; }
    else 
    { 
        my $chr = ref($rec) eq 'HASH' ? $$rec{CHROM} : 'undef';
        my $pos = ref($rec) eq 'HASH' ? $$rec{POS} : 'undef';
        $self->throw("Eh?: $chr:$pos .. $ref $allele\n");
    }

    if ( ref($rec) eq 'HASH' ) 
    {
        $$rec{_cached_events}{$allele} = [$type,$len,$ht];
    }
    return ($type,$len,$ht);
}

=head2 has_AGtags

    About   : Checks the header for the presence of tags with variable number of fields (Number=A or Number=G, such as GL)
    Usage   : $vcf->parse_header(); my $agtags = $vcf->has_AGtags();
    Args    : None
    Returns : Hash {fmtA=>[tags],fmtG=>[tags],infoA=>[tags],infoG=>[tags]} or undef if none is present

=cut

sub has_AGtags
{
    my ($self) = @_;
    my $out;
    if ( exists($$self{header}{FORMAT}) )
    {
        for my $tag (keys %{$$self{header}{FORMAT}})
        {
            if ( $$self{header}{FORMAT}{$tag}{Number} eq 'A' ) { push @{$$out{fmtA}},$tag; }
            if ( $$self{header}{FORMAT}{$tag}{Number} eq 'G' ) { push @{$$out{fmtG}},$tag; }
        }
    }
    if ( exists($$self{header}{INFO}) )
    {
        for my $tag (keys %{$$self{header}{INFO}})
        {
            if ( $$self{header}{INFO}{$tag}{Number} eq 'A' ) { push @{$$out{infoA}},$tag; }
            if ( $$self{header}{INFO}{$tag}{Number} eq 'G' ) { push @{$$out{infoG}},$tag; }
        }
    }
    if ( defined $out ) 
    {
        for my $key (qw(fmtA fmtG infoA infoG)) { if ( !exists($$out{$key}) ) { $$out{$key}=[] } }
    }
    return $out;
}

=head2 parse_AGtags

    About   : Breaks tags with variable number of fields (that is where Number is set to 'A' or 'G', such as GL) into hashes
    Usage   : my $x = $vcf->next_data_hash(); my $values = $vcf->parse_AGtags($x);
    Args    : VCF data line parsed by next_data_hash
            : Mapping between ALT representations based on different REFs [optional]
            : New REF [optional]
    Returns : Hash {Allele=>Value}

=cut

sub parse_AGtags
{
    my ($self,$rec,$ref_alt_map,$new_ref) = @_;

    if ( !exists($$rec{gtypes}) ) { return; }

    my (@atags,@gtags);
    for my $fmt (@{$$rec{FORMAT}})
    {
        # These have been listed explicitly for proper merging of v4.0  VCFs
        if ( $$self{fix_v40_AGtags} )
        {
            if ( $fmt eq 'GL' or $fmt eq 'PL' ) { push @gtags,$fmt; next; }
            if ( $fmt eq 'AC' or $fmt eq 'AF'  ) { push @atags,$fmt; next; }
        }
        if ( !exists($$self{header}{FORMAT}{$fmt}) ) { next; }
        if ( $$self{header}{FORMAT}{$fmt}{Number} eq 'A' ) { push @atags,$fmt; next; }
        if ( $$self{header}{FORMAT}{$fmt}{Number} eq 'G' ) { push @gtags,$fmt; next; }
    }
    my $missing = $$self{defaults}{default};
    if ( @atags )
    {
        # Parse Number=A tags
        my $alts;
        if ( defined $ref_alt_map )
        {
            $alts = [];
            for my $alt (@{$$rec{ALT}})
            {
                if ( !exists($$ref_alt_map{$new_ref}{$alt}) ) { $self->throw("FIXME: $new_ref $alt...?\n"); }
                push @$alts, $$ref_alt_map{$new_ref}{$alt};
            }
        }
        else
        {
            $alts = $$rec{ALT};
        }
        for my $tag (@atags)
        {
            for my $sample (values %{$$rec{gtypes}})
            {
                if ( !exists($$sample{$tag}) or $$sample{$tag} eq $missing ) { next; }
                my @values = split(/,/,$$sample{$tag});
                $$sample{$tag} = {};
                for (my $i=0; $i<@values; $i++)
                {
                    $$sample{$tag}{$$alts[$i]} = $values[$i];
                }
            }
        }
    }
    if ( @gtags )
    {
        # Parse Number=G tags
        my @alleles;
        if ( defined $ref_alt_map ) 
        {
            push @alleles, $new_ref;
            for my $alt (@{$$rec{ALT}})
            {
                if ( !exists($$ref_alt_map{$new_ref}{$alt}) ) { $self->throw("FIXME: [$new_ref] [$alt]...?\n", Dumper($ref_alt_map,$rec)); }
                push @alleles, $$ref_alt_map{$new_ref}{$alt};
            }
        }
        else
        {
            @alleles = ($$rec{REF},@{$$rec{ALT}});
            if ( @alleles==2 && $alleles[1] eq '.' ) { pop(@alleles); }
        }
        my @gtypes;
        for (my $i=0; $i<@alleles; $i++)
        {
            for (my $j=0; $j<=$i; $j++)
            {
                push @{$gtypes[1]}, $alleles[$i].'/'.$alleles[$j];
            }
            push @{$gtypes[0]}, $alleles[$i];
        }
        for my $tag (@gtags)
        {
            for my $name (keys %{$$rec{gtypes}})
            {
                my $sample = $$rec{gtypes}{$name};
                if ( !exists($$sample{$tag}) or $$sample{$tag} eq $missing ) { next; }
                my @values = split(/,/,$$sample{$tag});
                my $ploidy = $self->guess_ploidy(scalar @alleles, scalar @values) - 1;
                if ( $ploidy<0 ) 
                { 
                    my $nals  = scalar @alleles;
                    my $nvals = scalar @values;
                    my $ndip  = $nals*($nals+1)/2;
                    $self->throw(
                        "Wrong number of values in $name/$tag at $$rec{CHROM}:$$rec{POS} .. nAlleles=$nals, nValues=$nvals.\n".
                        "Expected $ndip values for diploid genotypes or $nals for haploid genotypes.\n");
                }
                if ( $ploidy>1 ) { $self->throw("Sorry, not ready for ploidy bigger than 2\n"); }
                if ( $ploidy!=1 ) { $$rec{_cached_ploidy}{$name} = $ploidy; }
                $$sample{$tag} = {};
                for (my $i=0; $i<@values; $i++)
                {
                    $$sample{$tag}{$gtypes[$ploidy][$i]} = $values[$i];
                }
            }
        }
    }
}

=head2 format_AGtag

    About   : Format tag with variable number of fields (that is where Number is set to 'A' or 'G', such as GL)
    Usage   : 
    Args    : 
            : 
            : 
    Returns : 

=cut

sub format_AGtag
{
    my ($self,$record,$sample,$tag_data,$tag) = @_;

    # The FORMAT field is checked only once and the results are cached.
    if ( !exists($$record{_atags}) )
    {
        $$record{_atags} = {};

        # Check if there are any A,G tags
        for my $fmt (@{$$record{FORMAT}})
        {
            # These have been listed explicitly for proper merging of v4.0  VCFs
            if ( $$self{fix_v40_AGtags} )
            {
                if ( $fmt eq 'GL' or $fmt eq 'PL' ) { $$record{_gtags}{$fmt}=1; next; }
                if ( $fmt eq 'AC' or $fmt eq 'AF'  ) { $$record{_atags}{$fmt}=1; next; }
            }
            if ( !exists($$self{header}{FORMAT}{$fmt}) ) { next; }
            if ( $$self{header}{FORMAT}{$fmt}{Number} eq 'A' ) { $$record{_atags}{$fmt}=1; next; }
            if ( $$self{header}{FORMAT}{$fmt}{Number} eq 'G' ) { $$record{_gtags}{$fmt}=1; next; }
        }
    }

    my @out;
    if ( exists($$record{_atags}{$tag}) )
    {
        for my $alt (@{$$record{ALT}})
        {
            push @out, exists($$tag_data{$alt}) ? $$tag_data{$alt} : $$self{defaults}{default};    
        }
    }

    if ( exists($$record{_gtags}{$tag}) )
    {
        my $gtypes  = $$record{_gtypes};
        my $gtypes2 = $$record{_gtypes2};
        if ( !defined $gtypes )
        {
            $gtypes  = [];
            $gtypes2 = [];

            my @alleles = ( $$record{REF}, @{$$record{ALT}} );
            for (my $i=0; $i<@alleles; $i++)
            {
                for (my $j=0; $j<=$i; $j++)
                {
                    push @{$$gtypes[1]}, $alleles[$i].'/'.$alleles[$j];
                    push @{$$gtypes2[1]}, $alleles[$j].'/'.$alleles[$i];
                }
                push @{$$gtypes[0]}, $alleles[$i];
            }
            
            $$record{_gtypes}  = $gtypes;
            $$record{_gtypes2} = $gtypes2;
        }

        my $ploidy = exists($$record{_cached_ploidy}{$sample}) ? $$record{_cached_ploidy}{$sample} : 1;
        for (my $i=0; $i<@{$$gtypes[$ploidy]}; $i++)
        {
            my $gt = $$gtypes[$ploidy][$i];
            if ( !exists($$tag_data{$gt}) ) { $gt = $$gtypes2[$ploidy][$i]; }
            push @out, exists($$tag_data{$gt}) ? $$tag_data{$gt} : $$self{defaults}{default};
        }
    }

    return join(',',@out);
}

=head2 parse_alleles

    About   : Deprecated, use parse_haplotype instead.
    Usage   : my $x = $vcf->next_data_hash(); my ($al1,$sep,$al2) = $vcf->parse_alleles($x,'NA00001');
    Args    : VCF data line parsed by next_data_hash
            : The genotype column name
    Returns : Alleles and the separator. If only one allele is present, $sep and $al2 will be an empty string.

=cut

sub parse_alleles
{
    my ($self,$rec,$column) = @_;
    if ( !exists($$rec{gtypes}) || !exists($$rec{gtypes}{$column}) ) { $self->throw("The column not present: '$column'\n"); }

    my $gtype = $$rec{gtypes}{$column}{GT};
    if ( !($gtype=~$$self{regex_gt}) ) { $self->throw("Could not parse gtype string [$gtype] [$$rec{CHROM}:$$rec{POS}]\n"); }
    my $al1 = $1;
    my $sep = $2;
    my $al2 = $3;

    if ( !$al1 ) { $al1 = $$rec{REF}; }
    elsif ( $al1 ne '.' ) 
    { 
        if ( !($al1=~/^\d+$/) ) { $self->throw("Uh, what is this? [$al1] $$rec{CHROM}:$$rec{POS}\n"); } 
        $al1 = $$rec{ALT}[$al1-1]; 
    }

    if ( !defined $al2 or $al2 eq '' )
    {
        $sep = '';
        $al2 = '';
    }
    else
    {
        if ( !$al2 ) { $al2 = $$rec{REF}; }
        elsif ( $al2 ne '.' ) { $al2 = $$rec{ALT}[$al2-1]; }
    }
    return ($al1,$sep,$al2);
}

=head2 parse_haplotype

    About   : Similar to parse_alleles, supports also multiploid VCFs. 
    Usage   : my $x = $vcf->next_data_hash(); my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x,'NA00001');
    Args    : VCF data line parsed by next_data_hash
            : The genotype column name
    Returns : Two array refs and two boolean flags: List of alleles, list of separators, and is_phased/empty flags. The values
                can be cashed and must be therefore considered read only!

=cut

sub parse_haplotype
{
    my ($self,$rec,$column) = @_;
    if ( !exists($$rec{gtypes}{$column}) ) { $self->throw("The column not present: '$column'\n"); }
    if ( !exists($$rec{gtypes}{$column}{GT}) ) { return (['.'],[],0,1); }

    my $gtype = $$rec{gtypes}{$column}{GT};
    if ( exists($$rec{_cached_haplotypes}{$gtype}) ) { return (@{$$rec{_cached_haplotypes}{$gtype}}); }

    my @alleles   = ();
    my @seps      = ();
    my $is_phased = 0;
    my $is_empty  = 1;

    my $buf = $gtype;
    while ($buf ne '')
    {
        if ( !($buf=~m{^(\.|\d+)([|/]?)}) ) { $self->throw("Could not parse gtype string [$gtype] .. $$rec{CHROM}:$$rec{POS} $column\n"); }
        $buf = $';

        if ( $1 eq '.' ) { push @alleles,'.'; }
        else
        {
            $is_empty = 0;
            if ( $1 eq '0' ) { push @alleles,$$rec{REF}; }
            elsif ( exists($$rec{ALT}[$1-1]) ) { push @alleles,$$rec{ALT}[$1-1]; }
            else { $self->throw(qq[The haplotype indexes in "$gtype" do not match the ALT column .. $$rec{CHROM}:$$rec{POS} $column\n]); }
        }
        if ( $2 )
        {
            if ( $2 eq '|' ) { $is_phased=1; }
            push @seps,$2;
        }
    }
    $$rec{_cached_haplotypes}{$gtype} = [\@alleles,\@seps,$is_phased,$is_empty];
    return (@{$$rec{_cached_haplotypes}{$gtype}});
}

=head2 format_haplotype

    Usage   : my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x,'NA00001'); print $vcf->format_haplotype($alleles,$seps);

=cut

sub format_haplotype
{
    my ($self,$alleles,$seps) = @_;
    if ( @$alleles != @$seps+1 ) { $self->throw(sprintf("Uh: %d vs %d\n",scalar @$alleles,scalar @$seps),Dumper($alleles,$seps)); }
    my $out = $$alleles[0];
    for (my $i=1; $i<@$alleles; $i++)
    {
        $out .= $$seps[$i-1];
        $out .= $$alleles[$i];
    }
    return $out;
}


=head2 format_genotype_strings

    Usage   : my $x = { REF=>'A', gtypes=>{'NA00001'=>{'GT'=>'A/C'}}, FORMAT=>['GT'], CHROM=>1, POS=>1, FILTER=>['.'], QUAL=>-1 };
              $vcf->format_genotype_strings($x); 
              print $vcf->format_line($x);
    Args 1  : VCF data line in the format as if parsed by next_data_hash with alleles written as letters.
         2  : Optionally, a subset of columns can be supplied. See also format_line.
    Returns : Modifies the ALT array and the genotypes so that ref alleles become 0 and non-ref alleles 
                numbers starting from 1. If the key $$vcf{trim_redundant_ALTs} is set, ALT alleles not appearing
                in any of the sample column will be removed.

=cut

sub format_genotype_strings
{
    my ($self,$rec,$columns) = @_;

    if ( !exists($$rec{gtypes}) ) { return; }

    my $ref = $$rec{REF};
    my $nalts = 0;
    my %alts  = ();

    if ( !$columns ) { $columns = [keys %{$$rec{gtypes}}]; }

    for my $key (@$columns)
    {
        my $gtype = $$rec{gtypes}{$key}{GT};
        my $buf = $gtype;
        my $out = '';
        while ($buf ne '')
        {
            $buf=~m{^([^/|]+)([/|]?)};
            $buf = $';

            my $al  = $1;
            my $sep = $2;
            if ( $al eq $ref or $al eq '0' or $al eq '*' ) { $al=0; }
            else
            {
                if ( $al=~/^\d+$/ ) 
                { 
                    if ( !exists($$rec{ALT}[$al-1]) ) { $self->throw("Broken ALT, index $al out of bounds\n"); }
                    $al = $$rec{ALT}[$al-1]; 
                }

                if ( exists($alts{$al}) ) { $al = $alts{$al} }
                elsif ( $al ne '.' )
                {
                    $alts{$al} = ++$nalts;
                    $al = $nalts;
                }

            }
            $out .= $al;
            if ( $sep ) { $out .= $sep; }
        }
        $$rec{gtypes}{$key}{GT} = $out;
    }
    if ( !$$self{trim_redundant_ALTs} && exists($$rec{ALT}) && @{$$rec{ALT}} )
    {
        for my $alt (@{$$rec{ALT}}) 
        { 
            if ( !exists($alts{$alt}) ) { $alts{$alt} = ++$nalts;  }
        }
    }
    $$rec{ALT} = [ sort { $alts{$a}<=>$alts{$b} } keys %alts ];
}

sub fill_ref_alt_mapping
{
    my ($self,$map) = @_;
    
    my $new_ref;
    for my $ref (keys %$map)
    {
        $new_ref = $ref;
        if ( $ref ne $new_ref ) { $self->warn("The reference prefixes do not agree: $ref vs $new_ref\n"); return undef; }
        for my $alt (keys %{$$map{$ref}})
        {
            $$map{$ref}{$alt} = $alt;
        }
    }
    $$map{$new_ref}{$new_ref} = $new_ref;
    return $new_ref;
}


=head2 format_header_line

    Usage   : $vcf->format_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'})
    Args    : 
    Returns : 

=cut

sub format_header_line
{
    my ($self,$rec) = @_;
    my $line = "##$$rec{key}";
    $line .= "=$$rec{value}" unless !exists($$rec{value});
    $line .= "=$$rec{ID}" unless !exists($$rec{ID});
    $line .= ",$$rec{Number}" unless !exists($$rec{Number});
    $line .= ",$$rec{Type}" unless !exists($$rec{Type});
    $line .= qq[,"$$rec{Description}"] unless !exists($$rec{Description});
    $line .= "\n";
    return $line;
}

=head2 remove_columns

    Usage   : my $rec=$vcf->next_data_hash(); $vcf->remove_columns($rec,remove=>['NA001','NA0002']);
    Args    : VCF hash pointer
            : list of columns to remove or a lookup hash with column names to keep (remove=>[] or keep=>{})
    Returns : 

=cut

sub remove_columns
{
    my ($self,$rec,%args) = @_;
    if ( ref($rec) ne 'HASH' ) { $self->throw("TODO: rec for array"); }
    if ( exists($args{keep}) )
    {
        for my $col (keys %{$$rec{gtypes}})
        {
            if ( !exists($args{keep}{$col}) ) { delete($$rec{gtypes}{$col}); }
        }
    }
    if ( exists($args{remove}) )
    {
        for my $col (@{$args{remove}})
        {
            if ( exists($$rec{gtypes}{$col}) ) { delete($$rec{gtypes}{$col}); }
        }
    }
}

=head2 add_columns

    Usage   : $vcf->add_columns('NA001','NA0002');
    Args    : 
    Returns : 

=cut

sub add_columns
{
    my ($self,@columns) = @_;
    if ( !$$self{columns} ) 
    { 
        # The columns should be initialized de novo. Figure out if the @columns contain also the mandatory
        #   columns and if FORMAT should be present (it can be absent when there is no genotype column present).
        my $has_other = 0;
        for my $col (@columns)
        {
            if ( !exists($$self{reserved}{cols}{$col}) ) { $has_other=1; last; }
        }

        $$self{columns} = [ @{$$self{mandatory}} ]; 
        if ( $has_other ) { push @{$$self{columns}},'FORMAT'; }

        for my $col (@{$$self{columns}}) { $$self{has_column}{$col}=1; }
    }
    my $ncols = @{$$self{columns}};
    for my $col (@columns)
    {
        if ( $$self{has_column}{$col} ) { next; }
        $ncols++;
        push @{$$self{columns}}, $col;
    }
}

=head2 add_format_field

    Usage   : $x=$vcf->next_data_hash(); $vcf->add_format_field($x,'FOO'); $$x{gtypes}{NA0001}{FOO}='Bar'; print $vcf->format_line($x);
    Args    : The record obtained by next_data_hash
            : The field name
    Returns : 

=cut

sub add_format_field
{
    my ($self,$rec,$field) = @_;

    if ( !$$rec{FORMAT} ) { $$rec{FORMAT}=[]; }

    for my $key (@{$$rec{FORMAT}})
    {
        if ( $key eq $field ) { return; } # already there
    }
    push @{$$rec{FORMAT}}, $field;
}


=head2 remove_format_field

    Usage   : $x=$vcf->next_data_hash(); $vcf->remove_format_field($x,'FOO'); print $vcf->format_line($x);
    Args    : The record obtained by next_data_hash
            : The field name
    Returns : 

=cut

sub remove_format_field
{
    my ($self,$rec,$field) = @_;

    if ( !$$rec{FORMAT} ) { $$rec{FORMAT}=[]; }

    my $i = 0;
    for my $key (@{$$rec{FORMAT}})
    {
        if ( $key eq $field ) { splice @{$$rec{FORMAT}},$i,1; }
        $i++;
    }
}


=head2 add_info_field

    Usage   : $x=$vcf->next_data_array(); $$x[7]=$vcf->add_info_field($$x[7],'FOO'=>'value','BAR'=>undef,'BAZ'=>''); print join("\t",@$x)."\n";
    Args    : The record obtained by next_data_array
            : The INFO field name and value pairs. If value is undef and the key is present in $$x[7],
                it will be removed. To add fields without a value, use empty string ''.
    Returns : The formatted INFO.

=cut

sub add_info_field
{
    my ($self,$info,%fields) = @_;

    my @out = ();

    # First handle the existing values, keep everything unless in %fields
    for my $field (split(/;/,$info))
    {
        my ($key,$value) = split(/=/,$field);
        if ( $key eq '.' ) { next; }
        if ( !exists($fields{$key}) ) { push @out,$field; next; }
    }

    # Now add the new values and remove the unwanted ones
    while (my ($key,$value)=each %fields)
    {
        if ( !defined($value) ) { next; }       # this one should be removed
        if ( $value eq '' ) { push @out,$key; } # this one is of the form HM2 in contrast to DP=3
        else { push @out,"$key=$value"; }       # this is the standard key=value pair
    }
    if ( !@out ) { push @out,'.'; }
    return join(';',@out);
}


=head2 add_filter

    Usage   : $x=$vcf->next_data_array(); $$x[6]=$vcf->add_filter($$x[6],'SnpCluster'=>1,'q10'=>0); print join("\t",@$x)."\n";
    Args    : The record obtained by next_data_array or next_data_hash
            : The key-value pairs for filter to be added. If value is 1, the filter will be added. If 0, the filter will be removed.
    Returns : The formatted filter field.

=cut

sub add_filter
{
    my ($self,$filter,%filters) = @_;

    my @out = ();
    my @filters = ref($filter) eq 'ARRAY' ? @$filter : split(/;/,$filter);

    # First handle the existing filters, keep everything unless in %filters
    for my $key (@filters)
    {
        if ( $key eq '.' or $key eq 'PASS' ) { next; }
        if ( !exists($filters{$key}) ) { push @out,$key; next; }
    }

    # Now add the new filters and remove the unwanted ones
    while (my ($key,$value)=each %filters)
    {
        if ( !$value ) { next; }                # this one should be removed
        push @out,$key;                         # this one should be added
    }
    if ( !@out ) { push @out,'PASS'; }
    return ref($filter) eq 'ARRAY' ? return \@out : join(';',@out);
}


=head2 validate_filter_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_filter_field($$x{FILTER});
    Args    : The FILTER arrayref
    Returns : Error message in case of an error.

=cut

sub validate_filter_field
{
    my ($self,$values) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }
    
    my @errs;
    my @missing;
    for my $item (@$values)
    {
        if ( $item eq $$self{filter_passed} ) { next; }
        if ( $item=~/,/ ) { push @errs,"Expected semicolon as a separator."; }
        if ( exists($$self{reserved}{FILTER}{$item}) ) { return qq[The filter name "$item" cannot be used, it is a reserved word.]; }
        if ( exists($$self{header}{FILTER}{$item}) ) { next; }
        push @missing, $item;
        $self->add_header_line({key=>'FILTER',ID=>$item,Description=>'No description'});
    }
    if ( !@errs && !@missing ) { return undef; }
    if ( $$self{version}<3.3 ) { return undef; }
    return join(',',@errs) .' '. 'The filter(s) [' . join(',',@missing) . '] not listed in the header.';
}


sub _add_unknown_field
{
    my ($self,$field,$key,$nargs) = @_;
    $self->add_header_line({key=>$field,ID=>$key,Number=>$nargs,Type=>'String',Description=>'No description'});
}

=head2 validate_header

    About   : Version specific header validation code.
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header(); $vcf->validate_header();
    Args    :

=cut

sub validate_header
{
    my ($self) = @_;
}

=head2 validate_line

    About   : Version specific line validation code.
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header(); $x = $vcf->next_data_hash; $vcf->validate_line($x);
    Args    :

=cut

sub validate_line
{
    my ($self,$x) = @_;

    # Is the ID composed of alphanumeric chars
    if ( !($$x{ID}=~/^[\w;\.]+$/) ) { $self->warn("Expected alphanumeric ID at $$x{CHROM}:$$x{POS}, but got [$$x{ID}]\n"); }
}

=head2 validate_info_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_info_field($$x{INFO},$$x{ALT});
    Args    : The INFO hashref
    Returns : Error message in case of an error.

=cut

sub validate_info_field
{
    my ($self,$values,$alts) = @_;

    if ( !defined $values ) { return 'Empty INFO field.'; }

    # First handle the empty INFO field (.)
    if ( scalar keys %$values == 1 && exists($$values{'.'}) ) { return undef; }

    # Expected numbers
    my $ng = -1;
    my $na = -1;
    my $nr = -1;
    if ( $$self{version}>4.0 )
    {
        if ( $$alts[0] eq '.' ) { $ng=1; $na=1; }
        else
        {
            $na = @$alts;
            $ng = (1+$na+1)*($na+1)/2;
            $nr = $na+1;
        }
    }

    my @errs;
    while (my ($key,$value) = each %$values)
    {
        if ( !exists($$self{header}{INFO}{$key}) )
        {
            push @errs, "INFO tag [$key] not listed in the header" unless $$self{version}<3.3;
            my $nargs = defined $value ? -1 : 0;
            $self->_add_unknown_field('INFO',$key,$nargs);
            next;
        }
        my $type = $$self{header}{INFO}{$key};

        my @vals = defined $value ? split(/,/, $value) : ();
        if ( $$type{Number} eq 'G' )
        {
            if ( $ng != @vals && !(@vals==1 && $vals[0] eq '.') ) { push @errs, "INFO tag [$key=$value] expected different number of values (expected $ng, found ".scalar @vals.")"; }
        }
        elsif ( $$type{Number} eq 'A' )
        {
            if ( $na != @vals && !(@vals==1 && $vals[0] eq '.') ) { push @errs, "INFO tag [$key=$value] expected different number of values (expected $na, found ".scalar @vals.")"; }
        }
        elsif ( $$type{Number} eq 'R' )
        {
            if ( $nr != @vals && !(@vals==1 && $vals[0] eq '.') ) { push @errs, "INFO tag [$key=$value] expected different number of values (expected $nr, found ".scalar @vals.")"; }
        }
        elsif ( $$type{Number}==0 ) 
        {
            if ( defined($value) ) { push @errs, "INFO tag [$key] did not expect any parameters, got [$value]"; }
            next; 
        }
        elsif ( $$type{Number}!=-1 && @vals!=$$type{Number} )
        {
            if ( !(@vals==1 && $vals[0] eq '.') ) { push @errs, "INFO tag [$key=$value] expected different number of values ($$type{Number})"; }
        }
        if ( !$$type{handler} ) { next; }
        for my $val (@vals)
        {
            my $err = &{$$type{handler}}($self,$val,$$type{default});
            if ( $err ) { push @errs, $err; }
        }
    }
    if ( !@errs ) { return undef; }
    return join(',',@errs);
}

=head2 validate_gtype_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_gtype_field($$x{gtypes}{NA00001},$$x{ALT},$$x{FORMAT});
    Args    : The genotype data hashref
              The ALT arrayref
    Returns : Error message in case of an error.

=cut

sub guess_ploidy
{
    my ($self, $nals, $nvals) = @_;
    if ( $nvals==$nals ) { return 1; }
    if ( $nvals==binom(1+$nals,2) ) { return 2; }
    return -1;
}

sub binom
{
    my ($n, $k) = @_;
    my $b = 1;
    if ( $k > $n-$k ) { $k = $n-$k; }
    if ( $k < 1 ) { return 1; }
    for (my $i=1; $i<=$k; $i++) { $b *= ($n-$k+$i)/$i; }
    return $b;
}

sub validate_gtype_field
{
    my ($self,$data,$alts,$format) = @_;

    my @errs;
    my $ploidy = 2; 
    if ( !exists($$data{GT}) ) { push @errs, "The mandatory tag GT not present." unless $$self{ignore_missing_GT}; }
    else
    {
        my (@als) = $self->split_by($$data{GT},@{$$self{gt_sep}});
        for my $al (@als)
        {
            if ( $al eq '.' or $al eq '0' ) { next; }
            if ( !($al=~/^[0-9]+$/) ) { push @errs, "Unable to parse the GT field [$$data{GT}], expected integers"; }
            if ( !exists($$alts[$al-1]) ) { push @errs, "Bad ALT value in the GT field, the index [$al] out of bounds [$$data{GT}]."; last; }
        }
        $ploidy = @als;
    }

    # Expected numbers
    my $ng = -1;
    my $na = -1;
    my $nr = -1;
    if ( $$self{version}>4.0 )
    {
        if ( $$alts[0] eq '.' ) { $ng=1; $na=1; $nr=1; }
        else
        {
            $na = @$alts;
            $ng = binom($ploidy+$na,$ploidy);
            $nr = $na+1;
        }
    }

    while (my ($key,$value) = each %$data)
    {
        if ( !exists($$self{header}{FORMAT}{$key}) )
        {
            push @errs, "FORMAT tag [$key] not listed in the header" unless $$self{version}<3.3;
            $self->_add_unknown_field('FORMAT',$key,-1);
            next;
        }
        my $type = $$self{header}{FORMAT}{$key};

        my @vals = split(/,/, $value);
        if ( $$type{Number} eq 'G' )
        {
            if ( $ng != @vals && !(@vals==1 && $vals[0] eq '.') ) { push @errs, "FORMAT tag [$key] expected different number of values (expected $ng, found ".scalar @vals.")"; }
        }
        elsif ( $$type{Number} eq 'A' )
        {
            if ( $na != @vals && !(@vals==1 && $vals[0] eq '.') ) { push @errs, "FORMAT tag [$key] expected different number of values (expected $na, found ".scalar @vals.")"; }
        }
        elsif ( $$type{Number} eq 'R' )
        {
            if ( $nr != @vals && !(@vals==1 && $vals[0] eq '.') ) { push @errs, "FORMAT tag [$key] expected different number of values (expected $nr, found ".scalar @vals.")"; }
        }
        elsif ( $$type{Number}!=-1 && @vals!=$$type{Number} )
        {
            if ( !(@vals==1 && $vals[0] eq '.') ) { push @errs, "FORMAT tag [$key] expected different number of values ($$type{Number})"; }
        }
        if ( !$$type{handler} ) { next; }
        for my $val (@vals)
        {
            my $err = &{$$type{handler}}($self,$val,$$type{default});
            if ( $err ) { push @errs, $err; }
        }
    }
    if ( !@errs ) { return undef; }
    return join(',',@errs);
}


sub validate_ref_field
{
    my ($self,$ref) = @_;
    if ( !($ref=~/^[ACGTN]$/) ) { return "Expected one of A,C,G,T,N, got [$ref]\n"; }
    return undef;
}

sub validate_int
{
    my ($self,$value,$default) = @_;

    if ( defined($default) && $value eq $default ) { return undef; }
    if ( $value =~ /^-?\d+$/ ) { return undef; }
    return "Could not validate the int [$value]";
}

sub validate_float
{
    my ($self,$value,$default) = @_;
    if ( defined($default) && $value eq $default ) { return undef; }
    if ( $value =~ /^-?\d+(?:\.\d*)$/ ) { return undef; }
    if ( $value =~ /^-?\d*(?:\.\d+)$/ ) { return undef; }
    if ( $value =~ /^-?\d+$/ ) { return undef; }
    if ( $value =~ /^-?\d*(?:\.?\d+)(?:[Ee][-+]?\d+)?$/ ) { return undef; }
    return "Could not validate the float [$value]";
}

sub validate_char
{
    my ($self,$value,$default) = @_;

    if ( defined($default) && $value eq $default ) { return undef; }
    if ( length($value)==1) { return undef; }
    return "Could not validate the char value [$value]";
}


=head2 run_validation

    About   : Validates the VCF file.
    Usage   : my $vcf = Vcf->new(file=>'file.vcf'); $vcf->run_validation('example.vcf.gz');
    Args    : File name or file handle.

=cut

sub run_validation
{
    my ($self) = @_;

    $self->parse_header();
    $self->validate_header();

    if ( !exists($$self{header}) ) { $self->warn(qq[The header not present.\n]); }
    elsif ( !exists($$self{header}{fileformat}) ) 
    {
        $self->warn(qq[The "fileformat" field not present in the header, assuming VCFv$$self{version}\n]);
    }
    elsif ( $$self{header_lines}[0]{key} ne 'fileformat' ) 
    {
        $self->warn(qq[The "fileformat" not the first line in the header\n]);
    }
    if ( !exists($$self{columns}) ) { $self->warn("No column descriptions found.\n"); }

    my $default_qual = $$self{defaults}{QUAL};
    my $warn_sorted=1;
    my $warn_duplicates = exists($$self{warn_duplicates}) ? $$self{warn_duplicates} : 1;
    my ($prev_chrm,$prev_pos);
    while (my $line=$self->next_data_array()) 
    {
        for (my $i=0; $i<@$line; $i++)
        {
            if (!defined($$line[$i]) or $$line[$i] eq '' ) 
            {
                my $colname = $i<@{$$self{columns}} ? $$self{columns}[$i] : $i+1;
                $self->warn("The column $colname is empty at $$line[0]:$$line[1].\n");
            }
        }

        my $x = $self->next_data_hash($line);
        $self->validate_line($x);

        # Is the position numeric?
        if ( !($$x{POS}=~/^\d+$/) ) { $self->warn("Expected integer for the position at $$x{CHROM}:$$x{POS}\n"); }

        if ( $warn_duplicates )
        {
            if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos eq $$x{POS} )
            {
                $self->warn("Warning: Duplicate entries, for example $$x{CHROM}:$$x{POS}\n");
                $warn_duplicates = 0;
            }
        }

        # Is the file sorted?
        if ( $warn_sorted )
        {
            if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos > $$x{POS} ) 
            { 
                $self->warn("Warning: The file is not sorted, for example $$x{CHROM}:$$x{POS} comes after $prev_chrm:$prev_pos\n");
                $warn_sorted = 0;
            }
            $prev_chrm = $$x{CHROM};
            $prev_pos  = $$x{POS};
        }

        # The reference base: one of A,C,G,T,N, non-empty.
        my $err = $self->validate_ref_field($$x{REF});
        if ( $err ) { $self->warn("$$x{CHROM}:$$x{POS} .. $err\n"); }
        
        # The ALT field (alternate non-reference base)
        $err = $self->validate_alt_field($$x{ALT},$$x{REF});
        if ( $err ) { $self->warn("$$x{CHROM}:$$x{POS} .. $err\n"); }

        # The QUAL field
        my $ret = $self->validate_float($$x{QUAL},$default_qual);
        if ( $ret ) { $self->warn("QUAL field at $$x{CHROM}:$$x{POS} .. $ret\n"); }
        elsif ( $$x{QUAL}=~/^-?\d+$/ && $$x{QUAL}<-1 ) { $self->warn("QUAL field at $$x{CHROM}:$$x{POS} is negative .. $$x{QUAL}\n"); }

        # The FILTER field
        $err = $self->validate_filter_field($$x{FILTER});
        if ( $err ) { $self->warn("FILTER field at $$x{CHROM}:$$x{POS} .. $err\n"); }

        # The INFO field
        $err = $self->validate_info_field($$x{INFO},$$x{ALT});
        if ( $err ) { $self->warn("INFO field at $$x{CHROM}:$$x{POS} .. $err\n"); } 

        while (my ($gt,$data) = each %{$$x{gtypes}})
        {
            $err = $self->validate_gtype_field($data,$$x{ALT},$$x{FORMAT});
            if ( $err ) { $self->warn("column $gt at $$x{CHROM}:$$x{POS} .. $err\n"); }
        }

        if ( scalar keys %{$$x{gtypes}} && (exists($$x{INFO}{AN}) || exists($$x{INFO}{AC})) )
        {
            my $nalt = scalar @{$$x{ALT}};
            if ( $nalt==1 && $$x{ALT}[0] eq '.' ) { $nalt=0; }
            my ($an,$ac) = $self->calc_an_ac($$x{gtypes},$nalt);    # Allow alleles in ALT which are absent in samples
            if ( exists($$x{INFO}{AN}) && $an ne $$x{INFO}{AN} ) 
            { 
                $self->warn("$$x{CHROM}:$$x{POS} .. AN is $$x{INFO}{AN}, should be $an\n"); 
            }
            if ( exists($$x{INFO}{AC}) && $ac ne $$x{INFO}{AC} ) 
            { 
                $self->warn("$$x{CHROM}:$$x{POS} .. AC is $$x{INFO}{AC}, should be $ac\n"); 
            }
        }
    }
}


=head2 get_chromosomes

    About   : Get list of chromosomes from the VCF file. Must be bgzipped and tabix indexed.
    Usage   : my $vcf = Vcf->new(); $vcf->get_chromosomes();
    Args    : none

=cut

sub get_chromosomes
{
    my ($self) = @_;
    if ( !$$self{file} ) { $self->throw(qq[The parameter "file" not set.\n]); }
    my (@out) = `tabix -l '$$self{file}'`;
    if ( $? ) 
    { 
        my @has_tabix = `which tabix`;
        if ( !@has_tabix ) { $self->throw(qq[The command "tabix" not found, please add it to your PATH\n]); }
        $self->throw(qq[The command "tabix -l $$self{file}" exited with an error. Is the file tabix indexed?\n]); 
    }
    for (my $i=0; $i<@out; $i++) { chomp($out[$i]); }
    return \@out;
}


=head2 get_samples

    About   : Get list of samples.
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header(); my (@samples) = $vcf->get_samples();
    Args    : none

=cut

sub get_samples
{
    my ($self) = @_;
    my $n = @{$$self{columns}} - 1;
    return (@{$$self{columns}}[9..$n]);
}


=head2 get_column

    About   : Convenient way to get data for a sample
    Usage   : my $rec = $vcf->next_data_array(); my $sample_col = $vcf->get_column($rec, 'NA0001');
    Args 1  : Array pointer returned by next_data_array
         2  : Column/Sample name

=cut

sub get_column
{
    my ($self,$line,$column) = @_;
    if ( !exists($$self{has_column}{$column}) ) { $self->throw("No such column: [$column]\n"); }
    my $idx = $$self{has_column}{$column};
    return $$line[$idx-1];
}

=head2 get_column_name

    About   : Mapping between zero-based VCF column and its name
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header(); my $name = $vcf->get_column_name(1); # returns POS
    Args    : Index of the column (0-based)

=cut

sub get_column_name
{
    my ($self,$idx) = @_;
    if ( $idx >= @{$$self{columns}} ) { $self->throw("The index out of bounds\n"); }
    return $$self{columns}[$idx];
}

=head2 get_column_index

    About   : Mapping between VCF column name and its zero-based index
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header(); my $name = $vcf->get_column_index('POS'); # returns 1
    Args    : Name of the column

=cut

sub get_column_index
{
    my ($self,$column) = @_;
    if ( !exists($$self{has_column}{$column}) ) { $self->throw("No such column: [$column]\n"); }
    return $$self{has_column}{$column}-1;
}


#------------------------------------------------
# Version 3.2 specific functions

package Vcf3_2;
use base qw(VcfReader);

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    bless $self, ref($class) || $class;

    $$self{_defaults} =
    {
        version => '3.2',
        drop_trailings => 1,
        filter_passed  => 0,

        defaults =>
        {
            QUAL    => '-1',
            default => '.',
            Flag    => undef,
            GT      => '.',
        },

        handlers =>
        {
            Integer   => \&VcfReader::validate_int,
            Float     => \&VcfReader::validate_float,
            Character => \&VcfReader::validate_char,
            String    => undef,
            Flag      => undef,
        },

        regex_snp   => qr/^[ACGTN]$/i,
        regex_ins   => qr/^I[ACGTN]+$/,
        regex_del   => qr/^D\d+$/,
        regex_gtsep => qr{[\\|/]},
        regex_gt    => qr{^(\.|\d+)([\\|/]?)(\.?|\d*)$},
        regex_gt2   => qr{^(\.|[0-9ACGTNIDacgtn]+)([\\|/]?)}, 
    };

    for my $key (keys %{$$self{_defaults}}) 
    { 
        $$self{$key}=$$self{_defaults}{$key}; 
    }


    return $self;
}


#------------------------------------------------
# Version 3.3 specific functions

package Vcf3_3;
use base qw(VcfReader);

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    bless $self, ref($class) || $class;

    $$self{_defaults} = 
    {
        version => '3.3',
        drop_trailings => 0,
        filter_passed  => 0,

        defaults =>
        {
            QUAL      => '-1',
            Integer   => '-1',
            Float     => '-1',
            Character => '.',
            String    => '.',
            Flag      => undef,
            GT        => './.',
            default   => '.',
        },

        handlers =>
        {
            Integer   => \&VcfReader::validate_int,
            Float     => \&VcfReader::validate_float,
            Character => \&VcfReader::validate_char,
            String    => undef,
            Flag      => undef,
        },

        regex_snp   => qr/^[ACGTN]$/i,
        regex_ins   => qr/^I[ACGTN]+$/,
        regex_del   => qr/^D\d+$/,
        regex_gtsep => qr{[\\|/]},
        regex_gt    => qr{^(\.|\d+)([\\|/]?)(\.?|\d*)$},
        regex_gt2   => qr{^(\.|[0-9ACGTNIDacgtn]+)([\\|/]?)}, # . 0/1 0|1 A/A A|A D4/IACGT
        gt_sep => [qw(\ | /)],
    };

    for my $key (keys %{$$self{_defaults}}) 
    { 
        $$self{$key}=$$self{_defaults}{$key}; 
    }

    return $self;
}


#------------------------------------------------
# Version 4.0 specific functions

=head1 VCFv4.0

VCFv4.0 specific functions

=cut

package Vcf4_0;
use base qw(VcfReader);

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    bless $self, ref($class) || $class;

    $$self{_defaults} = 
    {
        version => '4.0',
        drop_trailings => 1,
        filter_passed  => 'PASS',

        defaults => 
        {
            QUAL    => '.',
            Flag    => undef,
            GT      => '.',
            default => '.',
        },
        reserved => 
        {
            FILTER  => { 0=>1 },
        },

        handlers =>
        {
            Integer    => \&VcfReader::validate_int,
            Float      => \&VcfReader::validate_float,
            Character  => \&VcfReader::validate_char,
            String     => undef,
            Flag      => undef,
        },

        regex_snp   => qr/^[ACGTN]$|^<[\w:.]+>$/i,
        regex_ins   => qr/^[ACGTN]+$/,
        regex_del   => qr/^[ACGTN]+$/,
        regex_gtsep => qr{[|/]},                     # | /
        regex_gt    => qr{^(\.|\d+)([|/]?)(\.?|\d*)$},   # . ./. 0/1 0|1
        regex_gt2   => qr{^(\.|[0-9ACGTNacgtn]+|<[\w:.]+>)([|/]?)},   # . ./. 0/1 0|1 A/A A|A 0|<DEL:ME:ALU>
        gt_sep => [qw(| /)],
    };

    for my $key (keys %{$$self{_defaults}}) 
    { 
        $$self{$key}=$$self{_defaults}{$key}; 
    }

    return $self;
}

sub Vcf4_0::format_header_line
{
    my ($self,$rec) = @_;

    my %tmp_rec = ( %$rec );
    if ( exists($tmp_rec{Number}) && $tmp_rec{Number} eq '-1' ) { $tmp_rec{Number} = '.' }
    my $value;
    if ( exists($tmp_rec{ID}) or $tmp_rec{key} eq 'PEDIGREE' )
    {
        my %has = ( key=>1, handler=>1, default=>1 );   # Internal keys not to be output
        my @items;
        for my $key (qw(ID Number Type Description), sort keys %tmp_rec)
        {
            if ( !exists($tmp_rec{$key}) or $has{$key} ) { next; }
            my $quote = ($key eq 'Description' or $tmp_rec{$key}=~/\s/) ? '"' : '';
            push @items, "$key=$quote$tmp_rec{$key}$quote";
            $has{$key}=1;
        }
        $value = '<' .join(',',@items). '>';
    }
    else { $value = $tmp_rec{value}; }

    my $line = "##$tmp_rec{key}=".$value."\n";
    return $line;
}

=head2 parse_header_line

    Usage   : $vcf->parse_header_line(q[##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">])
              $vcf->parse_header_line(q[reference=1000GenomesPilot-NCBI36])
    Args    : 
    Returns : 

=cut

sub Vcf4_0::parse_header_line
{
    my ($self,$line) = @_;

    chomp($line);
    $line =~ s/^##//;
    
    if ( !($line=~/^([^=]+)=/) ) { $self->throw("Expected key=value pair in the header: $line\n"); }
    my $key   = $1;
    my $value = $';

    if ( !($value=~/^<(.+)>\s*$/) ) 
    { 
        # Simple sanity check for subtle typos
        if ( $key eq 'INFO' or $key eq 'FILTER' or $key eq 'FORMAT' or $key eq 'ALT' )
        {
            $self->throw("Hmm, is this a typo? [$key] [$value]");
        }
        return { key=>$key, value=>$value }; 
    }

    my $rec = { key=>$key };
    my $tmp = $1;
    my ($attr_key,$attr_value,$quoted);
    while ($tmp ne '')
    {
        if ( !defined $attr_key )
        {
            if ( $tmp=~/^([^=]+)="/ ) { $attr_key=$1; $quoted=1; $tmp=$'; next; }
            elsif ( $tmp=~/^([^=]+)=/ ) { $attr_key=$1; $quoted=0; $tmp=$'; next; }
            else { $self->throw(qq[Could not parse header line: $line\nStopped at [$tmp].\n]); }
        }

        if ( $tmp=~/^[^,\\"]+/ ) { $attr_value .= $&; $tmp = $'; }
        if ( $tmp=~/^\\\\/ ) { $attr_value .= '\\\\'; $tmp = $'; next; }
        if ( $tmp=~/^\\"/ ) { $attr_value .= '\\"'; $tmp = $'; next; }
        if ( $tmp eq '' or ($tmp=~/^,/ && !$quoted) or $tmp=~/^"/ )
        {
            if ( $attr_key=~/^\s+/ or $attr_key=~/\s+$/ or $attr_value=~/^\s+/ or $attr_value=~/\s+$/ ) 
            { 
                $self->warn("Leading or trailing space in attr_key-attr_value pairs is discouraged:\n\t[$attr_key] [$attr_value]\n\t$line\n"); 
                $attr_key =~ s/^\s+//;
                $attr_key =~ s/\s+$//;
                $attr_value =~ s/^\s+//;
                $attr_value =~ s/\s+$//;
            }
            $$rec{$attr_key} = $attr_value;
            $tmp = $';
            if ( $quoted && $tmp=~/^,/ ) { $tmp = $'; }
            $attr_key = $attr_value = $quoted = undef;
            next;
        }
        if ( $tmp=~/^,/ ) { $attr_value .= $&; $tmp = $'; next; }
        $self->throw(qq[Could not parse header line: $line\nStopped at [$tmp].\n]);
    }

    if ( $key eq 'INFO' or $key eq 'FILTER' or $key eq 'FORMAT' )
    {
        if ( $key ne 'PEDIGREE' && !exists($$rec{ID}) ) { $self->throw("Missing the ID tag in $line\n"); }
        if ( !exists($$rec{Description}) ) { $self->warn("Missing the Description tag in $line\n"); }
    }
    if ( exists($$rec{Number}) && $$rec{Number} eq '-1' ) { $self->warn("The use of -1 for unknown number of values is deprecated, please use '.' instead.\n\t$line\n"); }
    if ( exists($$rec{Number}) && $$rec{Number} eq '.' ) { $$rec{Number}=-1; }

    return $rec;
}

sub Vcf4_0::validate_ref_field
{
    my ($self,$ref) = @_;
    if ( !($ref=~/^[ACGTN]+$/) ) 
    {
        my $offending = $ref;
        $offending =~ s/[ACGTN]+//g;
        return "Expected combination of A,C,G,T,N for REF, got [$ref], the offending chars were [$offending]\n"; 
    }
    return undef;
}

sub Vcf4_0::validate_alt_field
{
    my ($self,$values,$ref) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }

    my $ret = $self->_validate_alt_field($values,$ref);
    if ( $ret ) { return $ret; }

    my $ref_len = length($ref);
    my $ref1 = substr($ref,0,1);

    my @err;
    my $msg = '';
    for my $item (@$values)
    {
        if ( !($item=~/^[ACTGN]+$|^<[^<>\s]+>$/) ) { push @err,$item; next; }
        if ( $item=~/^<[^<>\s]+>$/ ) { next; }
        if ( $ref_len==length($item) ) { next; }
        if ( substr($item,0,1) ne $ref1 ) { $msg=', first base does not match the reference.'; push @err,$item; next; }
    }
    if ( !@err ) { return undef; }
    return 'Could not parse the allele(s) [' .join(',',@err). ']' . $msg;
}


=head2 fill_ref_alt_mapping

    About   : A tool for merging VCFv4.0 records. The subroutine unifies the REFs and creates a mapping
                from the original haplotypes to the haplotypes based on the new REF. Consider the following
                example:
                    REF ALT
                    G    GA
                    GT   G
                    GT   GA
                    GT   GAA
                    GTC  G
                    G    <DEL>
                my $map={G=>{GA=>1},GT=>{G=>1,GA=>1,GAA=>1},GTC=>{G=>1},G=>{'<DEL>'=>1}};
                my $new_ref=$vcf->fill_ref_alt_mapping($map);
                
              The call returns GTC and $map is now
                    G    GA     ->      GTC  GATC
                    GT   G      ->      GTC  GC
                    GT   GA     ->      GTC  GAC
                    GT   GAA    ->      GTC  GAAC
                    GTC  G      ->      GTC  G
                    G    <DEL>  ->      GTC  <DEL>
    Args    : 
    Returns : New REF string and fills the hash with appropriate ALT or undef on error.

=cut

sub Vcf4_0::fill_ref_alt_mapping
{
    my ($self,$map) = @_;
    
    my $max_len = 0;
    my $new_ref;
    for my $ref (keys %$map)
    {
        my $len = length($ref);
        if ( $max_len<$len ) 
        { 
            $max_len = $len; 
            $new_ref = $ref;
        }
        $$map{$ref}{$ref} = 1;
    }
    for my $ref (keys %$map)
    {
        my $rlen = length($ref);
        if ( substr($new_ref,0,$rlen) ne $ref ) { $self->warn("The reference prefixes do not agree: $ref vs $new_ref\n"); return undef; }
        for my $alt (keys %{$$map{$ref}})
        {
            # The second part of the regex is for VCF>4.0, but does no harm for v<=4.0
            if ( $alt=~/^<.+>$/ or $alt=~/\[|\]/ ) { $$map{$ref}{$alt} = $alt; next; }
            my $new = $alt;
            if ( $rlen<$max_len ) { $new .= substr($new_ref,$rlen); }
            $$map{$ref}{$alt} = $new;
        }
    }
    return $new_ref;
}

=head2 normalize_alleles

    About   : Makes REF and ALT alleles more compact if possible (e.g. TA,TAA -> T,TA)
    Usage   : my $line = $vcf->next_data_array();
              ($ref,@alts) = $vcf->normalize_alleles($$line[3],$$line[4]);

=cut

sub Vcf4_0::normalize_alleles
{
    my ($self,$ref,$alt) = @_;

    my $rlen = length($ref);
    if ( $rlen==1 or length($alt)==1 )  { return ($ref,split(/,/,$alt)); }

    my @als = split(/,/,$alt);
    my $i = 1;
    my $done = 0;
    while ( $i<$rlen )
    {
        my $r = substr($ref,$rlen-$i,1);
        for my $al (@als)
        {
            my $len = length($al);
            if ( $i>=$len ) { $done = 1; }
            my $c = substr($al,$len-$i,1);
            if ( $c ne $r ) { $done = 1; last; }
        }
        if ( $done ) { last; }
        $i++;
    }
    if ( $i>1 )
    {
        $i--;
        $ref = substr($ref,0,$rlen-$i);
        for (my $j=0; $j<@als; $j++) { $als[$j] = substr($als[$j],0,length($als[$j])-$i); }
    }
    return ($ref,@als);
}

sub Vcf4_0::normalize_alleles_pos
{
    my ($self,$ref,$alt) = @_;
    my @als;
    ($ref,@als) = $self->normalize_alleles($ref,$alt);

    my $rlen = length($ref);
    if ( $rlen==1 ) { return (0,$ref,@als); }
    my $i = 0;
    my $done = 0;
    while ( $i+1<$rlen )
    {
        my $r = substr($ref,$i,1);
        for my $al (@als)
        {
            my $len = length($al);
            if ( $i+1>=$len ) { $done = 1; last; }
            my $c = substr($al,$i,1);
            if ( $c ne $r ) { $done = 1; last; }
        }
        if ( $done ) { last; }
        $i++;
    }
    if ( $i<0 ) { $i = 0; }
    if ( $i>0 )
    {
        substr($ref,0,$i,'');
        for (my $j=0; $j<@als; $j++) { substr($als[$j],0,$i,''); }
    }
    return ($i,$ref,@als);
}

sub Vcf4_0::event_type
{
    my ($self,$rec,$allele) = @_;

    my $ref = $rec;
    if ( ref($rec) eq 'HASH' ) 
    { 
        if ( exists($$rec{_cached_events}{$allele}) ) { return (@{$$rec{_cached_events}{$allele}}); }
        $ref = $$rec{REF};
    }

    if ( $allele=~/^<[^>]+>$/ ) 
    { 
        if ( ref($rec) eq 'HASH' ) { $$rec{_cached_events}{$allele} = ['u',0,$allele]; }
        return ('u',0,$allele); 
    }
    if ( $allele eq '.' )
    {
        if ( ref($rec) eq 'HASH' ) { $$rec{_cached_events}{$allele} = ['r',0,$ref]; }
        return ('r',0,$ref);
    }

    my $reflen = length($ref);
    my $len = length($allele);
    my $ht;
    my $type;
    if ( $len==$reflen )
    {
        # This can be a reference, a SNP, or multiple SNPs
        my $mism = 0;
        for (my $i=0; $i<$len; $i++)
        {
            if ( substr($ref,$i,1) ne substr($allele,$i,1) ) { $mism++; }
        }
        if ( $mism==0 ) { $type='r'; $len=0; }
        else { $type='s'; $len=$mism; }
    }
    else
    {
        ($len,$ht)=$self->is_indel($ref,$allele);
        if ( $len )
        {
            # Indel
            $type = 'i';
            $allele = $ht;
        }
        else 
        {
            $type = 'o'; $len = $len>$reflen ? $len-1 : $reflen-1;
        }
    }

    if ( ref($rec) eq 'HASH' )
    {
        $$rec{_cached_events}{$allele} = [$type,$len,$allele];
    }
    return ($type,$len,$allele);
}

# The sequences start at the same position, which simplifies things greatly.
# Returns length of the indel (+ insertion, - deletion), the deleted/inserted sequence
#   and the position of the first base after the shared sequence
sub is_indel
{
    my ($self,$seq1,$seq2) = @_;

    my $len1 = length($seq1);
    my $len2 = length($seq2);
    if ( $len1 eq $len2 ) { return (0,'',0); }

    my ($del,$len,$LEN);
    if ( $len1<$len2 )
    {
        $len = $len1;
        $LEN = $len2;
        $del = 1;
    }
    else
    {
        $len = $len2;
        $LEN = $len1;
        $del = -1;
        my $tmp=$seq1; $seq1=$seq2; $seq2=$tmp;
    }

    my $ileft;
    for ($ileft=0; $ileft<$len; $ileft++)
    {
        if ( substr($seq1,$ileft,1) ne substr($seq2,$ileft,1) ) { last; }
    }
    if ( $ileft==$len )
    {
        return ($del*($LEN-$len), substr($seq2,$ileft), $ileft);
    }

    my $iright;
    for ($iright=0; $iright<$len; $iright++)
    {
        if ( substr($seq1,$len-$iright,1) ne substr($seq2,$LEN-$iright,1) ) { last; }
    }
    if ( $iright+$ileft<=$len ) { return (0,'',0); }

    return ($del*($LEN-$len),substr($seq2,$ileft,$LEN-$len),$ileft);
}


#------------------------------------------------
# Version 4.1 specific functions

=head1 VCFv4.1

VCFv4.1 specific functions

=cut

package Vcf4_1;
use base qw(Vcf4_0);

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    bless $self, ref($class) || $class;

    $$self{_defaults} = 
    {
        version => '4.1',
        drop_trailings => 1,
        filter_passed  => 'PASS',

        defaults => 
        {
            QUAL    => '.',
            Flag    => undef,
            GT      => '.',
            default => '.',
        },
        reserved => 
        {
            FILTER  => { 0=>1 },
        },

        handlers =>
        {
            Integer    => \&VcfReader::validate_int,
            Float      => \&VcfReader::validate_float,
            Character  => \&VcfReader::validate_char,
            String     => undef,
            Flag      => undef,
        },

        regex_snp   => qr/^[ACGTN]$|^<[\w:.]+>$/i,
        regex_ins   => qr/^[ACGTN]+$/i,
        regex_del   => qr/^[ACGTN]+$/i,
        regex_gtsep => qr{[|/]},                     # | /
        regex_gt    => qr{^(\.|\d+)([|/]?)(\.?|\d*)$},   # . ./. 0/1 0|1
        regex_gt2   => qr{^(\.|[0-9ACGTNacgtn]+|<[\w:.]+>)([|/]?)},   # . ./. 0/1 0|1 A/A A|A 0|<DEL:ME:ALU>
        gt_sep => [qw(| /)],
    };

    $$self{ignore_missing_GT} = 1;

    for my $key (keys %{$$self{_defaults}}) 
    { 
        $$self{$key}=$$self{_defaults}{$key}; 
    }

    return $self;
}

sub Vcf4_1::validate_header
{
    my ($self) = @_;
    my $lines = $self->get_header_line(key=>'reference');
    if ( !@$lines ) { $self->warn("The header tag 'reference' not present. (Not required but highly recommended.)\n"); }
}

sub Vcf4_1::validate_line
{
    my ($self,$line) = @_;

    if ( !$$self{_contig_validated}{$$line{CHROM}} )
    {
        my $lines = $self->get_header_line(key=>'contig',ID=>$$line{CHROM});
        if ( !@$lines ) { $self->warn("The header tag 'contig' not present for CHROM=$$line{CHROM}. (Not required but highly recommended.)\n"); }
        $$self{_contig_validated}{$$line{CHROM}} = 1;
    }

    if ( index($$line{CHROM},':')!=-1 ) { $self->warn("Colons not allowed in chromosome names: $$line{CHROM}\n"); }

    # Is the ID composed of alphanumeric chars
    if ( !($$line{ID}=~/^\S+$/) ) { $self->warn("Expected non-whitespace ID at $$line{CHROM}:$$line{POS}, but got [$$line{ID}]\n"); }
}

sub Vcf4_1::validate_alt_field
{
    my ($self,$values,$ref) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }

    my $ret = $self->_validate_alt_field($values,$ref);
    if ( $ret ) { return $ret; }

    my $ref_len = length($ref);
    my $ref1 = substr($ref,0,1);

    my @err;
    my $msg = '';
    for my $item (@$values)
    {
        if ( $item=~/^(.*)\[(.+)\[(.*)$/ or $item=~/^(.*)\](.+)\](.*)$/ )
        {
            if ( $1 ne '' && $3 ne '' ) { $msg=', two replacement strings given (expected one)'; push @err,$item; next; }
            my $rpl;
            if ( $1 ne '' )
            {
                $rpl  = $1;
                if ( $rpl ne '.' )
                {
                    my $rref = substr($rpl,0,1);
                    if ( $rref ne $ref1 ) { $msg=', the first base of the replacement string does not match the reference'; push @err,$item; next; }
                }
            }
            else
            {
                $rpl  = $3;
                if ( $rpl ne '.' )
                {
                    my $rref = substr($rpl,-1,1);
                    if ( $rref ne $ref1 ) { $msg=', the last base of the replacement string does not match the reference'; push @err,$item; next; }
                }
            }
            my $pos = $2;
            if ( !($rpl=~/^[ACTGNacgtn]+$/) && $rpl ne '.' ) { $msg=', replacement string not valid (expected [ACTGNacgtn]+)'; push @err,$item; next; }
            if ( !($pos=~/^\S+:\d+$/) ) { $msg=', cannot parse sequence:position'; push @err,$item; next; }
            next;
        }
        if ( $item=~/^\.[ACTGNactgn]*([ACTGNactgn])$/ ) { next; }
        elsif ( $item=~/^([ACTGNactgn])[ACTGNactgn]*\.$/ ) { next; }
        if ( !($item=~/^[ACTGNactgn]+$|^<[^<>\s]+>$/) ) { push @err,$item; next; }
    }
    if ( !@err ) { return undef; }
    return 'Could not parse the allele(s) [' .join(',',@err). ']' . $msg;
}

sub Vcf4_1::next_data_hash
{
    my ($self,@args) = @_;

    my $out = $self->SUPER::next_data_hash(@args);
    if ( !defined $out or $$self{assume_uppercase} ) { return $out; }

    # Case-insensitive ALT and REF bases
    $$out{REF} = uc($$out{REF});
    my $nalt = @{$$out{ALT}};
    for (my $i=0; $i<$nalt; $i++)
    {
        if ( $$out{ALT}[$i]=~/^</ ) { next; }
        $$out{ALT}[$i] = uc($$out{ALT}[$i]);
    }

    return $out;
}

sub Vcf4_1::next_data_array
{
    my ($self,@args) = @_;

    my $out = $self->SUPER::next_data_array(@args);
    if ( !defined $out or $$self{assume_uppercase} ) { return $out; }

    # Case-insensitive ALT and REF bases
    $$out[3] = uc($$out[3]);
    my $alt  = $$out[4];
    $$out[4] = '';
    my $pos = 0;
    while ( $pos<length($alt) && (my $start=index($alt,'<',$pos))!=-1 )
    {
        my $end = index($alt,'>',$start+1);
        if ( $end==-1 ) { $self->throw("Could not parse ALT [$alt]\n") }
        if ( $start>$pos )
        {
            $$out[4] .= uc(substr($alt,$pos,$start-$pos));
        }
        $$out[4] .= substr($alt,$start,$end-$start+1);
        $pos = $end+1;
    }
    if ( $pos<length($alt) )
    {
        $$out[4] .= uc(substr($alt,$pos));
    }
    return $out;
}

sub Vcf4_1::event_type
{
    my ($self,$rec,$allele) = @_;

    my $len = length($allele);
    if ( $len==1 ) { return $self->SUPER::event_type($rec,$allele); }

    my $c = substr($allele,0,1);
    if ( $c eq '<' ) { return ('u',0,$allele); }
    elsif ( $c eq '[' or $c eq ']' or $c eq '.' ) { return 'b'; }

    $c = substr($allele,-1,1);
    if ( $c eq '[' or $c eq ']' or $c eq '.' ) { return 'b'; }
    elsif ( index($allele,'[')!=-1 or index($allele,']')!=-1 ) { return 'b'; }

    return $self->SUPER::event_type($rec,$allele);
}

#------------------------------------------------
# Version 4.2 specific functions

=head1 VCFv4.2

VCFv4.2 specific functions

=cut

package Vcf4_2;
use base qw(Vcf4_1);

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    bless $self, ref($class) || $class;

    $$self{version} = '4.2';
    return $self;
}

1;

