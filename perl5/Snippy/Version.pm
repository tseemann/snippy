package Snippy::Version;

use base Exporter;
@EXPORT_OK = qw(version);
%EXPORT_TAGS = ( 'all' =>  [ @EXPORT_OK ] );

our $VERSION = "4.5.0";

use strict;
use File::Basename;

sub version {
  return $VERSION;
}

sub command {
  return basename($0);
}

sub version_string {
  my $self = shift;
  return $self->command.' '.$self->version;
}

1;
