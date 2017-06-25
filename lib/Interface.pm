#!/usr/bin/env perl
use 5.10.0;

package Interface;

use File::Basename;

use Mouse;

use Path::Tiny;
# use Types::Path::Tiny qw/Path File AbsFile AbsPath/;
use Mouse::Util::TypeConstraints;

use namespace::autoclean;

use DDP;

use YAML::XS qw/LoadFile/;


use Getopt::Long::Descriptive;

use Seq;
with 'MouseX::Getopt';

##########Parameters accepted from command line#################
has input_file => (
  is        => 'ro',
  isa         => 'Str',
  required      => 1,
  metaclass => 'Getopt',
  cmd_aliases   => [qw/input i in/],
  documentation => qq{Input file path.},
);

has output_file_base => (
  is          => 'ro',
  isa         => 'Str',
  cmd_aliases   => [qw/o out/],
  metaclass => 'Getopt',
  documentation => qq{Where you want your output.},
);

has config => (
  is          => 'ro',
  isa         => 'Str',
  coerce      => 1,
  required    => 1,
  metaclass => 'Getopt',
  cmd_aliases   => [qw/c config/],
  documentation => qq{Yaml config file path.},
);

has overwrite => (
  is          => 'ro',
  isa         => 'Int',
  default     => 0,
  required    => 0,
  metaclass => 'Getopt',
  documentation => qq{Overwrite existing output file.},
);

has debug => (
  is          => 'ro',
  isa         => 'Num',
  default     => 0,
  required    => 0,
  metaclass   => 'Getopt',
 );

has verbose => (
  is          => 'ro',
  isa         => 'Int',
  required    => 0,
  metaclass   => 'Getopt',
 );

has compress => (
  is => 'ro', 
  isa => 'Bool',
  metaclass   => 'Getopt',
  documentation =>
    qq{Compress the output?},
  default => 0,
);

has archive => (
  is => 'ro',
  isa => 'Bool',
  metaclass   => 'Getopt',
  documentation =>
    qq{Place all outputs into a tarball?},
  default => 0,
);

has run_statistics => (
  is => 'ro', 
  isa => 'Int',
  metaclass   => 'Getopt',
  documentation =>
    qq{Create per-sample feature statistics (like transition:transversions)?},
  default => 1,
);

has delete_temp => (
  is => 'ro',
  isa => 'Int',
  documentation =>
    qq{Delete the temporary directory made during annotation},
  default => 1,
);

has wantedChr => (
  is => 'ro',
  isa => 'Str',
  metaclass => 'Getopt',
  cmd_aliases => [qw/chr wanted_chr/],
  documentation =>
    qq{Annotate a single chromosome},
);

has max_threads => (
  is => 'ro',
  isa => 'Int',
  metaclass => 'Getopt',
  documentation =>
    qq{Number of CPU threads to use (optional)},
);

subtype HashRefJson => as 'HashRef'; #subtype 'HashRefJson', as 'HashRef', where { ref $_ eq 'HASH' };
coerce HashRefJson => from 'Str' => via { from_json $_ };
subtype ArrayRefJson => as 'ArrayRef';
coerce ArrayRefJson => from 'Str' => via { from_json $_ };

has publisher => (
  is => 'ro',
  isa => 'HashRefJson',
  coerce => 1,
  required => 0,
  metaclass   => 'Getopt',
  documentation => 
    qq{Tell Seqant how to send messages to a plugged-in interface 
      (such as a web interface) }
);

has ignore_unknown_chr => (
  is          => 'ro',
  isa         => 'Bool',
  default     => 1,
  required    => 0,
  metaclass   => 'Getopt',
  documentation =>
    qq{Don't quit if we find a non-reference chromosome (like ChrUn)}
);

sub annotate {
  my $self = shift;
  
  my $args = {
    config => $self->config,
    input_file => $self->input_file,
    output_file_base => $self->output_file_base,
    debug => $self->debug,
    wantedChr => $self->wantedChr,
    ignore_unknown_chr => $self->ignore_unknown_chr,
    overwrite => $self->overwrite,
    publisher => $self->publisher,
    compress => $self->compress,
    archive => $self->archive,
    run_statistics => !!$self->run_statistics,
    delete_temp => !!$self->delete_temp,
  };

  if(defined $self->verbose) {
    $args->{verbose} = $self->verbose;
  }

  if(defined $self->max_threads) {
    $args->{max_threads} = $self->max_threads;
  }

  my $annotator = Seq->new_with_config($args);
  $annotator->annotate();
}

__PACKAGE__->meta->make_immutable;

1;

=item messanger

Contains a hash reference (also accept json representation of hash) that 
tells Seqant how to send data to a plugged interface.

Example: {
      room: jobObj.userID,
      message: {
        publicID: jobObj.publicID,
        data: tData,
      },
    };
=cut


# sub _run {
#   my $self = shift;

#   if ( $self->isProkaryotic ) {
#     my $args = "--vcf " . $self->snpfile . " --gb " . $self->genBankAnnotation;

#     system( $self->_prokAnnotatorPath . " " . $args );
#   }
#   else {
#     my $aInstance = Seq->new( $self->_annotatorArgsHref );
#     $aInstance->annotate_snpfile();
#   }
# }

###optional

# has genBankAnnotation => (
#   metaclass   => 'Getopt',
#   is          => 'ro',
#   isa         => 'Str',
#   cmd_aliases => [qw/gb g gen_bank_annotation/],
#   required    => 0,
#   documentation =>
#     qq{GenBank Annotation file path. Required for prokaryotic annotations. Type Str.},
#   predicate => 'isProkaryotic'
# );


# has serverMode  => (
#   metaclass => 'Getopt',
#   is => 'ro',
#   isa => 'Bool',
#   cmd_aliases => 'qw/s server/',
#   required => 0,
#   default => 0,
#   documentation => qq{Enables persistent server mode}
# );

#private vars

# has _prokAnnotatorPath => (
#   is       => 'ro',
#   isa      => AbsFile,
#   required => 1,
#   init_arg => undef,
#   default  => sub {
#     return path( abs_path(__FILE__) )->absolute('/')
#       ->parent->parent->child('./bin/prokaryotic_annotator/vcf-annotator');
#   }
# );
