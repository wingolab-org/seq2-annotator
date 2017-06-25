#!/usr/bin/env perl

use 5.10.0;
use strict;
use warnings;
use lib './lib';
use SeqElastic;
use Getopt::Long;
use DDP;

use Path::Tiny qw/path/;
use Pod::Usage;

my (
  $indexName,   $verbose, $dryRunInsertions, $logDir, $debug, $annotatedFilePath,
  $typeName, $configPath,
);

# usage
GetOptions(
  'n|index_name=s'     => \$indexName,
  't|index_type=s'     => \$typeName,
  'v|verbose'    => \$verbose,
  'd|debug=i'      => \$debug,
  'a|annotated_file_path=s' => \$annotatedFilePath,
  'd|dry_run_insertions|dry|dryRun' => \$dryRunInsertions,
  'l|log_dir=s' => \$logDir,
  'config|config=s' => \$configPath,
);

unless ($indexName && $typeName && $annotatedFilePath) {
  Pod::Usage::pod2usage();
}

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

$year += 1900;
#   # set log file
my $log_name = join '.', 'index', $indexName, "$mday\_$mon\_$year\_$hour\:$min\:$sec", 'log';

my $logPath = path($logDir || "/mnt/annotator_databases/logs/")->child($log_name)->absolute->stringify;

my $app = SeqElastic->new({
  indexName => $indexName,
  indexType  => $typeName,
  verbose => $verbose,
  debug => $debug || 0,
  annotatedFilePath => $annotatedFilePath,
  dryRun => $dryRunInsertions,
  logPath => $logPath,
  config => $configPath,
});

$app->go;
=head1 NAME

snpfile_annotate_mongo_command_line.pl

=head1 DESCRIPTION

This program annotates a snpfile using a binary index of a genome and annotation stored in a mongodb instance. The binary index of the genome is created by `create_mongo_genome.pl`, which creates the indexed genome and stores the annotation in a mongo database.

=head1 VALID_FILES

1. Seqant snp file (typically has .snp extension, but we accept any extension, as long as file is properly formatted (see below)

2. vcf file (typically has .vcf extension, but we accept any extension, as long as file is properly formatted (see below))

=head1 VALID_FORMATS

1. Seqant snp file format: tab delimited, with headers: Fragment Position Reference Alleles Allele_Counts Type  SampleID1 SampleID1 SampleID2 SampleID2 (every other sample header corresponds to the preceeding sample's allele probability, and can be blank) 

2. vcf file format: http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41

=head1 EXAMPLES

  snpfile_annotate_mongo_command_line.pl [-v] [-j] -f my_snpfile -a hg38 

Example output files (non-prokaryotic):

  1. my_snpfile.annotaiton.txt          => line-by-line annotation of sites 
  2. my_snpfile.annotation.json         => JSON data structure of annotations
  3. my_snpfile.annotation.summary.txt  => summary stats for ids in snpfile
  4. my_snpfile.annotation.summary.json => json of summary stats for ids

Example output files (non-prokaryotic):
  N/A : Prints annotation to screen

=head1 AUTHOR

Thomas Wingo <thomas.wingo@emory.edu>
https://wingolab.org

=head1 MAINTAINER

Thomas Wingo <thomas.wingo@emory.edu>

=head1 SEE ALSO

SeqAnt by Amol Shetty
http://seqant.sourceforge.net/

=head1 BUGS

All software has bugs. If you think you have found one please email the 
author.

=head1 COPYRIGHT

Copyright (c) 2014 Thomas S. Wingo (<thomas.wingo@emory.edu>). All rights
reserved.

=head1 LICENSE

This program is free software; you can distribute it and/or modify it under the same terms as GNU GPL 3.0 as published by the Free Software Foundation. You should have received a copy of the GNU Lesser General Public License along with this software; if not, write to: 

  Free Software Foundation, Inc., 
  59 Temple Place, Suite 330
  Boston, MA 02111-1307 USA

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut

# use Carp;
# use Getopt::Long;
# use File::Spec;
# use Pod::Usage;
# use Type::Params qw/ compile /;
# use Types::Standard qw/ :type /;
# use Log::Any::Adapter;
# use YAML::XS qw/ Dump LoadFile /;
# use Try::Tiny;

# use Data::Dump qw/ dump pp /;

# use Seq;

# my (
#   $snpfile,     $yaml_config, $verbose,
#   $help,        $out_file,    $overwrite,
#   $no_skip_chr, $debug,       $snp_file_version
# );
# $snp_file_version = 'snp_2';
# $debug            = 0;

# # TODO: read directly from argument_format.json

# # usage
# GetOptions(
#   'c|config=s'  => \$yaml_config,
#   's|snpfile=s' => \$snpfile,
#   'v|verbose'   => \$verbose,
#   'h|help'      => \$help,
#   'o|out=s'     => \$out_file,
#   'overwrite'   => \$overwrite,
#   'no_skip_chr' => \$no_skip_chr,
#   'd|debug=n'   => \$debug,
#   't|type=s'    => \$snp_file_version,
# );

# if ($help) {
#   Pod::Usage::pod2usage(1);
#   exit;
# }

# unless ( $yaml_config
#   and $snpfile
#   and $out_file
#   and $snp_file_version )
# {
#   Pod::Usage::pod2usage();
# }

# try {
#   # sanity checking
#   if ( -f $out_file && !$overwrite ) {
#     say "ERROR: '$out_file' already exists. Use '--overwrite' switch to over write it.";
#     exit(1);
#   }

#   if ( $debug < 0 || $debug > 2 ) {
#     say "ERROR: debug out of range (0-2); found $debug";
#     exit(1);
#   }

#   unless ( $snp_file_version eq 'snp_1' or $snp_file_version eq 'snp_2' ) {
#     say "ERROR: Expected snp file version to be either snp_1 or snp_2 but found . "
#       . $snp_file_version;
#     exit;
#   }

#   # get absolute path
#   $snpfile     = File::Spec->rel2abs($snpfile);
#   $out_file    = File::Spec->rel2abs($out_file);
#   $yaml_config = File::Spec->rel2abs($yaml_config);
#   say "writing annotation data here: $out_file" if $verbose;

#   # read config file to determine genome name for loging and to check validity of config
#   my $config_href = LoadFile($yaml_config)
#     || die "ERROR: Cannot read YAML file - $yaml_config: $!\n";

#   say pp($config_href) if $debug;

#   # set log file
#   my $log_name = join '.', $out_file, 'annotation', $config_href->{genome_name}, 'log';
#   my $log_file = File::Spec->rel2abs( ".", $log_name );
#   say "writing log file here: $log_file" if $verbose;
#   Log::Any::Adapter->set( 'File', $log_file );

#   # create the annotator
#   my $annotate_instance = Seq->new(
#     {
#       file_type          => $snp_file_version,
#       config_file        => $yaml_config,
#       debug              => $debug,
#       overwrite          => $overwrite,
#       ignore_unknown_chr => ( !$no_skip_chr ),
#       out_file           => $out_file,
#       snpfile            => $snpfile,
#     }
#   );

#   my $href = $annotate_instance->annotate_snpfile;
#   my $fh = IO::File->new( "$out_file.stats.txt", 'w' ) || die "$!";
#   print {$fh} Dump($href);
# }
# catch {
#   say $_;
# }
__END__

=head1 NAME

annotate_snpfile - annotates a snpfile using a given genome assembly specified
in a configuration file

=head1 SYNOPSIS

annotate_snpfile.pl --config <assembly config> --snp <snpfile> --out <file_ext> --type <snp_1, snp_2>

=head1 DESCRIPTION

C<annotate_snpfile.pl> takes a yaml configuration file and snpfile and gives
the annotations for the sites in the snpfile.

=head1 OPTIONS

=over 8

=item B<-s>, B<--snp>

Snp: snpfile

=item B<-r>, B<--type>

Type: version of snpfile: snp_1 or snp_2

=item B<-c>, B<--config>

Config: A YAML genome assembly configuration file that specifies the various
tracks and data associated with the assembly. This is the same file that is also
used by the Seq Package to build the binary genome without any alteration.

=item B<-o>, B<--out>

Output directory: This is the output director.

=item B<--overwrite>

Overwrite: Overwrite the annotation file if it exists.

=item B<--no_skip_chr>

No_Skip_Chr: Try to annotate all chromosomes in the snpfile and die if unable
to do so.


=back

=head1 AUTHOR

Thomas Wingo

=head1 SEE ALSO

Seq Package

=cut
