use 5.10.0;
use strict;
use warnings;

package Seq;

our $VERSION = '0.001';

# ABSTRACT: Annotate a snp file

# TODO: make temp_dir handling more transparent
use Mouse 2;
use Types::Path::Tiny qw/AbsPath AbsFile AbsDir/;

use namespace::autoclean;

use DDP;

# IO::FDPass recommended for MCE::Shared
# https://github.com/marioroy/mce-perl/issues/5
use IO::FDPass;
use MCE::Loop;
use MCE::Shared;

use Seq::InputFile;
use Seq::Output;
use Seq::Headers;

use Seq::DBManager;
use Path::Tiny;
use File::Which qw/which/;
use Carp qw/croak/;
use Scalar::Util qw/looks_like_number/;
use List::MoreUtils qw/first_index/;

use Cpanel::JSON::XS;

extends 'Seq::Base';

# We  add a few of our own annotation attributes
# These will be re-used in the body of the annotation processor below
# Users may configure these
has input_file => (is => 'rw', isa => AbsFile, coerce => 1, required => 1,
  handles  => { inputFilePath => 'stringify' }, writer => 'setInputFile');

# Maximum (signed) size of del allele
has maxDel => (is => 'ro', isa => 'Int', default => -32, writer => 'setMaxDel');

# Defines most of the properties that can be configured at run time
# Needed because there are variations of Seq.pm, ilke SeqFromQuery.pm
# Requires logPath to be provided (currently found in Seq::Base)
with 'Seq::Definition', 'Seq::Role::Validator';

# To initialize Seq::Base with only getters
has '+gettersOnly' => (init_arg => undef, default => 1);

# TODO: further reduce complexity
sub BUILD {
  my $self = shift;

  if($self->maxDel > 0) {
    $self->setMaxDel(-$self->maxDel);
  }

  ########### Create DBManager instance, and instantiate track singletons #########
  # Must come before statistics, which relies on a configured Seq::Tracks
  #Expects DBManager to have been given a database_dir
  $self->{_db} = Seq::DBManager->new();

  # Set the lmdb database to read only, remove locking
  # We MUST make sure everything is written to the database by this point
  # Disable this if need to rebuild one of the meta tracks, for one run
  $self->{_db}->setReadOnly(1);

  # We separate out the reference track getter so that we can check for discordant
  # bases, and pass the true reference base to other getters that may want it (like CADD)
  # Store these references as hashes, to avoid accessor penalty
  $self->{_refTrackGetter} = $self->tracksObj->getRefTrackGetter();
  $self->{_trackGettersExceptReference} = $self->tracksObj->getTrackGettersExceptReference();

  ######### Build the header, and write it as the first line #############
  my $headers = Seq::Headers->new();

  # Seqant has a single pseudo track, that is always present, regardless of whether
  # any tracks exist
  # Note: Field order is required to stay in the follwoing order, because
  # current API allows use of constant as array index:
  # these to be configured:
  # idx 0:  $self->chromField,
  # idx 1: $self->posField,
  # idx 2: $self->typeField,
  # idx 3: $self->discordantField,
  # index 4: $self->altField,
  # index 5: $self->heterozygotesField,
  # index 6: $self->homozygotesField
  $headers->addFeaturesToHeader([
    #index 0
    $self->chromField,
    #index 1
    $self->posField,
    #index 2
    $self->typeField,
    #index 3
    $self->discordantField,
    #index 4
    $self->altField,
    #index 5
    $self->heterozygotesField,
    #index 6
    $self->homozygotesField,
    #index 7
    $self->missingField,
  ], undef, 1);

  $self->{_lastHeaderIdx} = $#{$headers->get()};

  $self->{_trackIdx} = $headers->getParentFeaturesMap();

  ################### Creates the output file handler #################
  # Used in makeAnnotationString
  $self->{_outputter} = Seq::Output->new();

  ################## Make the full output path ######################
  # The output path always respects the $self->output_file_base attribute path;
  $self->{_outPath} = $self->_workingDir->child($self->outputFilesInfo->{annotation});
}

sub annotate {
  my $self = shift;

  $self->log( 'info', 'Checking input file format' );

  my $err;
  my $fh = $self->get_read_fh($self->input_file);
  my $firstLine = <$fh>;

  $err = $self->setLineEndings($firstLine);

  if($err) {
    $self->_errorWithCleanup($err);
    return ($err, undef);
  }

  # TODO: For now we only accept tab separated files
  # We could change this, although comma separation causes may cause with our fields
  # And is slower, since we cannot split on a constant
  # We would also need to take care with properly escaping intra-field commas
  # $err = $self->setDelimiter($firstLine);

  # if($err) {
  #   $self->_errorWithCleanup($err);
  #   return ($err, undef);
  # }

  # Calling in annotate allows us to error early
  ($err, $self->{_chunkSize}) = $self->getChunkSize($self->input_file, $self->max_threads);

  if($err) {
    $self->_errorWithCleanup($err);
    return ($err, undef);
  }

  $self->log('debug', "chunk size is $self->{_chunkSize}");

  #################### Validate the input file ################################
  # Converts the input file if necessary
  ($err, my $fileType) = $self->validateInputFile($self->input_file);

  if($err) {
    $self->_errorWithCleanup($err);
    return ($err, undef);
  }

  # TODO: Handle Sig Int (Ctrl + C) to close db, clean up temp dir
  # local $SIG{INT} = sub {
  #   my $message = shift;
  # };

  if($fileType eq 'snp') {
    $self->log( 'info', 'Beginning annotation' );
    return $self->annotateSnp();
  }

  if($fileType eq 'vcf') {
    $self->log( 'info', 'Beginning annotation' );
    return $self->annotateVcf();
  }

  # TODO: we don't really check for valid vcf, just assume it is
  # So this message is never reached
  $self->_errorWithCleanup("File type isn\'t vcf or snp. Please use one of these files");
  return ("File type isn\'t vcf or snp. Please use one of these files", undef);
}

# TODO: maybe clarify interface; do we really want to return stats and _outputFilesInfo
# or just make those public attributes
sub annotateSnp {
  my $self = shift;

  # File size is available to logProgressAndStatistics
  (my $err, my $inputFileCompressed, my $fh) = $self->get_read_fh($self->inputFilePath);

  if($err) {
    $self->_errorWithCleanup($!);
    return ($!, undef);
  }

  # Copy once to avoid accessor penalty
  my $taint_check_regex = $self->taint_check_regex;

  # Get the header fields we want in the output, and print the header to the output
  my $firstLine = <$fh>;

  my @firstLine;

  chomp $firstLine;

  if ( $firstLine =~ $taint_check_regex ) {
    #Splitting on literal character is much,much faster
    #time perl -e '$string .= "foo \t " for(1..150); for(1..100000) { split('\t', $string) }'
    #vs
    #time perl -e '$string .= "foo \t " for(1..150); for(1..100000) { split("\t", $string) }'
    @firstLine = split $self->delimiter, $1;
  } else {
    $self->_errorWithCleanup("First line of input file has illegal characters: '$firstLine'");
    return ("First line of input file has illegal characters: '$firstLine'", undef);
  }

  my $inputFileProcessor = Seq::InputFile->new();

  $err = $inputFileProcessor->checkInputFileHeader(\@firstLine);

  if($err) {
    $self->_errorWithCleanup($err);
    return ($err, undef);
  }

  ########## Gather the input fields we want to use ################
  # Names and indices of input fields that will be added as the first few output fields
  # Initialized after inputter reads first file line, to account for snp version diffs
  # The first 3 fields are always Fragment, Position, Reference,
  # The next 2 are either Alleles, Type or Type, Alleles, so account for that

  $self->{_alleleFieldIdx} = $inputFileProcessor->alleleFieldIdx;
  $self->{_typeFieldIdx} = $inputFileProcessor->typeFieldIdx;

  # If user specified a temp output path, use that
  my $outFh = $self->get_write_fh( $self->{_outPath} );

  # Stats may or may not be provided
  my $statsFh;

  # Write the header
  my $headers = Seq::Headers->new();
  my $outputHeader = $headers->getString();
  say $outFh $outputHeader;

  # TODO: error handling if fh fails to open
  if($self->run_statistics) {
    my $args = $self->_statisticsRunner->getStatsArguments();
    open($statsFh, "|-", $args);

    say $statsFh $outputHeader;
  }

  ############# Set the sample ids ###############
  # sample name by the genotype idx (the sample names are stored in the sam column as the genotype)
  $self->{_inputHeader} = \@firstLine;
  $self->{_sampleGenosIdx} = $inputFileProcessor->getSampleNamesGenos(\@firstLine);

  if(defined $self->{_sampleGenosIdx}) {
    $self->{_inputHeader} = $self->_normalizeSampleNames(\@firstLine, $self->{_sampleGenosIdx});
  }

  my $abortErr;
  MCE::Loop::init {
    max_workers => $self->max_threads || 8, use_slurpio => 1,
    # Small chunk_size progress messages take large portion of recipient's run time
    # Causing jobs to appear to take much longer to complete.
    # So in the end, we optimize for files at least ~65k lines in size, which
    # is of course a very small file size
    #parallel_io => 1,
    chunk_size => $self->{_chunkSize} . 'K',
    gather => $self->makeLogProgressAndPrint(\$abortErr, $outFh, $statsFh),
  };

  # Avoid hash lookups when possible, those are slow too
  my $typeFieldIdx = $self->{_typeFieldIdx};

  # Store a reference; Moose/Mouse accessors aren't terribly fast
  # And this is used millions of times
  my %wantedChromosomes = %{ $self->{_refTrackGetter}->chromosomes };
  my $refTrackGetter = $self->{_refTrackGetter};

  # If the file isn't compressed, MCE::Loop is faster when given a string
  # instead of a file handle
  # https://github.com/marioroy/mce-perl/issues/5
  if(!$inputFileCompressed) {
    close $fh;
    $fh = $self->inputFilePath;
  }

  # MCE::Mutex must come after any close operations on piped file handdles
  # or MCE::Shared::stop() must be called before they are closed
  # https://github.com/marioroy/mce-perl/issues/5
  # We want to check whether the program has any errors. An easy way is to
  # Return any error within the mce_loop_f
  # Doesn't work nicely if you need to return a scalar value (no export)
  # and need to call MCE::Shared::stop() to exit 

  my $m1 = MCE::Mutex->new;
  tie my $readFirstLine, 'MCE::Shared', 0;

  mce_loop_f {
    #my ($mce, $slurp_ref, $chunk_id) = @_;
    #    $_[0], $_[1],     $_[2]
    my @lines;

    # Reads: open my $MEM_FH, '<', $slurp_ref; binmode $MEM_FH, ':raw';
    open my $MEM_FH, '<', $_[1]; binmode $MEM_FH, ':raw';

    # If the file isn't compressed, we pass the file path to 
    # MCE, because it is faster that way...
    # Howver, when we do that, we need to skip the header
    # Reads:                                        && chunk_id == 1) {
    if(!$inputFileCompressed && $readFirstLine == 0 && $_[2] == 1) {
      #skip first line
      my $firstLine = <$MEM_FH>;
      $m1->synchronize(sub{ $readFirstLine = 1 });
    }

    my $skipCount = 0;
    # Note: We don't sanitize anything other than the header line, because
    # 1) We don't evaluate any user input
    # 2) We don't trust most user input on a per value basis, because:
    #  a) the chrom field must match an expected value
    #  b) the pos field must be a position found in our database
    #  c) the ref field is not stored
    #  d) we explicitily only allow certain Type field values
    # TODO: Currently someone "could" add unsafe values to the alt field
    while ( my $line = $MEM_FH->getline() ) {
      chomp $line;
      #Splitting on literal character is much,much faster
      #time perl -e '$string .= "foo \t " for(1..150); for(1..100000) { split('\t', $string) }'
      #vs
      #time perl -e '$string .= "foo \t " for(1..150); for(1..100000) { split("\t", $string) }'
      #$isTab ? split '\t', $line : split $delimiter, $line is even slower than variable split
      my @fields = split '\t', $line;

      if ( !$wantedChromosomes{$fields[0]} ) {
        # We expect snp fields to be be ucsc chr format, but if not, correct that
        $fields[0] = $refTrackGetter->normalizedWantedChr($fields[0]);

        if(!$fields[0]) {
          $skipCount++;
          next;
        }
      }

      # Explicitly only allow some user inputs. Our benchmarks suggest this is fastest
      # When variants are primarily SNP, which is what we expect
      if(!defined $typeFieldIdx || $fields[$typeFieldIdx] eq "SNP" ||
      $fields[$typeFieldIdx] eq "INS" || $fields[$typeFieldIdx] eq "DEL" ||
      $fields[$typeFieldIdx] eq "MULTIALLELIC" || $fields[$typeFieldIdx] eq "DENOVO_SNP" ||
      $fields[$typeFieldIdx] eq "DENOVO_INS" || $fields[$typeFieldIdx] eq "DENOVO_DEL" ||
      $fields[$typeFieldIdx] eq "DENOVO_MULTIALLELIC") {
        push @lines, \@fields;
        next;
      }

      $skipCount++;
    }

    close $MEM_FH;

    my ($err, $outAref);

    if(@lines) {
      #TODO: implement better error handling
      # ($err, $outString) = $self->makeAnnotationString(\@lines);
      ($err, $outAref) = $self->addTrackData(\@lines);

      if($err) {
        #           $annotated, $skippd, $error
        MCE->gather(0, $skipCount + @lines, $err);
      } else {
        MCE->gather(scalar @$outAref, @lines - @$outAref, undef, $self->{_outputter}->makeOutputString($outAref));
      }

      undef @lines;
    } else {
      MCE->gather(0, $skipCount);
    }
  } $fh;

  MCE::Loop::finish();

  # This removes the content of $abortErr
  # https://metacpan.org/pod/MCE::Shared
  # Needed to exit, and close piped file handles
  MCE::Shared::stop();

  # Unfortunately, MCE::Shared::stop() removes the value of $abortErr
  # according to documentation, and I did not see mention of a way
  # to copy the data from a scalar, and don't want to use a hash for this alone
  # So, using a scalar ref to abortErr in the gather function.
  if($abortErr) {
    say "Aborted job";

    # Database & tx need to be closed
    $self->{_db}->cleanUp();

    return ('Job aborted due to error', undef);
  }

  ################ Finished writing file. If statistics, print those ##########
  my $statsHref;
  if($statsFh) {
    # Force the stats program to write its outputs
    close $statsFh;
  }

  close $outFh;
  system('sync');

  ################ Compress if wanted ##########
  $err = $self->_moveFilesToOutputDir();

  # If we have an error moving the output files, we should still return all data
  # that we can
  if($err) {
    $self->log('error', $err);
  }

  $self->{_db}->cleanUp();

  return ($err || undef, $self->outputFilesInfo);
}

sub annotateVcf {
  #Inspired by T.S Wingo: https://github.com/wingolab-org/GenPro/blob/master/bin/vcfToSnp
  my $self = shift;

  my $errPath = $self->_workingDir->child($self->input_file->basename . '.vcf-log.log');
  my $inPath = $self->inputFilePath;
  my $echoProg = $self->isCompressedSingle($inPath) ? $self->gzip . ' -d -c' : 'cat';
  my $delim = $self->{_outputter}->delimiters->emptyFieldChar;
  # Retruns chr, pos, homozygotes, heterozygotes, alt, ref in that order, tab delim
  open(my $fh, '-|', "$echoProg $inPath | seq-vcf --emptyField $delim 2> $errPath");

  # If user specified a temp output path, use that
  my $outFh = $self->get_write_fh( $self->{_outPath} );

  ########################## Write the header ##################################
  my $headers = Seq::Headers->new();
  my $outputHeader = $headers->getString();

  say $outFh $outputHeader;

  ########################## Tell stats program about our annotation ##############
  # TODO: error handling if fh fails to open
  my $statsFh;
  if($self->run_statistics) {
    my $args = $self->_statisticsRunner->getStatsArguments();
    open($statsFh, "|-", $args);

    say $statsFh $outputHeader;
  }

  my $abortErr;
  MCE::Loop::init {
    max_workers => $self->max_threads || 8, use_slurpio => 1,
    # seq-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => $self->{_chunkSize} > 4192 ? "4192K" : $self->{_chunkSize}. "K",
    gather => $self->makeLogProgressAndPrint(\$abortErr, $outFh, $statsFh),
  };

  my $trackIndices = $self->{_trackIdx};
  my $refTrackIdx = $self->{_trackIdx}{$self->{_refTrackGetter}->name};
  my @trackGettersExceptReference = @{$self->{_trackGettersExceptReference}};
  my %wantedChromosomes = %{ $self->{_refTrackGetter}->chromosomes };
  my $maxDel = $self->maxDel;

  my $err = $self->setLineEndings("\n");

  if($err) {
    $self->_errorWithCleanup($err);
    return ($err, undef);
  }

  mce_loop_f {
    #my ($mce, $slurp_ref, $chunk_id) = @_;
    #    $_[0], $_[1],     $_[2]
    #open my $MEM_FH, '<', $slurp_ref; binmode $MEM_FH, ':raw';
    open my $MEM_FH, '<', $_[1]; binmode $MEM_FH, ':raw';

    my $total = 0;

    my @indelDbData;
    my @indelRef;
    my $chr;
    my $inputRef;
    my @lines;
    my $dataFromDbAref;
    # Each line is expected to be
    # chrom \t pos \t type \t inputRef \t alt \t hets \t homozygotes \n
    # the chrom is always in ucsc form, chr (the golang program guarantees it)
    while (my $line = $MEM_FH->getline()) {
      chomp $line;

      my @fields = split '\t', $line;
      $total++;

      if(!$wantedChromosomes{$fields[0]}) {
        next;
      }

      $dataFromDbAref = $self->{_db}->dbReadOne($fields[0], $fields[1] - 1, 1);

      if(length($fields[4]) > 1) {
        # INS or DEL
        if(looks_like_number($fields[4])) {
          # We ignore -1 alleles, treat them just like SNPs
          if($fields[4] < -1)  {
            # Grab everything from + 1 the already fetched position to the $pos + number of deleted bases - 1
            # Note that position_1_based - (negativeDelLength + 2) == position_0_based + (delLength - 1)
            if($fields[4] < $maxDel) {
              @indelDbData = ($fields[1] .. $fields[1] - ($maxDel + 2));
              # $self->log('info', "$fields[0]:$fields[1]: long deletion. Annotating up to $maxDel");
            } else {
              @indelDbData = ($fields[1] .. $fields[1] - ($fields[4] + 2));
            }

            #last argument: skip commit
            $self->{_db}->dbRead($fields[0], \@indelDbData, 1);

            #Note that the first position keeps the same $inputRef
            #This means in the (rare) discordant multiallelic situation, the reference
            #Will be identical between the SNP and DEL alleles
            #faster than perl-style loop (much faster than c-style)
            @indelRef = ($fields[3], map { $self->{_refTrackGetter}->get($_) } @indelDbData);

            #Add the db data that we already have for this position
            unshift @indelDbData, $dataFromDbAref;
          }
        } else {
          #It's an insertion, we always read + 1 to the position being annotated
          # which itself is + 1 from the db position, so we read  $out[1][0][0] to get the + 1 base
          # Read without committing by using 1 as last argument
          @indelDbData = ($dataFromDbAref, $self->{_db}->dbReadOne($fields[0], $fields[1], 1));

          #Note that the first position keeps the same $inputRef
          #This means in the (rare) discordant multiallelic situation, the reference
          #Will be identical between the SNP and DEL alleles
          @indelRef = ($fields[3], $self->{_refTrackGetter}->get($indelDbData[1]));
        }
      }

      if(@indelDbData) {
        ############### Gather all track data (besides reference) #################
        for my $posIdx (0 .. $#indelDbData) {
          for my $track (@trackGettersExceptReference) {
            $fields[$trackIndices->{$track->name}] //= [];

            $track->get($indelDbData[$posIdx], $fields[0], $indelRef[$posIdx], $fields[4], 0, $posIdx, $fields[$trackIndices->{$track->name}]);
          }

          $fields[$refTrackIdx][0][$posIdx] = $indelRef[$posIdx];
        }

        # If we have multiple indel alleles at one position, need to clear stored values
        @indelDbData = ();
        @indelRef = ();
      } else {
        for my $track (@trackGettersExceptReference) {
          $fields[$trackIndices->{$track->name}] //= [];
          $track->get($dataFromDbAref, $fields[0], $fields[3], $fields[4], 0, 0, $fields[$trackIndices->{$track->name}])
        }

        $fields[$refTrackIdx][0][0] = $self->{_refTrackGetter}->get($dataFromDbAref);
      }

       # 3 holds the input reference, we'll replace this with the discordant status
      $fields[3] = $self->{_refTrackGetter}->get($dataFromDbAref) ne $fields[3] ? 1 : 0;
      push @lines, \@fields;
    }

    if(@lines) {
      MCE->gather(scalar @lines, $total - @lines, undef, $self->{_outputter}->makeOutputString(\@lines));
    } else {
      MCE->gather(0, $total);
    }
  } $fh;

  MCE::Loop::finish();

  # This removes the content of $abortErr
  # https://metacpan.org/pod/MCE::Shared
  # Needed to exit, and close piped file handles
  MCE::Shared::stop();

  # Unfortunately, MCE::Shared::stop() removes the value of $abortErr
  # according to documentation, and I did not see mention of a way
  # to copy the data from a scalar, and don't want to use a hash for this alone
  # So, using a scalar ref to abortErr in the gather function.
  if($abortErr) {
    say "Aborted job";

    # Database & tx need to be closed
    $self->{_db}->cleanUp();

    return ('Job aborted due to error', undef);
  }

  ################ Finished writing file. If statistics, print those ##########
  # Sync to ensure all files written
  close $outFh;

  if($statsFh) {
    close $statsFh;
  }

  system('sync');

  $err = $self->_moveFilesToOutputDir();

  # If we have an error moving the output files, we should still return all data
  # that we can
  if($err) {
    $self->log('error', $err);
  }

  $self->{_db}->cleanUp();

  return ($err || undef, $self->outputFilesInfo);
}

sub makeLogProgressAndPrint {
  my ($self, $abortErrRef, $outFh, $statsFh) = @_;

  my $totalAnnotated = 0;
  my $totalSkipped = 0;

  my $publish = $self->hasPublisher;

  return sub {
    #my $annotatedCount, $skipCount, $err, $outputLines = @_;
    ##    $_[0],          $_[1],     $_[2], $_[3]
    if($_[2]) {
      $$abortErrRef = $_[2];
      return;
    }

    if($publish) {
      $totalAnnotated += $_[0];
      $totalSkipped += $_[1];
      $self->publishProgress($totalAnnotated, $totalSkipped);
    }

    if($_[3]) {
      if($statsFh) {
        print $statsFh $_[3];
      }

      print $outFh $_[3];
    }
  }
}

###Private genotypes: used to decide whether sample is het, hom, or compound###
my %hets = (K => 1,M => 1,R => 1,S => 1,W => 1,Y => 1,E => 1,H => 1);
my %homs = (A => 1,C => 1,G => 1,T => 1,D => 1,I => 1);
my %iupacArray = (A => ['A'], C => ['C'], G => ['G'], T => ['T'], R => ['A','G'],
  Y => ['C','T'],S => ['G','C'],W => ['A','T'],K => ['G','T'],M => ['A','C']);
my %indels = (E => '-',H => '+',D => '-',I => '+');

# The main annotation function
# @params <ArrayRef> $inputAref : the input lines
# @returns <ArrayRef> : the annotated lines
# Expects the first 3 fields in the input file to be chrom, pos, ref
sub addTrackData {
  my ($self, $inputAref) = @_;
  # Cache $alleles; for large jobs, or long-running environments
  # This will save us millions of split's and assignments
  # Careful, this could be a place for a subtle bug in a multi-process
  # long-running environment
  # ONLY WORKS IF WE ARE COMPARING INPUT'S Ref & Alleles columns
  state $cached = {};

  # We store data for each indel
  my (@indelDbData, @indelRef);

  # We accumulate all of our results, using alleleNumber
  # Each track get method accumulates values on its own, using this
  # to keep track of whether or not an array represents multiple alleles
  # or one allele's values for a particular field
  # We phase heterozygotes and homozygotes on allele
  my ($alleleIdx, $inputRef, $alleles, $strAlleles, $gtIdx, $geno);
  my ($pos, $chr, @finalOut, $dataFromDbAref);

  my $refTrackIdx = $self->{_trackIdx}{$self->{_refTrackGetter}->name};
  my $alleleFieldIdx = $self->{_alleleFieldIdx};
  my $sampleNames = $self->{_inputHeader};
  my $trackIndices = $self->{_trackIdx};

  my $maxDel = $self->maxDel;
  POSITION_LOOP: for my $inputRow (@$inputAref) {
    # The first 3 fields are expected to gives chrom, pos, ref in that order
    $chr = $inputRow->[0]; #$inputAref->[$i][$self->{_chrFieldIdx}];
    $pos = $inputRow->[1]; #$inputAref->[$i][$self->{_positionFieldIdx}];
    $inputRef = $inputRow->[2]; #$inputAref->[$i][$self->{_referenceFieldIdx}];

    # Read without committing (last argument 1)
    $dataFromDbAref = $self->{_db}->dbReadOne($chr, $pos - 1, 1);

    if(!defined $dataFromDbAref) {
      return ($self->_errorWithCleanup("Wrong assembly? $chr\: $pos not found."), 0, undef);
    }

    my @out;
    # Set array size, initializing every value to undef
    $#out = $self->{_lastHeaderIdx};

    ############# Store chr, position, alleles, type, and discordant status ###############
    $out[0][0][0] = $chr;

    $out[1][0][0] = $pos;

    if(defined $self->{_typeFieldIdx}) {
      $out[2][0][0] = $inputRow->[$self->{_typeFieldIdx}];
    }

    # Record whether discordant; boolean returns 1 if true
    # https://ideone.com/9POQD9
    $out[3][0][0] = $self->{_refTrackGetter}->get($dataFromDbAref) ne $inputRef || 0;

    ############### Get the minor alleles, cached to avoid re-work ###############
    # Calculate the minor alleles from the user's reference
    # It's kind of hard to know exactly what to do in discordant cases
    if(!defined $cached->{$inputRef}{$inputRow->[$alleleFieldIdx]}) {
      #https://ideone.com/cb3yvi
      if($inputRow->[$alleleFieldIdx] !~ /^[\+\-\,ACTG0-9]+$/) {
        $cached->{$inputRef}{$inputRow->[$alleleFieldIdx]} = 0;
        next POSITION_LOOP;
      }

      my @alleles = grep { $_ ne $inputRef } split(',', $inputRow->[$alleleFieldIdx]);

      if(@alleles) {
        # Single base, used for fast lookup in sample loop. Stores a + or - for indels
        $cached->{$inputRef}{$inputRow->[$alleleFieldIdx]} =
          @alleles > 1 ? [\@alleles, join('', map{length($_) > 1 ? substr($_, 0, 1) : $_ } @alleles)]
                       : [$alleles[0], length($alleles[0]) > 1 ? substr($alleles[0], 0, 1) : $alleles[0]]
      } else {
        $cached->{$inputRef}{$inputRow->[$alleleFieldIdx]} = 0;
        next POSITION_LOOP;
      }
    } elsif($cached->{$inputRef}{$inputRow->[$alleleFieldIdx]} == 0) {
      next POSITION_LOOP;
    }

    ($alleles, $strAlleles) = @{ $cached->{$inputRef}{$inputRow->[$alleleFieldIdx]} };
    ############ Store homozygotes, heterozygotes, compoundHeterozygotes ########
    # Homozygotes are index 4, heterozygotes 5

    #TODO: Inline -C. Also, we do not need the alleleIdx == -1 check if 
    #we always have the relationship that those genotypes with confidence < .95
    #are excluded from the list of Alleles
    #This may not be true for some vcf files we choose to annotate

    if(defined $self->{_sampleGenosIdx}) {
      SAMPLE_LOOP: for my $idx (@{$self->{_sampleGenosIdx}}) {
        #Does the sample genotype equal our assembly's reference?
        #Skip assignment for the 99.9% of cases that are expected to be ref
        if($inputRow->[$idx] eq $inputRef) {
          next SAMPLE_LOOP;
        }

        $geno = $inputRow->[$idx];

        #Bad samples
        if($geno eq 'N') {
          # In case of multiallelics, we will store this sample as missing
          # for each possible allele, since it is impossible to tell to which the allele applies
          if(ref $alleles) {
            for my $altIdx (0 .. $#$alleles) {
              push @{$out[7][$altIdx][0]}, $sampleNames->[$idx];
            }
          } else {
            push @{$out[7][0][0]}, $sampleNames->[$idx];
          }

          next SAMPLE_LOOP;
        }

        if ($hets{$geno}) {
          # Is this a bi-allelic sample? if so, call that homozygous
          # Indel hets are never bi-allelic, limitation of PECaller merge script
          if($indels{$geno}) {
            # Heterozygote is column 5
            $gtIdx = index($strAlleles, $indels{$geno});

            if($gtIdx == -1) {
              $self->log('warn', "Genotype $geno doesn't have a corresponding allele @ $chr:$pos. Skipping");
              next POSITION_LOOP;
            }

            push @{$out[5][$gtIdx][0]}, $sampleNames->[$idx];
          } else {
            #There can be bi-allelic SNPs, where both calls in a het are non-reference
            for my $genoAllele ( @{$iupacArray{$geno}} ) {
              if($genoAllele ne $inputRef) {
                # Heterozygote is column 5
                $gtIdx = index($strAlleles, $genoAllele);

                if($gtIdx == -1) {
                  $self->log('warn', "Genotype $geno doesn't have a corresponding allele @ $chr:$pos. Skipping");
                  next POSITION_LOOP;
                }

                push @{$out[5][$gtIdx][0]}, $sampleNames->[$idx];
              }
            }
          }

          next SAMPLE_LOOP;
        }

        # Check if the sample looks like a homozygote
        if($homs{$geno}) {
          $gtIdx = $indels{$geno} ? index($strAlleles, $indels{$geno}) : index($strAlleles, $geno);

          if($gtIdx == -1) {
            $self->log('warn', "Genotype $geno doesn't have a corresponding allele @ $chr:$pos. Skipping");
            next POSITION_LOOP;
          }

          push @{$out[6][$gtIdx][0]}, $sampleNames->[$idx];

          next SAMPLE_LOOP;
        }

        $self->log( 'warn', "$sampleNames->[$idx] wasn't homozygous or heterozygote: $chr\:$pos: $geno" );
      }

      if(!$out[5] && !$out[6]) {
        # We don't have any data, skip this posiiton
        next POSITION_LOOP;
      }
    }

    $alleleIdx = 0;
    for my $allele (ref $alleles ? @$alleles : $alleles) {
      # The minorAlleles column
      $out[4][$alleleIdx][0] = $allele;

      # In case of multiallelics, if we don't have homs, hets, or missing samples, output empty field value for each allele
      $out[5][$alleleIdx] //= [undef];
      $out[6][$alleleIdx] //= [undef];
      $out[7][$alleleIdx] //= [undef];

      if(length($allele) > 1) {
        # It's a deletion
        if( looks_like_number($allele) ) {
          # If the allele is == -1, it's a single base deletion, treat like a snp
          # with a weird genotype
          if($allele < -1)  {
            # Grab everything from + 1 the already fetched position to the $pos + number of deleted bases - 1
            # Note that position_1_based - (negativeDelLength + 2) == position_0_based + (delLength - 1)
            if($pos < $maxDel) {
              @indelDbData = ($pos .. $pos - ($maxDel + 2));
              # $self->log('info', "$chr:$pos: long deletion. Annotating up to $maxDel");
            } else {
              @indelDbData = ($pos .. $pos - ($allele + 2));
            }

            #last argument: skip commit
            $self->{_db}->dbRead($chr, \@indelDbData, 1);

            #Note that the first position keeps the same $inputRef
            #This means in the (rare) discordant multiallelic situation, the reference
            #Will be identical between the SNP and DEL alleles
            #faster than perl-style loop (much faster than c-style)
            @indelRef = ($inputRef, map { $self->{_refTrackGetter}->get($_) } @indelDbData);

            #Add the db data that we already have for this position
            unshift @indelDbData, $dataFromDbAref;
          }
        } else {
          #It's an insertion, we always read + 1 to the position being annotated
          # which itself is + 1 from the db position, so we read  $out[1][0][0] to get the + 1 base
          # Read without committing by using 1 as last argument
          @indelDbData = ($dataFromDbAref, $self->{_db}->dbReadOne($chr, $pos, 1));

          #Note that the first position keeps the same $inputRef
          #This means in the (rare) discordant multiallelic situation, the reference
          #Will be identical between the SNP and DEL alleles
          @indelRef =  ( $inputRef, $self->{_refTrackGetter}->get($indelDbData[1]) );
        }
      }

       if(@indelDbData) {
        ############### Gather all track data (besides reference) #################
        for my $posIdx (0 .. $#indelDbData) {
          for my $track (@{ $self->{_trackGettersExceptReference} }) {
            $out[$trackIndices->{$track->name}] //= [];

            $track->get($indelDbData[$posIdx], $chr, $indelRef[$posIdx], $allele,
              $alleleIdx, $posIdx, $out[$trackIndices->{$track->name}]);
          }

          $out[$refTrackIdx][$alleleIdx][$posIdx] = $indelRef[$posIdx];
        }

        # If we have multiple indel alleles at one position, need to clear stored values
        @indelDbData = ();
        @indelRef = ();
      } else {
        for my $track (@{ $self->{_trackGettersExceptReference} }) {
          $out[$trackIndices->{$track->name}] //= [];

          $track->get($dataFromDbAref, $chr, $inputRef, $allele,
            $alleleIdx, 0, $out[$trackIndices->{$track->name}])
        }

        $out[$refTrackIdx][$alleleIdx][0] = $inputRef;
      }

      $alleleIdx++;
    }

    #If we get to this point, we have data
    push @finalOut, \@out;
  }

  # 0 status indicates success
  return (0, \@finalOut);
}

# This function is expected to run within a fork, so itself it doesn't clean up 
# files
# TODO: better error handling
sub _errorWithCleanup {
  my ($self, $msg) = @_;

  $self->log('error', $msg);

  $self->{_db}->cleanUp();

  return $msg;
}

__PACKAGE__->meta->make_immutable;

1;
