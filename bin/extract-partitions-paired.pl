#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw/remove_tree/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $max_size = 1000000;
my $min_part_size = 5;
GetOptions(
	   "max-size=i" => \$max_size,						# max group size
	   "min-partition-size=i" => \$min_part_size,		# min partition size worth keeping
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a *part file\n" if ! $ARGV[0];


### MAIN
find_pairs_and_write($ARGV[0]);

### Subroutines
sub find_pairs_and_write{
	my $infile = shift;
	
	# I/O #
	(my $basename = $infile) =~ s/\.[^\.]+$|$//;
	
	remove_tree(".$basename") if -d ".$basename";
	mkdir ".$basename"; 
	
	open IN, $infile or die $!;
	
	# running through file #
	my %seq_names;
	my %parts;
	while(<IN>){
		my @line = split /\t|\//;		# header, pair, partition
		
		if(! exists($parts{$line[2]})){		# if need to make a file
			open $parts{$line[2]}, ">.$basename/$line[2].part" or die $!;
			
			# writing read #
			print {$parts{$line[2]}} $_;
			my $nline = <IN>;
			print {$parts{$line[2]}} $nline;
			}
		else{
			print {$parts{$line[2]}} $_;
			my $nline = <IN>;
			print {$parts{$line[2]}} $nline;
			}
		
		# checking for other pair #
		#if(exists $seq_names{$line[1]} ){		# if paired file found 		
		#	}
		
		
		}
	
	# closing #
	close IN;
	map{ close $parts{$_} } keys %parts; 
	}


__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

