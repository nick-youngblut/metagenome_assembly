#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw/remove_tree/;
use IPC::Open2;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $dump_connect);
my $max_size = 1000000;
my $min_part_size = 10;
my $status_int = 1000000;
my $inflation = 2;				# inflation param default
my $threads = 1;
GetOptions(
	   "max-size=i" => \$max_size,						# max group size
	   "min-partition-size=i" => \$min_part_size,		# min partition size worth keeping
	   "status=i" => \$status_int,						# mcl inflation param
	   "threads=i" => \$threads,						# how often to provide a status update
	   "inflation=f" => \$inflation,					# number of threads for mcl
	   "dump" => \$dump_connect,						# dumping connected graph (for running mcl outside of script)
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a *part file\n" if ! $ARGV[0];

### MAIN
my $basename = sort_pairs($ARGV[0]);										# parsing out pairs
my $connect_r = partition_connection_hash($basename);						# making a parition graph for loading into mcl
dump_connect_hash($connect_r, $basename) if $dump_connect;
my $mcl_out = call_mcl($connect_r, $basename, $inflation, $threads);		# clustering w/ mcl
#my $clusters_r = load_mcl_out($mcl_out);
#my $groups = parts_in_groups($clusters_r, $max_size, $min_part_size)		# placing paritions into groups
#write_groups($clusters_r, $basename);

### writing out groups ###
#my $parts_r = count_parts($basename);
#my $groups_r = parts_in_groups($parts_r, $max_size, $min_part_size);
#write_groups($groups_r, $basename);

### Subroutines
sub dump_connect_hash{
# writting out connected hash to file for running mcl outside of script #
	my ($connect_r, $basename) = @_;

	# making an output file name #
	my $dump_out = $basename . "_abc.txt";
	
	# status #
	print STDERR " Dumping partition connection hash: '$dump_out'\n";
	
	open OUT, ">$dump_out" or die $!;
	
	# running mcl #
	foreach my $c1 (keys %$connect_r){
		foreach my $c2 (keys %{$$connect_r{$c1}}){
			print OUT join(" ", $c1, $c2, $$connect_r{$c1}{$c2}), "\n";
			}
		}	
	close OUT;
	
	print STDERR " Connection hash dumped. Exiting\n";
	exit;
	}

sub parts_in_groups{
# placing partitions in groups #
	my ($parts_r, $max_size, $min_part_size) = @_;

 	# status #
	print STDERR " Grouping partitions\n";

	# number of groups depends on $max_size
	my %groups;
	my $group_cnt;			# tracking total
	foreach my $part (keys %$parts_r){
		# filtering low abundant reads #
		next if $$parts_r{$part} < $min_part_size;
		
		# loading groups #
		$group_cnt += $$parts_r{$part};
		$groups{$part} = int($group_cnt / $max_size);		# part -> group_number
		}
	
	return \%groups;
	}


sub load_mcl_out{
# loading mcl output into a hash #
	my ($mcl_out) = @_;
	
	# status #
	print STDERR " Loading mcl output\n";
	
	open IN, $mcl_out or die $!;
	my %clusters;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		map{ $clusters{$_} = $. } @line;			# all partitions in cluster pointing to clsuter ID
		}
	close IN;
	
		#print Dumper %clusters; exit;
	return \%clusters;
	}

sub call_mcl{
# calling mcl for clustering #
	my ($connect_r, $basename, $inflation, $threads) = @_;
	
	# status #
	print STDERR " Calling mcl\n";
	
	# making an output file name #
	my $mcl_out = $basename . "_mcl.txt";
	
	my $cmd = "mcl - --abc -I $inflation -te $threads -o $mcl_out";
	print STDERR $cmd, "\n";
	open PIPE, " | $cmd" or die $!;
	
	# running mcl #
	foreach my $c1 (keys %$connect_r){
		foreach my $c2 (keys %{$$connect_r{$c1}}){
			print PIPE join(" ", $c1, $c2, $$connect_r{$c1}{$c2}), "\n";
			}
		}
	close PIPE;

	return $mcl_out;
	}

sub partition_connection_hash{
# merging partitions based on paired-end reads #
	my ($basename) = @_;
	
	# status #
	print STDERR " Making connectedness hash\n";
	
	# reading pair file and making pair hash #
	open PAIR, "$basename-pair.fna" or die $!;

	# making connectedness hash #
	#my $max_connections = 0;		# for normalizing connections
	my %connect;
	while(<PAIR>){
		chomp;
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR " Read pairs processed: ", ($. - 1) / 2, "\n";
			}
	
		# loading lines #
		my @line1 = split /\t/;
		my $line2 = <PAIR>;
		my @line3 = split /\t|\n/, <PAIR>;
		my $line4 = <PAIR>;
		
		# adding to merge hash #
		$connect{$line1[1]}{$line3[1]}++ if $line1[1] ne $line3[1];
			
		}
	close PAIR;

		#print Dumper %connect; exit;

	return \%connect;
	}

sub sort_pairs{
# parsing out existing paired-end reads #
	my ($infile) = @_;
	
	# status #
	print STDERR " Splitting reads into pairs and singletons\n";
	
	# I/O #
	(my $basename = $infile) =~ s/\.[^\.]+$|$//;
	open IN, $infile or die $!;
	open PAIR, ">$basename-pair.fna" or die $!;
	open SING, ">$basename-single.fna" or die $!;

	my %pairs;
	while(<IN>){
		
		#last if $. > 10000;
		
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR " Reads processed: ", ($. - 1) / 2, "\n";
			}
		
		# load lines #
		my @line = split /\//;		# header, pair, partition
		die " ERROR: read names must be in old illumina format (>NAME/1 or >NAME/2)\n"
			if scalar @line != 2;
		my $nline = <IN>;

		# checking for pairs; writing files #
		if( exists $pairs{$line[0]} ){
			print PAIR join("", $pairs{$line[0]}, join("/", @line), $nline);
			delete $pairs{$line[0]};
			}
		else{
			$pairs{$line[0]} = join("", join("/", @line), $nline);
			}

		}
	close IN;
	
	# writing out all singletons #
	print STDERR " Writing out all singletons\n";
	print SING join("", values %pairs);
	
	close PAIR;
	close SING;
	
	return $basename;
	}



__END__

=pod

=head1 NAME

extract-partitions-paired.pl -- extract-partitions.py, but retaining paired-end reads

=head1 SYNOPSIS

extract-partitions-paired.pl [options] file.part

=head2 options

=over

=item --max-size 

Max group size. [1000000]

=item --min-partition-size	

Minimum partition size (number of reads) to retain. [10]

=item -s 

Status output interval (number of reads). [1000000]

=item -h

This help message

=back

=head2 For more information:

perldoc extract-partitions-paired.pl

=head1 DESCRIPTION

The script is basically the same as extract-partitions.py, but it retains paired-end
reads in 1 group. 

The reasoning is that paired-end reads should be from the same molecule and thus should
come from the same organism.

The reads are parsed by group and by singleton or paired status, so that velvet or idba_ud
can be used on the reads.


=head1 EXAMPLES

$ extract-partitions-paired.pl iowa-corn-50m.fa.gz.part
# output to ./iowa-corn-50m.fa.gz_groups/

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

