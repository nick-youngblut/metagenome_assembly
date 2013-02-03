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
my $basename = sort_pairs($ARGV[0]);													# parsing out pairs
my ($connect_r, $abunds_r) = partition_connection_hash($basename, $min_part_size);		# making a parition graph for loading into mcl
dump_connect_hash($connect_r, $basename) if $dump_connect;
my $mcl_out = call_mcl($connect_r, $basename, $inflation, $threads);					# clustering w/ mcl
my ($clusters_r, $cluster_abunds_r) = load_mcl_out($mcl_out, $abunds_r);
my ($groups_r, $group_abunds_r) = parts_in_groups($clusters_r, $cluster_abunds_r, $max_size);				# placing paritions into groups
write_groups($groups_r, $group_abunds_r, $basename);

### writing out groups ###
#my $parts_r = count_parts($basename);
#my $groups_r = parts_in_groups($parts_r, $max_size, $min_part_size);
#write_groups($groups_r, $basename);

### Subroutines
sub write_groups{
# writing out groups (group => clusters => partitions) #
	my ($groups_r, $group_abunds_r, $basename) = @_;
	
	# I/O #
	## read paired-end reads ##
	open PAIR, "$basename-pair.fna" or die $!;

	## making directory for group files ##
	remove_tree("$basename\_groups") if -d "$basename\_groups";
	mkdir "$basename\_groups" or die $!;
	my %outfh;
	map{ open $outfh{$_}, ">$basename\_groups/$basename\_g$_.fna" or die $! } keys %$group_abunds_r;
	
	# writing to group files #
	my %stats;
	while(<PAIR>){

		# loading lines #
		my @pairs;
		push( @pairs, split /\t|\n/ );
		for my $i (0..2){
			my $line = <PAIR>;
			push( @pairs, split /\t|\n/, $line);
			}
			#print Dumper @pairs; exit;
		
		if(exists $$groups_r{$pairs[1]} && $$groups_r{$pairs[4]}){		# if partitions of pair are both in groups (i.e. > min partition cutoff)
			if( $$groups_r{$pairs[1]} == $$groups_r{$pairs[4]} ){		# pair partitions are not in the same group
				print {$outfh{ $$groups_r{ $pairs[1]} }} 
					join("\n", join("\t", @pairs[0..1]), $pairs[2], join("\t", @pairs[3..4]), $pairs[5]), "\n";
				}
			else{
				print STDERR "Different groups: ", join(" <-> ", $$groups_r{$pairs[1]}, $$groups_r{$pairs[4]}), "\n";
				}
			}
		else{
			print STDERR " Not found in a group: ", join("\t", @pairs[0..1]), "\n"
				if ! exists $$groups_r{$pairs[1]};
			print STDERR " Not found in a group: ", join("\t", @pairs[3..4]), "\n"
				if ! exists $$groups_r{$pairs[4]};
			}
		#else{
		#	push(@{stats{
		#	}


		}
	
	# closing #
	close PAIR;
	map{ close $outfh{$_} or die $! } keys %$group_abunds_r;
	}

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
# placing partition clusters in groups based on max size #
	my ($clusters_r, $cluster_abunds_r, $max_size) = @_;

 	# status #
	print STDERR " Grouping partition clusters\n";

	# grouping #
	my %groups;
	my %group_abunds;
	my $group_cnt = 0;
	foreach my $cID (keys %$cluster_abunds_r){
		$group_cnt += $$cluster_abunds_r{$cID};						# summing cluster abundances, groupID defined by max_size
		
		foreach my $part (keys %{$$clusters_r{$cID}} ){				# all partitions in cluster
			my $gID = int($group_cnt / $max_size);
			$groups{$part} = $gID; 			# partition => group
			$group_abunds{$gID}++; 			# summing groups
			}
		}

	# status #
	print STDERR " Number of groups: ", scalar keys %group_abunds, "\n";
	die " ERROR: Number of groups exceeds 1000\n" if (scalar keys %group_abunds) > 1000;
	
		#print Dumper %groups; exit;
	return \%groups, \%group_abunds;			# partition => group
	}

sub load_mcl_out{
# loading mcl output into a hash #
	my ($mcl_out, $abunds_r) = @_;
	
	# status #
	print STDERR " Loading mcl output\n";
	
	open IN, $mcl_out or die $!;
	my %clusters;
	my %cluster_abunds;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		map{ $clusters{$.}{$_} = 1 } @line;			# all clusterID => partID
		map{ $cluster_abunds{$.} += $$abunds_r{$_} } @line;		# clusterID => sum(partitions abundances)
		}
	close IN;
	
		#print Dumper %cluster_abunds; exit;
		#print Dumper %clusters; exit;
	return \%clusters, \%cluster_abunds;				# partition => clusterID
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
# getting abundances of each partition #
	my ($basename, $min_part_size) = @_;
	
	# status #
	print STDERR " Getting partition abundances\n";
	
	# reading pair file and making pair hash #
	open PAIR, "$basename-pair.fna" or die $!;

	# making connectedness hash #
	#my $max_connections = 0;		# for normalizing connections
	my %abunds;

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
		
		# summing paritition abundances #
		$abunds{$line1[1]}++;					# number of reads per partition
		$abunds{$line3[1]}++;
		
		}
	
	seek(PAIR, 0, 0);

	# status #
	print STDERR " Making connectedness hash\n";

	my %connect;
	my $min_part_size_cnt = 0;
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
		
		# summing paritition abundances; not including partition if < min abundance #
		#$connect{$line1[1]}{$line3[1]}++ if $abunds{$line1[1]} >= $min_part_size && 
		#	$abunds{$line3[1]} >= $min_part_size && $line1[1] != $line3[1];
		if ($abunds{$line1[1]} >= $min_part_size && $abunds{$line3[1]} >= $min_part_size){
			$connect{$line1[1]}{$line3[1]}++ ;
			}
		else{ $min_part_size_cnt++; }		# counting partitions not included due to min_part_size
		}

	close PAIR;
	
	# stats #
	print STDERR " Number of read pairs not included due to '-min' cutoff: $min_part_size_cnt\n";
	
		#print Dumper %connect; exit;
		#print Dumper %abunds; exit;
	return \%connect, \%abunds;
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

