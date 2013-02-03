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

my ($verbose, $dump_connect, $write_abunds);
my $max_size = 1000000;
my $min_part_size = 5;
my $status_int = 1000000;
my $inflation = 5;				# inflation param default
my $threads = 1;
GetOptions(
	   "max-size=i" => \$max_size,						# max group size
	   "min-partition-size=i" => \$min_part_size,		# min partition size worth keeping
	   "status=i" => \$status_int,						# mcl inflation param
	   "threads=i" => \$threads,						# how often to provide a status update
	   "Inflation=f" => \$inflation,					# number of threads for mcl
	   "dump" => \$dump_connect,						# dumping connected graph (for running mcl outside of script)
	   "abundance" => \$write_abunds,						# writing cluster abundance distribution
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a *part file\n" if ! $ARGV[0];

### MAIN	
# finding paired-end reads and making connection hash #
my $basename = sort_pairs($ARGV[0]);													# parsing out pairs
my ($connect_r, $abunds_r) = partition_connection_hash($basename, $min_part_size);		# making a parition graph for loading into mcl
dump_connect_hash($connect_r, $basename) if $dump_connect;

# mcl clustering #
my $mcl_out = call_mcl($connect_r, $basename, $inflation, $threads);					# clustering w/ mcl
my ($clusters_r, $cluster_abunds_r) = load_mcl_out($mcl_out, $abunds_r);
write_cluster_abunds($cluster_abunds_r, $basename) if $write_abunds;

# grouping clustering for writing out read files #
my ($groups_r, $group_abunds_r) = parts_in_groups($clusters_r, $cluster_abunds_r, $max_size);				# placing paritions into groups
write_groups($groups_r, $group_abunds_r, $basename);


### Subroutines
sub write_groups{
# writing out groups (group => clusters => partitions) #
	my ($groups_r, $group_abunds_r, $basename) = @_;
	
	# status #
	print STDERR "...Writing group files\n";
	
	# I/O #
	## read paired-end reads ##
	open PAIR, "$basename-pair.fna" or die $!;

	## making directory for group files ##
	remove_tree("$basename\_groups") if -d "$basename\_groups";
	mkdir "$basename\_groups" or die $!;
	my %outfh;
	map{ open $outfh{$_}, ">$basename\_groups/$basename\_g$_.fna" or die $! } keys %$group_abunds_r;
	
	open LOG, ">$basename.log" or die $!;
	
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
		
		if(exists $$groups_r{$pairs[1]} && exists $$groups_r{$pairs[4]}){		# if partitions of pair are both in groups (i.e. > min partition cutoff)
			if( $$groups_r{$pairs[1]} == $$groups_r{$pairs[4]} ){		# pair partitions are not in the same group
				print {$outfh{ $$groups_r{ $pairs[1]} }} 
					join("\n", join("\t", @pairs[0..1]), $pairs[2], join("\t", @pairs[3..4]), $pairs[5]), "\n";
				}
			else{
				print LOG "Different groups:\t", join(" <-> ", $$groups_r{$pairs[1]}, $$groups_r{$pairs[4]}), "\n";
				}
			}
		else{
			print LOG "Not in a group; partition < '-min':\t", join("\t", @pairs[0..1]), "\n"
				if ! exists $$groups_r{$pairs[1]};
			print LOG "Not in a group; partition < '-min':\t", join("\t", @pairs[3..4]), "\n"
				if ! exists $$groups_r{$pairs[4]};
			}
		}
	
	# closing #
	close PAIR;
	map{ close $outfh{$_} or die $! } keys %$group_abunds_r;
	close LOG;
	
	print STDERR "\tGroup files written to: '$basename\_groups/'\n";
	print STDERR "\tLog file written to: '$basename.log'\n";
	}

sub dump_connect_hash{
# writting out connected hash to file for running mcl outside of script #
	my ($connect_r, $basename) = @_;

	# making an output file name #
	my $dump_out = $basename . "_abc.txt";
	
	# status #
	print STDERR "...Dumping partition connection hash: '$dump_out'\n";
	
	open OUT, ">$dump_out" or die $!;
	
	# running mcl #
	foreach my $c1 (keys %$connect_r){
		foreach my $c2 (keys %{$$connect_r{$c1}}){
			print OUT join(" ", $c1, $c2, $$connect_r{$c1}{$c2}), "\n";
			}
		}	
	close OUT;
	
	print STDERR "\tConnection hash dumped. Exiting\n";
	exit;
	}

sub parts_in_groups{
# placing partition clusters in groups based on max size #
	my ($clusters_r, $cluster_abunds_r, $max_size) = @_;

 	# status #
	print STDERR "...Grouping partition clusters\n";

	# grouping #
	my %groups;
	my %group_abunds;
	my $group_cnt = 0;
	foreach my $cID (sort {$$cluster_abunds_r{$b} <=> $$cluster_abunds_r{$a}} 
		keys %$cluster_abunds_r){
		
		$group_cnt += $$cluster_abunds_r{$cID};						# summing cluster abundances, groupID defined by max_size
		
		foreach my $part (keys %{$$clusters_r{$cID}} ){				# all partitions in cluster
			my $gID = int($group_cnt / $max_size);
			$groups{$part} = $gID; 									# partition => group
			$group_abunds{$gID}++; 									# summing groups
			}
		}

	# status #
	print STDERR "\tNumber of groups: ", scalar keys %group_abunds, "\n";
	die " ERROR: Number of groups exceeds 1000\n" if (scalar keys %group_abunds) > 1000;
	
		#print Dumper %groups; exit;
	return \%groups, \%group_abunds;			# partition => group
	}

sub write_cluster_abunds{
	my ($clusters_abunds_r, $basename) = @_;
	
	# status #
	print STDERR "...Writing cluster abundances distribution\n";
	
	open OUT, ">$basename\_clust-abund.txt";
	print OUT join("\t", qw/Cluster_id N_partitions/), "\n";
	foreach my $cID (sort {$$cluster_abunds_r{$b} <=> $$cluster_abunds_r{$a}} 
		keys %$cluster_abunds_r){
		print OUT join("\t", $cID, $$cluster_abunds_r{$cID}), "\n";
		}
	close OUT;
	
	print STDERR "\tFile written to: $basename\_clust-abund.txt\n";
	}

sub load_mcl_out{
# loading mcl output into a hash #
	my ($mcl_out, $abunds_r) = @_;
	
	# status #
	print STDERR "\n...Loading mcl output\n";
	
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
	
	# status #
	print STDERR "\tNumber of clusters: ", scalar keys %cluster_abunds, "\n";
		#print Dumper %cluster_abunds; exit;
		#print Dumper %clusters; exit;
	return \%clusters, \%cluster_abunds;				# partition => clusterID
	}

sub call_mcl{
# calling mcl for clustering #
	my ($connect_r, $basename, $inflation, $threads) = @_;
	
	# status #
	print STDERR "...Calling mcl\n\n";
	
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
	print STDERR "...Getting partition abundances\n";
	
	# reading pair file and making pair hash #
	open PAIR, "$basename-pair.fna" or die $!;

	# making connectedness hash #
	#my $max_connections = 0;		# for normalizing connections
	my %abunds;
	while(<PAIR>){
		chomp;
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR "\tRead pairs processed: ", ($. - 1) / 2, "\n";
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
	print STDERR "...Making connectedness hash\n";

	my %connect;
	my $line_cnt = 0;
	while(<PAIR>){
		chomp;
		# status #
		if(($line_cnt) % $status_int == 0){
			print STDERR "\tRead pairs processed: ", ($line_cnt) / 2, "\n";
			}
	
		# loading lines #
		my @line1 = split /\t/;
		my $line2 = <PAIR>;
		my @line3 = split /\t|\n/, <PAIR>;
		my $line4 = <PAIR>;
		
		# summing paritition abundances; not including partition if < min abundance #
		if ($abunds{$line1[1]} >= $min_part_size && $abunds{$line3[1]} >= $min_part_size){
			$connect{$line1[1]}{$line3[1]}++ ;
			}
		
		# reads processed
		$line_cnt += 4;
		}

	close PAIR;
	
	# getting number of partitions below cutoff #
	my $below_cnt = 0;
	foreach my $part (keys %abunds){
		$below_cnt++ if $abunds{$part} < $min_part_size;
		}
	
	# stats #
	print STDERR "\n### Partition filtering stats ###\n";
	print STDERR "Number of partitions: ", scalar keys %abunds, "\n";
	print STDERR "Number of partitions < '-min' cutoff: $below_cnt\n";
	print STDERR "% partitions remaining: ", 
		sprintf("%.0f", ((scalar keys %abunds) - $below_cnt) / (scalar keys %abunds) * 100), "%\n\n";
	
		#print Dumper %connect; exit;			# partID => number of paired-end connections
		#print Dumper %abunds; exit;
	return \%connect, \%abunds;
	}

sub sort_pairs{
# parsing out existing paired-end reads #
	my ($infile) = @_;
	
	# status #
	print STDERR "...Splitting reads into pairs and singletons\n";
	
	# I/O #
	(my $basename = $infile) =~ s/\.[^\.]+$|$//;
	open IN, $infile or die $!;
	open PAIR, ">$basename-pair.fna" or die $!;
	open SING, ">$basename-single.fna" or die $!;

	my %pairs;
	my $pair_cnt = 0;			# counting number of pairs
	my $read_cnt = 0;			# total number of reads
	while(<IN>){	
		$read_cnt++;
		
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR "\tReads processed: ", ($. - 1) / 2, "\n";
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
			$pair_cnt++;
			}
		else{
			$pairs{$line[0]} = join("", join("/", @line), $nline);
			}

		}
	close IN;
	
	# stats #
	print STDERR "\n### Paired-reads stats ###\n";
	print STDERR "Number of reads: $read_cnt\n";
	print STDERR "Number of paired reads: ", $pair_cnt * 2, "\n";
	print STDERR "% paired reads: ", sprintf("%.1f", ($pair_cnt * 2)/ $read_cnt * 100), "%\n\n";
	
	# writing out all singletons #
	print STDERR "...Writing out all singletons\n";
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

Minimum partition size (number of reads) to retain for mcl clustering. [5]

=item --dump

Write out the table used for mcl clustering?. [FALSE]

=item --threads

Number of threads for mcl clustering. [1]

=item --Inflation

Inflation value for mcl clustering. [5]

=item --abundance

Write out the distribution of partitions in each cluster? [FALSE] 

=item -s 

Status output interval (number of reads). [1000000]

=item -h

This help message

=back

=head2 For more information:

perldoc extract-partitions-paired.pl

=head1 DESCRIPTION

Extract the partitions and retain paired-end reads for assembly with an assembler that
requires paired-end reads.

Partitions are clustered by paired-end reads spanning partitions. mcl is used for 
the clusterin process with a high default inflation parameter.

The clusters of partitions are then grouped for writing to files (--min-partition-size).

=head2 Reasoning for clustering partitions

Paired-end reads should be from the same molecule and thus should
come from the same organism. Sequencing or khmer partitioning artifacts could cause
paired-end reads to span partitions seperating different organisms. If this occurs too often,
all of the partitions will agglomerate. You may have to adjust the inflation parameter if this
happens.

=head1 EXAMPLES

=head2 Basic usage

extract-partitions-paired.pl iowa-corn-50m.fa.gz.part

=head2 Dumping table for mcl clustering

extract-partitions-paired.pl -d iowa-corn-50m.fa.gz.part

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/metagenome_assembly/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

