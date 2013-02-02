#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw/remove_tree/;
use Graph;
use Graph::Undirected;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $max_size = 1000000;
my $min_part_size = 10;
my $status_int = 1000000;
GetOptions(
	   "max-size=i" => \$max_size,						# max group size
	   "min-partition-size=i" => \$min_part_size,		# min partition size worth keeping
	   "status=i" => \$status_int,						# how often to provide a status update
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a *part file\n" if ! $ARGV[0];

### MAIN
my $basename = sort_pairs($ARGV[0]);

### making graph ###
my $g = make_paired_graph($basename);
merge_connected($g, $basename);

### writing out Moth files ###
#my ($connect_r, $max_connections, $names_r) = paired_connectivity($basename);
#write_names($basename, $names_r);
#$connect_r = normalize_connections($connect_r, $max_connections);
#write_connect_table($basename, $connect_r);
exit;

### writing out groups ###
my $parts_r = count_parts($basename);
my $groups_r = parts_in_groups($parts_r, $max_size, $min_part_size);
write_groups($groups_r, $basename);

### Subroutines
sub merge_connected{
# merging pairtitions connected by paired-end reads #
	my ($g, $basename) = @_;

	print STDERR " Finding connected components\n";
	my @cc = $g->connected_components();
	
	# status #
	print STDERR " Number of groups (connected partitions): ", scalar @cc, "\n";
	
	# component size distribution #
	print STDERR " Writting size distribution\n";
	my %size_dist;
	map{ $size_dist{scalar @$_}++ } @cc;
	print "\n", join("\t", qw/N_partitions N_groups/), "\n";
	foreach my $size (sort {$a <=> $b} keys %size_dist){
		print join("\t", $size, $size_dist{$size}), "\n";
		}
	print "\n"; 
	
	}

sub make_paired_graph{
# merging partitions based on paired-end reads #
	my ($basename) = @_;
	
	# status #
	print STDERR " Making graph\n";
	
	# reading pair file and making pair hash #
	open PAIR, "$basename-pair.fna" or die $!;

	# making graph object #
	my $g = Graph::Undirected->new;

	my $max_connections = 0;		# for normalizing connections
	while(<PAIR>){
		chomp;
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR " lines processed: ", ($. - 1) / 2, "\n";
			}
	
		# loading lines #
		my @line1 = split /\t/;
		my $line2 = <PAIR>;
		my @line3 = split /\t|\n/, <PAIR>;
		my $line4 = <PAIR>;
		
		# adding to merge hash #
		if($line1[1] ne $line3[1]){		# if partitions for pairs don't match, add to graph
			$g->add_edge($line1[1], $line3[1]);
			}
			
		}
	close PAIR;

		#print "The graph is $g\n";
	# graph stats #
	print STDERR " The graph has ", scalar $g->vertices, " vertices & ", scalar $g->edges, " edges\n";	
	
	return $g;
	}




sub write_groups{
# writing out partitions in groups #
	my ($groups_r, $basename) = @_;
	
	# getting unique groups #
	my %ugroups = map{$_, 1} values %$groups_r; 
	
	# status #
	print STDERR " Number of groups: ", scalar keys %ugroups, "\n";

	# output dir #
	my $dirname = $basename . "_groups";
	remove_tree("$dirname") if -d "$dirname";		# removing old dir if present
	mkdir "$dirname"; 
	print STDERR " Writing files to: $dirname\n";

	# status #
	print STDERR " Writting out paired reads\n";
	
	# opening filehandles for each group #
	my %fh;
	foreach my $group (keys %ugroups){
		open $fh{$group}, ">$dirname/group$group-pair.fna" or die $!;
		}
	
	# going through paired reads and writing to group #
	open PAIR, "$basename-pair.fna" or die $!;
	while(<PAIR>){
		my @line = split /\t/;
		my $line2 = <PAIR>;
		my $line3 = <PAIR>;
		my $line4 = <PAIR>;
		
		# 2nd pair placed w/ first #
		next if ! exists $$groups_r{$line[1]}; 		# if read is < cutoff
		print {$fh{ $$groups_r{$line[1]} }} join("", join("\t", @line), $line2, $line3, $line4);
		}

	# closing filehandles #
	close PAIR;
	map{ close $fh{$_} } keys %fh;
	
	
	# status #
	print STDERR " Writting out singleton reads\n";
	
	# going through singletons #
	foreach my $group (keys %ugroups){
		open $fh{$group}, ">$dirname/group$group-single.fna" or die $!;
		}
	
	# going through paired reads and writing to group #
	open SING, "$basename-single.fna" or die $!;
	while(<SING>){
		my @line = split /\t/;
		my $line2 = <SING>;
		
		next if ! exists $$groups_r{$line[1]}; 		# if read is < cutoff
		print {$fh{ $$groups_r{$line[1]} }} join("", join("\t", @line), $line2);
		}
	
	
	# closing filehandles #
	close SING;
	map{ close $fh{$_} } keys %fh;
	
	
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

sub count_parts{
### counting partitions in each file ###
	my ($basename) = @_;

	# status #
	print STDERR " Counting partitions\n";

	# counting pairs #	
	open PAIR, "$basename-pair.fna" or die $!;

	my %parts;
	while(<PAIR>){
		
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR " lines processed: ", ($. - 1) / 2, "\n";
			}
		
		my @line = split /\t/;
		$parts{$line[1]} += 2;
		for my $i (0..2){		# skipping 2nd in pair
			my $line = <PAIR>;
			}
		}
	close PAIR;
	
	# counting singles #
	open SING, "$basename-single.fna" or die $!;
	while(<SING>){
		
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR " lines processed: ", ($. - 1) / 2, "\n";
			}
		
		my @line = split /\t/;
		my $nline = <SING>;
		$parts{$line[1]}++;
		}
	close SING;


	return \%parts;	
	}





sub merge_merges{
# reciprically merging hashes until no more linkages found #
	my ($merge_hash_r) = @_;

	# status #
	print STDERR " Merging partitions in merge hash\n";
	
	
	my $N_group = scalar keys %$merge_hash_r;		# saving last group size (1-step memory)
	
	$merge_hash_r = merge_recip($merge_hash_r, $N_group);
	
	print Dumper $merge_hash_r; exit;
	
	## sub ##
	sub merge_recip{
	# reciprical merging #
		my ($merge_hash_r, $N_group) = @_;
	
		
		# merging if hash-lev1 value == hash-lev2 value (within any lev1 hash) #
		foreach my $query (keys %$merge_hash_r){		# checking other hashes for $merge
			my @to_merge;
			foreach my $merge (keys %$merge_hash_r){				# going through each merge-group and checking for query
				if( exists  $$merge_hash_r{$merge}{$query} ){		# could apply count number required HERE
					push(@to_merge, $merge); 						# list of merge-groups to merge
					next;
					}
				}
			
			if(@to_merge){
				print Dumper scalar @to_merge; 
				}
			
			# merging all merge-groups into the query group #
			foreach my $merging (@to_merge){
				foreach my $part (keys %{$$merge_hash_r{$merging}}){		# adding to the other group
					$$merge_hash_r{$query}{$part} = $$merge_hash_r{$merging}{$part};
					}
				delete $$merge_hash_r{$merging}; 		# removing the original group
				}
			}
			
		# checking number of groups after merging #
		my $N2_group = scalar keys %$merge_hash_r;
		
		# status #
		print STDERR " B4_merge: $N_group; After_merge: $N2_group\n";
		
		# running again if less groups than before #
		if($N2_group < $N_group){
			merge_recip($merge_hash_r, $N2_group);
			}
		else{
			return $merge_hash_r;
			}
		}
	
	}

sub write_connect_table{
# writing connectivity table (paired-ends spanning partitions) #
	my ($basename, $connect_r) = @_;
	
	open OUT, ">$basename-pairnet.dist" or die $!;
	
	foreach my $part1 (keys %$connect_r){
		foreach my $part2 (keys %{$$connect_r{$part1}}){
			print OUT join("\t", $part1, $part2, $$connect_r{$part1}{$part2}), "\n";
			}
		}
	close OUT;

	}

sub normalize_connections{
#  normalizing connections by max connection #
	my ($connect_r, $max_connections) = @_;
	
	foreach my $part1 (keys %$connect_r){
		foreach my $part2 (keys %{$$connect_r{$part1}}){
			$$connect_r{$part1}{$part2} = $$connect_r{$part1}{$part2} / $max_connections;
			}
		}

		#print Dumper $connect_r; exit;
	return $connect_r; 
	}

sub write_names{
# writing names file #
	my ($basename, $names_r) = @_;
	
	open OUT, ">$basename.names" or die $!;
	
	foreach (keys %$names_r){
		print OUT join("\t", $_, $_), "\n";
		}
	
	close OUT;
	}

sub paired_connectivity{
# merging partitions based on paired-end reads #
	my ($basename) = @_;
	
	# status #
	print STDERR " Making merge hash\n";
	
	# reading pair file and making pair hash #
	open PAIR, "$basename-pair.fna" or die $!;

	
	my %connect;
	my %names;
	my $max_connections = 0;		# for normalizing connections
	while(<PAIR>){
		chomp;
		# status #
		if(($. -1) % $status_int == 0){
			print STDERR " lines processed: ", ($. - 1) / 2, "\n";
			}
	
		# loading lines #
		my @line1 = split /\t/;
		my $line2 = <PAIR>;
		my @line3 = split /\t|\n/, <PAIR>;
		my $line4 = <PAIR>;
		
		# adding to merge hash #
		if($line1[1] ne $line3[1]){		# comparing partitions for pairs and adding to partition merge
			$connect{$line1[1]}{$line3[1]}++;
			$max_connections = $connect{$line1[1]}{$line3[1]} if $connect{$line1[1]}{$line3[1]} > $max_connections;
			}
			
		# adding to %names  for making a names file #
		$names{$line1[1]} = 1;
		$names{$line3[1]} = 1;
		}
	close PAIR;

		#print Dumper $max_connections; exit;	
		#print Dumper %connect; exit;
	return \%connect, $max_connections, \%names;
	}

sub sort_pairs{
# counting number of reads in each partition #
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
			print STDERR " lines processed: ", ($. - 1) / 2, "\n";
			}
		
		# load lines #
		my @line = split /\//;		# header, pair, partition
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

