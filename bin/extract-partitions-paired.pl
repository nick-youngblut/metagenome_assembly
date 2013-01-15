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
my $parts_r = count_parts($basename);
my $groups_r = parts_in_groups($parts_r, $max_size, $min_part_size);
write_groups($groups_r, $basename);

### Subroutines
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

