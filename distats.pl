#!/usr/bin/perl

use strict;
use warnings;
use autodie;
use File::Temp;
use Data::Dumper;

# Implements getopt function GetOpt(), parses command line from @ARGV, recognizing
# and removing specified options and their possible values
use Getopt::Long;

# Declare option variables, set value to zero (=false) or one (=true)
my $real_values = 0; # default: distances in percent with 2 decimal places
my $treat_subspecies_as_separate = 1; # default: treat subspecies as distinct species
my $distance = 'p'; # defaults to p distance
my $num_threads = 1; # defaults to single-threaded
my $print_dist_matrix = 0; # default: distance matrix is not printed into in an extra file

# Set flag
GetOptions(
	'real_values'   => \$real_values,
	'subspecies!'   => \$treat_subspecies_as_separate,
	'distance=s'    => \$distance,
	'num_threads=i' => \$num_threads,
	'print_dist_matrix' 	=> \$print_dist_matrix,
) or die "Unknown option";

# Usage parameters
my $usage = "USAGE: perl $0 [OPTIONS] FASTA_FILE (.fasta) OUTPUTFILE (.txt or .csv)\n";
$usage .= "Options:\n";
$usage .= "  --real_values   original results with endless decimal places, no option: results in percent with two decimal places\n";
$usage .= "  --distance=X    use different distance matrix X. Possible values are 'p' and 'k2p'. Default: 'p'\n";
$usage .= "  --num_threads=N use N threads for calculating the distance matrix\n";
$usage .= "  --nosubspecies  do not treat subspecies as separate species; they are merged into one species\n";
$usage .= "  --print_dist_matrix	 print distance matrix into a file";
$usage .= "\n";

# If number of arguments in command line is not 2, show
# $usage message and exit script.
if (scalar @ARGV != 2) { print "\n", $usage and exit 1 };

# See whether Parallel::ForkManager is installed, if not print error
# message which prompts to install package

my $forkmanager = eval {
	require Parallel::ForkManager;
	Parallel::ForkManager->new( $num_threads );
};
unless ($forkmanager) {
	if ($num_threads != 1) {
		print "Error: Requested to multi-thread but Parallel::ForkManager not usable\n" and exit 1;	
	}
}

# Forces to exit script if distance is wrongly indicated. Currently, 
# distance can be 'p' or 'k2p'. Default is 'p' distance (see above).

if ($distance ne 'p' and $distance ne 'k2p') { print $usage and exit 1 };

# Declare files: input_file must be fasta, output_file a .txt or .csv file.
my $fasta_file = shift @ARGV;
my $output_file = shift @ARGV;


# fasta2hash: check if file is fasta, read sequence headers and sequences
# to keys and values of a hash; loop through keys of hash, check whether
# comparison between two keys already happened or not, if not: compare seqs
# and calculate distances, write into output file.


print "Reading Fasta file '$fasta_file'... ";
Seqload::Fasta::check_if_fasta($fasta_file) or die "Not a valid fasta file: $fasta_file\n";
my $header2seq = Seqload::Fasta::slurp_fasta($fasta_file);
print "Completed.\n";

my $fh_dist_matrix;

if ($print_dist_matrix){	
	# open output file for writing distance "matrix"
	open $fh_dist_matrix, '>', 'distance_matrix_' . $output_file;
	
	print {$fh_dist_matrix} join ("\t",
		"Species_1",
		"Header_1",
		"Species_2",
		"Header_2",
		"Distance",
	), "\n";
}
print "Printed distance matrix.\n";
my $data = {};

# setup the callback routine for parallel processes
# this adds its child data structure to the global hashref $data
if ($forkmanager) {
	$forkmanager->run_on_finish( sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
		if (defined $data_structure_reference) {
			$data->{$ident} = $data_structure_reference;
		}
	});
}

# do the pairwise calculations
print "Calculating pairwise distances... ";

# prepare progress reporter
my $number_of_sequences = scalar keys %$header2seq;
my $last_time = time;
my $current_time = $last_time;
my $counter = 0;

foreach my $header_1 (sort { $a cmp $b } keys %$header2seq) {

	# report on progress
	$counter++;
	$current_time = time;
	printf "%.2f%%\rCalculating pairwise distances... ", $counter / $number_of_sequences * 100;
	$last_time = $current_time;

	# fork this
	if ($forkmanager) { $forkmanager->start($header_1) and next }

	# data structure for collecting child data
	my $child_data = { };

	foreach my $header_2 (sort { $a cmp $b } keys %$header2seq) {

		next if ($header_1 eq $header_2);
		# even if "next if $data->{header2}->{header1}" already present
		

		my $dist = dist($header2seq->{$header_1}, $header2seq->{$header_2});

		$child_data->{$header_2} = $dist;
		
		
		
		# do this only if $dist_matrix = 1
		if ($print_dist_matrix){
			
			my ($genus_1, $species_1) = parse_header($header_1);
			my ($genus_2, $species_2) = parse_header($header_2);
			
			unless (defined $data->{$header_2}->{$header_1}) {
				print {$fh_dist_matrix} join ("\t",
					$genus_1 . '_' . $species_1,
					$header_1,
					$genus_2 . '_' . $species_2,
					$header_2,
					$dist,
				), "\n";
			}
		}	
	}

	# catch fork
	if ($forkmanager) { $data->{$header_1} = $forkmanager->finish(0, $child_data) }
	else { $data->{$header_1} = $child_data }

}

# wait for all children
if ($forkmanager) { $forkmanager->wait_all_children }

if ($print_dist_matrix){
	close $fh_dist_matrix;
}

# data collection complete

print "Completed.\n";
print "Categorizing... ";

# open output file for writing statistics table
open my $fh_out, '>', $output_file;

# print header in output file
print {$fh_out} join ("\t",
	"Taxon",
	"#_of_intrasp_distances",
	"Intrasp_dist_Mean",
	"Intrasp_dist_Min",
	"Intrasp_dist_Median",
	"Intrasp_dist_Max",
	"Closest_Sp",
	"Closest_Sp_Min_Dist",
	"Closest_Sp_Max_Dist",
	"Closest_Sp_Dist_Median",
	"#_distances_to_Closest_Sp",
	"most_distant_congener",
	"most_distant_congener_Max_Dist",
), "\n";

my $distances;

# loop through large data structure, allocate each distance value to 
# corresponding category and save it in hash

foreach my $header_1 (sort {$a cmp $b} keys %$data){
	
#~ while (my ($header_1, $data_for_header_1) = each %$data) {

	# get species out of sample ID
	my ($genus_1, $species_1) = parse_header($header_1);
	
	# we go through each distance to this species
	foreach my $header_2 (sort {$a cmp $b} keys %{$data->{$header_1}}){
			
	#~ while (my ($header_2, $dist) = each %$data_for_header_1) {
		
		# get species out of sample ID
		my ($genus_2, $species_2) = parse_header($header_2);
		
		# determine category
		my $cat = categorize($genus_1, $species_1, $genus_2, $species_2);
		
		# save
		$distances->{$genus_1 . '_' . $species_1}->{$cat}->{$genus_2 . '_' . $species_2}->{$header_1 . '_' . $header_2} = $data->{$header_1}->{$header_2};
		
	}
	
}

print "Completed.\n";

my $min_intrasp_dist;
my $median_intrasp_dist;
my $mean_intrasp_dist;
my $max_intrasp_dist;
my $intrasp_nr_of_dist;

my $name_closest_sp;
my $min_dist_closest_sp;
my $max_dist_closest_sp;
my $median_dist_closest_sp;
my $mean_dist_closest_sp;


print "Calculating distance statistics... ";
foreach my $species (sort {$a cmp $b} keys %$distances) {

	$name_closest_sp = '';
	$min_intrasp_dist = '';
	$median_intrasp_dist = '';
	$mean_intrasp_dist = '';
	$max_intrasp_dist = '';

	my $dist_cat_1 = $distances->{$species}->{1};
	my $dist_cat_2 = $distances->{$species}->{2};
	my $dist_cat_3 = $distances->{$species}->{3};	
	$intrasp_nr_of_dist = 0;
	
	# determine intraspecific distances
	if (defined $dist_cat_1){

		my $list_of_intrasp_dist = get_distances($dist_cat_1);

		$intrasp_nr_of_dist  = scalar @$list_of_intrasp_dist;
		$min_intrasp_dist    = min($list_of_intrasp_dist);
		$median_intrasp_dist = median($list_of_intrasp_dist);
		$mean_intrasp_dist   = mean($list_of_intrasp_dist);
		$max_intrasp_dist    = max($list_of_intrasp_dist);

	}

	# determine congeners and largest intra-generic distance
	my $list_of_congener_and_dist = [];
	my ($most_distant_congener, $most_distant_congener_dist);

	if (defined $dist_cat_2){

		($most_distant_congener, $most_distant_congener_dist) = get_most_distant_congener($dist_cat_2);

	}


	# closest species in cat 2 and 3
	# combine the hashes of cat 2 and 3
	my $other_distances = { };
	if (defined $dist_cat_2) {
		$other_distances = { %$other_distances, map { $_ => $dist_cat_2->{$_} } keys %$dist_cat_2 };
	}
	if (defined $dist_cat_3) {
		$other_distances = { %$other_distances, map { $_ => $dist_cat_3->{$_} } keys %$dist_cat_3 };
	}

	# extract distances for closest species
	my $distances_for_closest_species = get_closest_species($other_distances);

	foreach my $closest_species (sort {$a cmp $b} keys %$distances_for_closest_species){
		
		# get distances for closest species

		my $distances_of_closest_sp = get_distances_for_closest_sp($distances_for_closest_species->{$closest_species});

		$name_closest_sp = $closest_species;
		
		# determine min, max, median
		$min_dist_closest_sp = min($distances_of_closest_sp);
		
		$max_dist_closest_sp = max($distances_of_closest_sp);
		
		$median_dist_closest_sp = median($distances_of_closest_sp);
		
		# print output!
		print_table_row(
			$species,
			$intrasp_nr_of_dist,
			$min_intrasp_dist,
			$max_intrasp_dist,
			$mean_intrasp_dist,
			$median_intrasp_dist,
			$name_closest_sp,
			$min_dist_closest_sp,
			$max_dist_closest_sp,
			$median_dist_closest_sp,
			$distances_of_closest_sp,
			$most_distant_congener,
			$most_distant_congener_dist
		); 
	}
	
}

close $fh_out;

print "Completed.\n";

exit;

# We all live in a yellow subroutine, yellow subroutine, yellow subroutine... 

# parse a fasta header, return (genus, species, ID)
sub parse_header {
	# parse out species name
	# headers are formatted like "species|foo bar baz"
		
	my $header = shift;
	my $fields;
	
	my ($taxon, $ID) = split /\|/, $header;
	
	if ($treat_subspecies_as_separate) { $fields = 2 } else { $fields = 3 }
	
	my ($genus, $species) = split "\_", $taxon, $fields;
	
		
	return($genus, $species, $ID);
	
}
# categorize two species,
# input: taxonA genus, taxonA species, taxonB genus, taxonB species
# output: 1, 2, 3, or 0 (could not be categorized; never happens)
sub categorize {
	my $taxonAgenus = shift;
	my $taxonAspecies = shift;
	my $taxonBgenus = shift;
	my $taxonBspecies = shift;
	
		
	if (($taxonAgenus eq $taxonBgenus) and ($taxonAspecies eq $taxonBspecies)){
		return 1;
	}
	
	elsif (($taxonAgenus eq $taxonBgenus) and ($taxonAspecies ne $taxonBspecies)){
		return 2;
	}
	
	elsif ($taxonAgenus ne $taxonBgenus){
		return 3;
	}
	else{
		return 0;
	}
}


# return minimal value of a list
# argument: array reference
# usage:
# my $minimum = min( $array_ref );

sub min {
	
	my $distances = shift @_;
	my $min = $$distances[0];
 for (@$distances){
	 $min = $_ if $_ < $min
	 }
	
	 return $min;
 }
# return maximal value of a list
# argument: array reference
# usage:
# my $maximum = max( $array_ref );
sub max {
	
	 my $distances = shift @_;
	 my $max = $$distances[0];
	
	 for (@$distances){
		 $max = $_ if $_ > $max;
	 }
	
	 return $max;
 }
 
# return mean of a list of values 
# argument: array reference
# usage:
# my $mean = mean( $array_ref );
 sub mean {
	
	 my $distances = shift @_;
	 my $sum;
	 my $mean;
	
	 for( @$distances ){
		 $sum += $_;
	 }
	
	 $mean = ($sum/scalar @$distances);

	 return $mean;
 }

# return median of a list of values
# argument: array reference
# usage:
# my $median = median( $array_ref )
sub median {

	my $distances = shift @_;
	
	$distances = [ sort {$a <=> $b} @$distances ];
	
	my $nr_of_dist = scalar @$distances;
	
	if ($nr_of_dist == 0) {
		return undef;
	}

	my $odd = $nr_of_dist%2;
	
	if( $odd ){
			
		return $distances->[($nr_of_dist-1)/2];
	}

	else{

		return (( $distances->[$nr_of_dist/2]) + ($distances->[($nr_of_dist/2)-1]) )/2;
	}

}

# return distances, either Hamming (p) or Kimura 2-parameter (k2p)
sub dist {
	# corrected Hamming (p) distance 
	if ($distance eq 'p') {
		return pdist( $_[0], $_[1] );
	}
	# Kimura 2-parameter distance
	elsif ($distance eq 'k2p') {
		return k2pdist( $_[0], $_[1] );
	}
}

# Calculate and return k2p distance
sub k2pdist {
	my $seqA = shift;
	my $seqB = shift;
	
	my $shorter_seq;
	
	if (length $seqA > length $seqB){
		$shorter_seq = $seqB;
	}
	else{
		$shorter_seq = $seqA;
	}
	
	# die "Unequal sequence lengths in k2p comparison\n" if length($seqA) != length($seqB);
	my $s_equal    = 0;
	my $s_transit  = 0;
	my $s_transv   = 0;
	my $s_ambig    = 0;
	for (my $i = 0; $i < length $seqA; $i++) {
		my $nucA = substr $seqA, $i, 1;
		my $nucB = substr $seqB, $i, 1;
		# equal sites
		if (lc $nucA eq lc $nucB) {
			++$s_equal;
		}
		# transition (A <-> G, C <-> T)
		elsif (transition($nucA, $nucB)) { 
			++$s_transit;
		}
		# transversion (T <-> G, A <-> C, A <-> T, C <-> G)
		elsif (transversion($nucA, $nucB)) { 
			++$s_transv;
		}
		# ambiguous site (N or -)
		elsif ($nucA =~ /[nN-]/ or $nucB =~ /[nN-]/) {
			$s_ambig++;
		}
	}
	# transition frequency
	my $p = $s_transit / length $shorter_seq;
	# transversion frequency
	my $q = $s_transv  / length $shorter_seq;
	# kimura two-parameter distance
	# return -0.5 * log(1 - 2 * $p - $q) - 0.25 * log(1 - 2 * $q);
	return -0.5 * log((1 - 2 * $p - $q) * sqrt(1 - 2 * $q));
}

# Calculate and return number of transversions in sequence
sub transversion {
	my $n1 = shift;
	my $n2 = shift;
	if (lc $n1 eq 'g' and lc $n2 eq 't') { return 1 }
	elsif (lc $n1 eq 't' and lc $n2 eq 'g') { return 1}
	elsif (lc $n1 eq 'a' and lc $n2 eq 'c') { return 1}
	elsif (lc $n1 eq 'c' and lc $n2 eq 'a') { return 1}
	elsif (lc $n1 eq 'a' and lc $n2 eq 't') { return 1}
	elsif (lc $n1 eq 't' and lc $n2 eq 'a') { return 1}
	elsif (lc $n1 eq 'c' and lc $n2 eq 'g') { return 1}
	elsif (lc $n1 eq 'g' and lc $n2 eq 'c') { return 1}
	else { return 0 }
}

# Calculate and return number of transitions in sequence

sub transition {
	my $n1 = shift;
	my $n2 = shift;
	if (lc $n1 eq 'a' and lc $n2 eq 'g') { return 1 }
	elsif (lc $n1 eq 'g' and lc $n2 eq 'a') { return 1}
	elsif (lc $n1 eq 'c' and lc $n2 eq 't') { return 1 }
	elsif (lc $n1 eq 't' and lc $n2 eq 'c') { return 1}
	else { return 0 }
}

# Calculate and return p distance
sub pdist {
	my $seqA = shift;
	my $seqB = shift;
	my $s_equal = 0;
	my $s_diff  = 0;
	my $s_ambig = 0;
	
	my $shorter_seq;
	
	if (length $seqA > length $seqB){
		$shorter_seq = $seqB;
	}
	else{
		$shorter_seq = $seqA;
	}
		
	for (my $i = 0; $i < length $shorter_seq; $i++) {
		my $nucA = substr $seqA, $i, 1;
		my $nucB = substr $seqB, $i, 1;
		if (lc $nucA eq lc $nucB) {
			$s_equal++;
		}
		elsif ($nucA =~ /[NRYWSMKHBVD-]/i or $nucB =~ /[NRYWSMKHBVD-]/i) {
			$s_ambig++;
		}
		else {
			$s_diff++;
		}
	}
	return $s_diff/(length($shorter_seq)-$s_ambig);
}

# printing subroutine

sub print_table_row {
	# input values
	my ($taxon_1, $intrasp_nr_of_dist, $min_intrasp_dist, $max_intrasp_dist, $mean_intrasp_dist, $median_intrasp_dist, $closest_sp, $min_dist, $max_dist, $median_dist, $closest_sp_dist, $most_distant_congener, $dist_to_most_distant_congener) = @_;
	my $factor = 100;
	my $places = 2;
	if ($real_values){
		$factor = 1;
		$places = 18;
	}
	printf {$fh_out} join( "\t", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s"  ) . "\n", 
		$taxon_1, 
		$intrasp_nr_of_dist     ? $intrasp_nr_of_dist                                    : 0, 
		$intrasp_nr_of_dist     ? sprintf("%.${places}f", $mean_intrasp_dist*$factor)    : '',
		$intrasp_nr_of_dist     ? sprintf("%.${places}f", $min_intrasp_dist*$factor)     : '',
		$intrasp_nr_of_dist     ? sprintf("%.${places}f", $median_intrasp_dist*$factor)  : '',
		$intrasp_nr_of_dist     ? sprintf("%.${places}f", $max_intrasp_dist*$factor)     : '',
		# closest sp
		defined $closest_sp ? $closest_sp : '',
		# closest sp dist = min dist
		( defined $closest_sp && defined $min_dist ) ? sprintf("%.${places}f", $min_dist*$factor) : '',
		# closest sp max dist
		( defined $closest_sp && defined $max_dist ) ? sprintf("%.${places}f", $max_dist*$factor) : '',
		# closest sp median dist
		( defined $closest_sp && defined $median_dist ) ? sprintf("%.${places}f", $median_dist*$factor) : '',
		# nr of dist to closest sp
		defined $closest_sp ? scalar( @$closest_sp_dist ) : '',
		# most dist congener
		$most_distant_congener ? $most_distant_congener : '',
		# max dist of most dist congener
		$most_distant_congener ? sprintf("%.${places}f", $dist_to_most_distant_congener*$factor) : '',
	;
}

# get most distant congener (= of same genus but not same species)

sub get_most_distant_congener {
	my $hashref_of_congeners_and_dist = shift;
	my $most_distant_congener = '';
	my $dist_to_most_distant_congener = 0; 
	
	foreach my $congener (keys %$hashref_of_congeners_and_dist) {

		while (my ($header, $dist) = each %{$hashref_of_congeners_and_dist->{$congener}}) {

			if ($dist > $dist_to_most_distant_congener) {
				$dist_to_most_distant_congener = $dist;
				$most_distant_congener = $congener;
			}

		}

	}
	# print "most distant congener: $most_distant_congener, dist $dist_to_most_distant_congener\n";
	return ($most_distant_congener, $dist_to_most_distant_congener);
}

# input: hashref
# return an arrayref of all distances (i.e., of all species in this hashref)
sub get_distances {

	my $hashref = shift;
	
	my $list_of_distances = [ ];
	
	foreach my $species (keys %{$hashref}) {
		while (my ($header, $dist) = each %{$hashref->{$species}}){
		
			push @$list_of_distances, $dist;
		
		}
	}
	
	return $list_of_distances;

}

# get closest species

sub get_closest_species {

	my $distances_for = shift;
	my $closest_species = { };
	my $min_dist = 99;

	# determine minimum distance
	foreach my $species (keys %$distances_for){
		
		while (my ($header, $dist) = each %{$distances_for->{$species}}) {
			if ($dist < $min_dist) {
				$min_dist = $dist;
			}
		}
	}
	
	# pick out all that have this min distance
	foreach my $species (keys %$distances_for){
		while (my ($header, $dist) = each %{$distances_for->{$species}}) {
			if ($dist == $min_dist) {
				$closest_species->{$species} = $distances_for->{$species};
				
			}
		}
	}
	
	return $closest_species;

}

# Calculate and return distances for closest species

sub get_distances_for_closest_sp {
	my $hashref = shift;
	

	my $list_of_distances = [ ];
	
	while (my ($header, $dist) = each %$hashref){
	
		push @$list_of_distances, $dist;
		
	}
	return $list_of_distances;
}



# Copyright 2015 Malte Petersen <mptrsen@uni-bonn.de>
#
# This module is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This module is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this module. If not, see http://www.gnu.org/licenses/.
#--------------------------------------------------
package Seqload::Fasta;
use strict;
use warnings;
use Carp;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( fasta2csv check_if_fasta );

# Constructor. Returns a sequence database object.
sub open {
  my ($class, $filename) = @_;
  open (my $fh, '<', $filename)
    or confess "Fatal: Could not open $filename\: $!\n";
  my $self = {
    'filename' => $filename,
    'fh'       => $fh
  };
  bless($self, $class);
  return $self;
}

# Returns the next sequence as an array (hdr, seq). 
# Useful for looping through a seq database.
sub next_seq {
  my $self = shift;
  my $fh = $self->{'fh'};
	# this is the trick that makes this work
  local $/ = "\n>"; # change the line separator
  return unless defined(my $item = readline($fh));  # read the line(s)
  chomp $item;
  
  if ($. == 1 and $item !~ /^>/) {  # first line is not a header
    croak "Fatal: " . $self->{'filename'} . " is not a FASTA file: Missing descriptor line\n";
  }

	# remove the '>'
  $item =~ s/^>//;

	# split to a maximum of two items (header, sequence)
  my ($hdr, $seq) = split(/\n/, $item, 2);
	$hdr =~ s/\s+$//;	# remove all trailing whitespace
  $seq =~ s/>//g if defined $seq;
  $seq =~ s/\s+//g if defined $seq; # remove all whitespace, including newlines

  return($hdr, $seq);
}

# Closes the file and undefs the database object.
sub close {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $filename = $self->{'filename'};
  close($fh) or carp("Warning: Could not close $filename\: $!\n");
  undef($self);
}

# Destructor. This is called when you undef() an object
sub DESTROY {
  my $self = shift;
  $self->close;
}


# validates a fasta file by looking at the FIRST (header, sequence) pair
# arguments: scalar string path to file
# returns: true on validation, false otherwise
sub check_if_fasta {
	my $infile = shift;
	my $infh = Seqload::Fasta->open($infile);
	my ($h, $s) = $infh->next_seq() or return 0;
	return 1;
}

# loads a Fasta file into a hashref
# arguments: scalar string path to file
# returns: hashref (header => sequence)
sub slurp_fasta {
	my $infile = shift;
	my $sequences = {};
	my $infh = Seqload::Fasta->open($infile);
	while (my ($h, $s) = $infh->next_seq()) {
		$sequences->{$h} = $s;
	}
	undef $infh;
	return $sequences;
}
