# DiStats - calculate distance statistics for a DNA sequence alignment

## SYNOPSIS

	distats.pl [OPTIONS] ... INPUT_FILE OUTPUT_FILE

## DESCRIPTION

DiStats calculates distance statistics for all sequences in an input FASTA file. By default, DiStats calculates the p distance and prints a table with the results as percent values with two decimal places. 

## INPUT AND OUTPUT

The input file must be a FASTA file containing a nucleotide sequence alignment with sequences of the same length. The header must have the following format: `>Genus_species|Sequence0001` - i.e., first the species name with genus and species separated by an underscore, and after the pipe, a unique sequence identifier (a plain unique number or whatever) is necessary. 

The output file is a tab-delimited table and should be named accordingly, e.g., output.txt. 


## COMMAND LINE OPTIONS

	--real_values

Output distances with 15 decimal places. Default: rounded to two decimal places. 

	--distance[=X]

Use different distance matrix X. Possible values are 'p' (corrected Hamming distance) and 'k2p' (Kimura 2-parameter distance). Default: p distance.

	--num_threads[=N]

Use N threads for calculating the distance matrix. Default: single-threaded.

	--nosubspecies

Do not treat subspecies as distinct species; they are merged into one species. Default: treat subspecies as distinct species. 

	--print_dist_matrix

Distance matrix will be printed into an output file named "distance_matrix_OUTPUT_FILE". Warning: because there are n*(n-1) comparisons, the file can be very large! Default: no distance matrix file output.
 

## SYSTEM REQUIREMENTS

DiStats was developed and tested with Perl version 5.20.2. Downward compatibility has not been tested.

If you run Linux or some flavor of UNIX, your Perl version should be fine. If you work on a Windows machine, possible Perl distributions to download are:

* [Strawberry Perl](http://strawberryperl.com/)
* [ActiveState Perl](http://www.activestate.com/activeperl)
* [DWIM Perl](http://www.dwimperl.com/windows.html)

If you want to use parallel threads, you need the module Parallel::ForkManager (available from CPAN).

If you get an error like "Can't locate Some/Module.pm in @INC", you need to update your Perl installation; see below for more information.

### Required Perl modules

DiStats uses core modules that are present in any standard Perl distribution. The modules listed below have to be installed additionally. If you get the error message "Error: Requested to multi-thread but Parallel::ForkManager not usable", it means that Perl can't find the Parallel::ForkManager module. 

* Parallel::ForkManager (if you want to run the program with more than 1 thread)

## DATA REQUIREMENTS

DiStats only needs a regular Fasta file as input file. The header of the Fasta file must contain the genus and species name separated by an underscore and, after a pipe, an unambiguously unique ID for that sequence. Here is an example: 

	>Drosophila_melanogaster|Sequence0001
	ctctatatttaatatttggggtttggtcagctataatagggactgctataagagtattaa
	ttcgaatagaattaggaaatcctgggagattgttaggagatgatcatttatataatgtta
	tagttactgctcatgcttttgtaataattttttttatagtaataccaattcttattggag
	>Drosophila_yakuba|Sequence0002
	ctttatatttaatatttggggtttggtcagctataatagggactgctataagagtattaa
	ttcgaatagaattaggaaatcctgggagattgttacgagatgatcatttatataatgtta
	tagttactgctcatgcttttgtaataattaattttatagtaataccaattcttattggag


## AUTHORS

Written by Hannah Petersen and Malte Petersen, Zoological Research Museum Alexander Koenig, Bonn, Germany.

Cite as: Astrin, Höfer, Spelda, Holstein, Bayer, Hendrich, Huber, Kielhorn, Krammer, Lemke, Monje, Morinière, Rulik, Petersen, Janssen, Muster (2016). Towards a DNA Barcode Reference Database for Spiders and Harvestmen of Germany.  PLoS ONE 11(9): e0162624.  DOI:10.1371/journal.pone.0162624

## REPORTING BUGS

Report bugs to h.petersen.zfmk at uni-bonn.de or mptrsen at uni-bonn.de

DiStats home page: <http://github.com/mptrsen/distats/>

## COPYRIGHT

Copyright 2016 Hannah Petersen & Malte Petersen
