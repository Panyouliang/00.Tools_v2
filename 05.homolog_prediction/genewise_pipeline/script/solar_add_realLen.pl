#!/usr/bin/perl -w
use strict;
die "Usage: <solar result> <query files(one or more)>\n" unless @ARGV >= 2;
my $solar_file = shift @ARGV;

my %seqLen;
foreach my $fasta_file (@ARGV) {
	if ( $fasta_file =~ /\.gz$/ ) {
        open IN, "gzip -dc $fasta_file|";
    }else{
        open IN, $fasta_file;
    }
    #open IN, $fasta_file;
	$/ = ">";
	<IN>;
	while (<IN>) {
		/(.+)\n/;
		my $id = (split /\s+/)[0];
		s/.+\n//;
		s/\s+|>//g;
		my $len = length($_);
		$seqLen{$id} = $len;	
	}
	$/ = "\n";
	close IN;
}

open IN, $solar_file;
while (<IN>) {
	my @info = split /\s+/;
	$info[1] = $seqLen{$info[0]};
	$info[6] = $seqLen{$info[5]} if $seqLen{$info[5]};
	my $out = join "\t", @info;
	print "$out\n";
}
close IN;
