#!/usr/bin/perl -w
use strict;
die "Usage: <in dirs> <gene list>\n" unless @ARGV >= 2;
my $list = pop;

my %mark;
open IN, $list;
while (<IN>) {
	my $gene = (split /\s+/)[0];
	$mark{$gene} ++;
}
close IN;

foreach my $in_dir (@ARGV) {
	opendir DH, $in_dir;
	while (my $file = readdir DH) {
		next unless $file =~ /^(.+)\.fa\.gw$/;
		my $gene = $1;
		next unless $mark{$gene};
		my $in_file = "$in_dir/$file";
		open IN, $in_file;
		my @data = <IN>;
		close IN;
		print @data;
	}
	closedir DH;
}
