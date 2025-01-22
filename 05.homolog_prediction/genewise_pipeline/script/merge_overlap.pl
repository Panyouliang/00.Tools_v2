#!/usr/bin/perl -w
use strict;
die "Usage: <alg file> <overlap cutoff(0~1)>\n" unless @ARGV == 2;
my $overlapCutoff = $ARGV[1];

## get and sort the data by chr and start site
my @data;
open IN, $ARGV[0];
while (<IN>) {
	next if /^#/;
	push @data, $_;
}
#my @data = <IN>;
close IN;

my %line_to_chr;
my %line_to_bg;
foreach (@data) {
	my ($chr, $bg) = (split /\s+/)[5, 6];
	$line_to_chr{$_} = $chr;
	$line_to_bg{$_} = $bg;
}
@data = sort {$line_to_chr{$a} cmp $line_to_chr{$b} or $line_to_bg{$a} <=> $line_to_bg{$b}} @data;


## merge overlaping genes. First by genewise score than by identity.
for (my $i = 0; $i < @data-1; $i ++) {
	my @info1 = split /\s+/, $data[$i];
	my @info2 = split /\s+/, $data[$i+1];
	my ($gene1, $chr1, $strand1, $bg1, $ed1, $score1, $id1) = ($info1[0], $info1[5], $info1[4], $info1[6], $info1[7], $info1[9], $info1[10]);
	my ($gene2, $chr2, $strand2, $bg2, $ed2, $score2, $id2) = ($info2[0], $info2[5], $info2[4], $info2[6], $info2[7], $info2[9], $info2[10]);
	if ($chr1 eq $chr2) {
		if ($ed2 >= $bg1 && $bg2 <= $ed1) {
			if ($strand1 eq $strand2 || $strand1 ne $strand2) {
				my ($overlap_bg, $overlap_ed) = (sort {$a <=> $b} ($bg1, $ed1, $bg2, $ed2))[1,2];
				my $overlap = $overlap_ed - $overlap_bg + 1;
				my $gene1_len = $ed1 - $bg1 + 1;
				my $gene2_len = $ed2 - $bg2 + 1;
				my $coverRate1 = $overlap / $gene1_len;
				my $coverRate2 = $overlap / $gene2_len;
				if ($coverRate1 > $overlapCutoff || $coverRate2 > $overlapCutoff) {
					#print $data[$i], $data[$i+1] , "$overlap_bg\t$overlap_ed\t$overlap\t$coverRate1\t$coverRate2\n\n";
					if ($score1 > $score2) {
						$data[$i+1] = $data[$i]; 
						$data[$i] = "";
					} elsif ($score2 > $score1) {
						$data[$i] = "";
					} elsif ($score1 == $score2) {
						if ($id1 > $id2) {
							$data[$i+1] = $data[$i]; 
							$data[$i] = "";
						} elsif ($id1 <= $id2) {
							$data[$i] = "";
						}
					}
				}
			}
		}
	}
}
## output resutl
foreach (@data) {
	next unless $_;
	print;
}
