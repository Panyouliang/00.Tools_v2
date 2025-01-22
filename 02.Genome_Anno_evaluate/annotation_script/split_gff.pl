#!/usr/bin/perl -w
use strict;
my$num=$ARGV[0];
my$file=$ARGV[1];
open IN,"$file";
my$tag;
my%gff;
while(<IN>){
	chomp;
	my@A=split(/\t/);
	if($A[2] eq 'mRNA'){$tag=$_}
	if($A[2] eq 'CDS'){push @{$gff{$tag}},$_}
}
close IN;
$tag=0;
foreach my$key(keys %gff){
	if($tag % $num ==0){
		my$name=int($tag/$num);
		mkdir "part$name";
		open OUT,">part$name/part$name.gff";
	}
	print OUT "$key\n";
	foreach my$ele(@{$gff{$key}}){print OUT "$ele\n"}
	$tag++;
}
