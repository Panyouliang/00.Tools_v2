#!/usr/bin/perl -w
#Zhang Pei: zhangpei@genomics.cn
use strict;
if(!$ARGV[0]){die "perl $0 junctions.bed1 junctions.bed2 ..... >all.junctions.bed\n"}
my%junctions;
for(my$i=0;$i<=$#ARGV;$i++){
	open IN,"$ARGV[$i]";
	while(<IN>){
#		print "$_\n";
		chomp;
		if(/track/){next}
		my@A=split(/\t/);
		my@B=split(/,/,$A[10]);
		my$start=$A[1]+$B[0];
		my$end=$A[2]-$B[1]+1;
		if(!$junctions{$A[0]}{$start}{$end}{$A[5]}){$junctions{$A[0]}{$start}{$end}{$A[5]}=0}
		$junctions{$A[0]}{$start}{$end}{$A[5]}+=$A[4];
	}
	close IN;
}
#print "step\n";
foreach my$key(keys %junctions){
	my@sort=sort{$a<=>$b}keys%{$junctions{$key}};
	for(my$i=0;$i<=$#sort;$i++){
		my@sort2=sort{$a<=>$b}keys%{$junctions{$key}{$sort[$i]}};
		for(my$j=0;$j<=$#sort2;$j++){
			foreach my$strand(keys %{$junctions{$key}{$sort[$i]}{$sort2[$j]}}){
				print "$key\t$sort[$i]\t$sort2[$j]\t$strand\t$junctions{$key}{$sort[$i]}{$sort2[$j]}{$strand}\n";
			}
		}
	}
}
