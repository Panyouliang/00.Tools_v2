#!/usr/bin/perl -w
use strict;
if(!$ARGV[0]){die "perl $0 chr_list genome.fa alignment.bam> junc.bed\n"}
my%list;
open IN,"$ARGV[0]";
while(<IN>){
	chomp;
	$list{$_}++;
}
close IN;
open IN,"$ARGV[1]";
our%seq;
my$chr;
my$seq="";
while(<IN>){
	chomp;
	my@A=split(/\s+/);
	if($A[0]=~/>(.+?)$/){
		if($seq && $list{$chr}){
			$seq{$chr}=$seq
		}
		$seq="";
		$chr=$1
	}else{
		$seq=$seq.$_;
	}
}
if($list{$chr}){
	$seq{$chr}=$seq;
}
close IN;
my%junc;
foreach $chr(keys %list){
	system "/hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/zhangpei/bin/samtools-0.1.18/samtools view $ARGV[2] $chr >$chr.sam";
	open IN,"$chr.sam";
	while(<IN>){
		chomp;
		my@A=split(/\t/);
#		if($A[3] ==1 && $0 !~ /XS:A/){next}
		my$pair=$A[1] & 2;
		if($pair==0){next}
		my$tag=0;
		foreach my$ele(@A){
			if($ele eq 'NH:i:1'){$tag=1}
		}
		if($tag == 0){next}
		if($A[5] !~ /N/){next}
		if($A[5] =~ /I/){next}
		if($A[5] =~ /D/){next}
		my$seq=&MiSeq($A[5],$A[3],$A[2]);
		my@value=&TagCheck($A[5],$seq,$A[9]);
		if($value[0] eq 'NA'){next}
		my$length=length($A[9]);
		for(my$i=0;$i<=$#value;$i+=2){
			my$start=$A[3]+$value[$i]-1;
			my$end=$A[3]+$value[$i+1]-1;
			my$startChr=substr($seq{$A[2]},$start,2);
			$startChr=uc($startChr);
			my$endChr=substr($seq{$A[2]},$end-3,2);
			$endChr=uc($endChr);
#			print "$startChr\t$endChr\n";
			my$strand;
			if($startChr eq 'GT' && $endChr eq 'AG'){
				$strand = "+"
			}elsif($startChr eq 'GC' && $endChr eq 'AG'){
				$strand = "+"
			}elsif($startChr eq 'AT' && $endChr eq 'AC'){
				$strand = "+"
			}elsif($startChr eq 'CT' && $endChr eq 'AC'){
				$strand = "-"
			}elsif($startChr eq 'CT' && $endChr eq 'GC'){
				$strand = "-"
			}elsif($startChr eq 'GT' && $endChr eq 'AT'){
				$strand = "-"
			}else{$strand="."}
			if(!$junc{$A[2]}{$start}{$end}{$strand}[0]){$junc{$A[2]}{$start}{$end}{$strand}[0]=$value[$i]}
			if(!$junc{$A[2]}{$start}{$end}{$strand}[1]){$junc{$A[2]}{$start}{$end}{$strand}[1]=$length-$value[$i]}
			if($value[$i]>$junc{$A[2]}{$start}{$end}{$strand}[0]){$junc{$A[2]}{$start}{$end}{$strand}[0]=$value[$i]}
			if($length-$value[$i]>$junc{$A[2]}{$start}{$end}{$strand}[1]){$junc{$A[2]}{$start}{$end}{$strand}[1]=$length-$value[$i]}
			$junc{$A[2]}{$start}{$end}{$strand}[2]++;
		}
	}
	close IN;
	system "rm $chr.sam";
}
my$JuncNum=0;
foreach my$chr(keys %junc){
	my@start=sort{$a<=>$b}keys%{$junc{$chr}};
	foreach my$ele(@start){
		my@end=sort{$a<=>$b}keys%{$junc{$chr}{$ele}};
		foreach my$ele2(@end){
			foreach my$strand(keys %{$junc{$chr}{$ele}{$ele2}}){
				$JuncNum++;
				my$start=$ele-$junc{$chr}{$ele}{$ele2}{$strand}[0];
				my$end=$ele2+$junc{$chr}{$ele}{$ele2}{$strand}[1]-1;
				my$intron=$ele2-$start;
#				print "$chr\t$ele\t$ele2\t$strand\t$junc{$chr}{$ele}{$ele2}{$strand}[2]\n";
				print "$chr\t$start\t$end\tJunc$JuncNum\t$junc{$chr}{$ele}{$ele2}{$strand}[2]\t$strand\t$start\t$end\t255,0,0\t2\t$junc{$chr}{$ele}{$ele2}{$strand}[0],$junc{$chr}{$ele}{$ele2}{$strand}[1]\t0,$intron\n";
			}
		}
	}
}
sub MiSeq{
	my@tag=split(/([A-Z-])/,$_[0]);
	my$loc=$_[1]-1;
	my$seq="";
	for(my$i=0;$i<=$#tag;$i++){
		if($tag[$i] eq 'M'){
#			die "$loc\t$tag[$i-1]\n";
			$seq=$seq.substr($seq{$_[2]},$loc,$tag[$i-1]);
			$loc+=$tag[$i-1];
		}
		if($tag[$i] eq "D"){$loc+=$tag[$i-1]}
		if($tag[$i] eq 'N'){$loc+=$tag[$i-1]}
		if($tag[$i] eq 'I'){
			for(my$j=1;$j<=$tag[$i-1];$j++){
				$seq=$seq.'I';
			}
		}
	}
	return($seq);
}
sub TagCheck{
	my%mis;
	my@seq1=split(//,uc($_[1]));
	my@seq2=split(//,uc($_[2]));
#	die "@seq1\n@seq2\n";
	for(my$i=0;$i<=$#seq1;$i++){
		if($seq1[$i] ne $seq2[$i]){
			$mis{$i+1}++;
#			print $i+1 ."\n"
		}
	}
#	die "\n";
	my@tag=split(/([A-Z-])/,$_[0]);
#	print "@tag\n";
	my$arg=1;
	my@loc;
	my$loc=0;
	for(my$i=0;$i<=$#tag;$i++){
		if($tag[$i] eq 'M'){$loc+=$tag[$i-1]}
		if($tag[$i] eq 'D'){$loc+=$tag[$i-1]}
		if($tag[$i] eq 'N'){
			if($tag[$i-2] ne 'M' || $tag[$i+2] ne 'M'){$arg=0; last}
			if($tag[$i-3] <5 || $tag[$i+1] <5){$arg=0; last}
			push @loc,($loc, $loc+$tag[$i-1]+1);
			foreach my$ele($loc-4..$loc){
				if($mis{$ele}){$arg=0;last}
			}
			$loc=$loc+1;
			foreach my$ele($loc..$loc+4){
				if($mis{$ele}){$arg=0;last}
			}
			$loc+=$tag[$i-1]-1;
		}
	}
	if($arg==0){$loc[0]="NA"}
	return(@loc);
}
