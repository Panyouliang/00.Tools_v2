#!/usr/bin/perl -w
#Zhang Pei: zhangpei@genomics.org.cn
use strict;
use Getopt::Long;
GetOptions(
		'GffFile=s'=>\our$GffFile,
		'JuncFile=s'=>\our$JuncFile,
		'BamFile=s'=>\our$BamFile,
		'GenomeFile=s'=>\our$GenomeFile,
		'ResultFile=s'=>\our$ResultFile,
		'Extension=i'=>\our$Extension,
		'Unique'=>\our$Unique,
		'StrandSpecific'=>\our$StrandSpecific,
		'JuncDepth=i'=>\our$JuncDepth,
		'CovDepth=i'=>\our$CovDepth
		);
my$help=<<"HELP";
	--GffFile:	<str>	Gff format gene annotation file input
	--JuncFile:	<str>	Bed format junction annotation file input
	--BamFile:	<str>	Transcriptome alignment in indexed Bam format
	--GenomeFile:	<str>	Genome file input in fasta format
	--ResultFile:	<str>	Output result in Gff format
	--Extension:	<int>	Upstream and Downstream extension length for every gene; default: 50000
	--Unique:	Only uniquely mapped reads used
	--StrandSpecific:	StrandSpecific Sequence
	--JuncDepth:	<int> Junctions reads threshold
	--CovDepth:	<int> CovDepth
HELP
if(!$GffFile || !$BamFile || !$ResultFile){die "$help\n"}
if(!$GenomeFile && !$JuncFile){die "$help\n"}
if(!$Extension){$Extension=50000}
if(!$JuncDepth){$JuncDepth=5}
if(!$CovDepth){$CovDepth=5}
use lib '/hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/zhangpei/perl/pm';
use Overlap;
my%start;
my%end;
my%strand;
my%loc;
our%seq;
if(!$JuncFile){
	my$chr;
	my$seq="";
	if($GenomeFile=~/\.gz$/){
		open IN, "gzip -cd $GenomeFile |"
	}else{
		open IN,"$GenomeFile";
	}
	while(<IN>){
		chomp;
		my@A=split(/\s+/);
		if($A[0]=~/>(.+?)$/){
			if($seq){$seq{$chr}=$seq}
			$seq="";
			$chr=$1
		}else{
			$seq=$seq.$_;
		}
	}
	$seq{$chr}=$seq;
	close IN;
	undef($chr);
	undef($seq);
}
open IN,"$GffFile";
while(<IN>){
	my@A=split(/\s+/);
#	if($A[2] eq 'mRNA' || $A[2] eq 'transcript'){
#		my$id;
#		if($A[8]=~/ID=(.+?);/){$id=$1}
#		$loc{$A[0]}{$id}[0]=$A[3];
#		$loc{$A[0]}{$id}[1]=$A[4];
#	}
	if($A[2] ne 'CDS'){next}
	my$id;
	if($A[8]=~/Parent=(.+?);/){$id=$1}
	push@{$start{$A[0]}{$id}},$A[3];
	push@{$end{$A[0]}{$id}},$A[4];
	$strand{$id}=$A[6];
}
close IN;
foreach my$key(keys %start){
	foreach my$ke(keys %{$start{$key}}){
		my@sort=sort{$a<=>$b}@{$start{$key}{$ke}};
		@{$start{$key}{$ke}}=@sort;
		$loc{$key}{$ke}[0]=$sort[0];
		@sort=sort{$a<=>$b}@{$end{$key}{$ke}};
		$loc{$key}{$ke}[1]=$sort[-1];
		@{$end{$key}{$ke}}=@sort;
	}
}
my%temp_junc;
my%junctions;
my%junctions_num;
if($JuncFile){
	open IN,"$JuncFile";
	while(<IN>){
		if(/track/){next}
		chomp;
		my@A=split(/\t/);
		if($A[4]<$JuncDepth){next}
		my@B=split(/,/,$A[10]);
		my$start=$A[1]+$B[0];
		my$end=$A[2]-$B[1]+1;
		$temp_junc{$A[0]}{$start}{$end}{$A[5]}=$A[4];
	}
	close IN;
	foreach my$key(keys %temp_junc){
		if(!$junctions_num{$key}){$junctions_num{$key}=0}
		my@start=sort{$a<=>$b}keys%{$temp_junc{$key}};
		foreach my$start(@start){
			my@end=sort{$a<=>$b}keys%{$temp_junc{$key}{$start}};
			foreach my$end(@end){
				foreach my$strand(keys %{$temp_junc{$key}{$start}{$end}}){
					$junctions{$key}[$junctions_num{$key}][0]=$start;
					$junctions{$key}[$junctions_num{$key}][1]=$end;
					$junctions{$key}[$junctions_num{$key}][2]=$strand;
					$junctions{$key}[$junctions_num{$key}][3]=$temp_junc{$key}{$start}{$end}{$strand};
					$junctions_num{$key}++;
				}
			}
		}
	}
}
open OUT, ">$ResultFile";
foreach my$key(keys %start){
	foreach my$ke(keys %{$start{$key}}){
		my$start=$start{$key}{$ke}[0];
		my$end=$end{$key}{$ke}[-1];
		my$temp_start;
		my$temp_end;
		if($start-$Extension>0){$temp_start=$start-$Extension}else{$temp_start=1}
		$temp_end=$end+$Extension;
		my%coverage_depth;
		my%coverage_end;
		&sam($key,$temp_start,$temp_end,$strand{$ke});
		&coverage(\%coverage_depth,\%coverage_end);
		if(!$JuncFile){
			undef(%junctions);
			undef(%junctions_num);
			&Junction(\%junctions,\%junctions_num);
		}
=cut	for(my$i=0;$i<=$#{$start{$key}{$ke}};$i++){
			$length=$length+$end{$key}{$ke}[$i]-$start{$key}{$ke}[$i]+1;
			foreach my$k(keys %{$coverage_depth{$key}}){
				my@overlap=&overlap($start{$key}{$ke}[$i],$end{$key}{$ke}[$i],$k,$coverage_end{$key}{$k});
				if($overlap[0]>0){
					$totalbase=$totalbase+$overlap[0]*$coverage_depth{$key}{$k};
				}
			}
		}
		my$average=$totalbase/$length;
=cut
		my@start;
		my@end;
		my@sort=sort{$a<=>$b}keys%{$coverage_depth{$key}};
		my$arg=0;
		foreach my$i(0..$#sort){
			if($coverage_depth{$key}{$sort[$i]}<$CovDepth){
				if($arg==1){push @end,$coverage_end{$key}{$sort[$i-1]}}
				$arg=0;
			}else{
				if($arg==0){push @start,$sort[$i]}
				$arg=1;
			}
		}
		if($arg==1){push @end,$coverage_end{$key}{$sort[-1]}}
		my@over;
		my$start2=$start;
		my$end2=$end;
		for(my$i=0;$i<=$#start;$i++){
			my@overlap=&overlap($start,$end,$start[$i],$end[$i]);
			if($overlap[0]>0){
				push @over,$i;
				if($start[$i]<$start2){$start2=$start[$i]}
				if($end[$i]>$end2){$end2=$end[$i]}
			}
		}
		my$juncstart;
		my$juncend;
		if(!exists $junctions{$key}){
			print OUT "$key\tTrans\tmRNA\t$loc{$key}{$ke}[0]\t$loc{$key}{$ke}[1]\t.\t$strand{$ke}\t.\tID=$ke;\n";
			for(my$i=0;$i<=$#{$start{$key}{$ke}};$i++){
			    print OUT "$key\tTrans\tCDS\t$start{$key}{$ke}[$i]\t$end{$key}{$ke}[$i]\t.\t$strand{$ke}\t.\tParent=$ke;\n";
			}
			next;
		}
		for(my$i=0;$i<=$#{$junctions{$key}};$i++){
			if($junctions{$key}[$i][0]<=$temp_start){$juncstart=$i}
			if($junctions{$key}[$i][1]<=$temp_end){$juncend=$i}
			if($junctions{$key}[$i][0] >$temp_end){last}
		}
		if(!$juncstart){$juncstart=0}
		if(!$juncend){$juncend=0}
		for(my$i=$juncend;$i>=$juncstart;$i--){
			if($junctions{$key}[$i][2] ne $strand{$ke}){next}
			if($junctions{$key}[$i][1]>=$start2-2 && $junctions{$key}[$i][1]<=$end2+2){
				for(my$j=0;$j<=$#start;$j++){
					if($junctions{$key}[$i][0]>=$start[$j]-2 && $junctions{$key}[$i][0]<=$end[$j]+2){
						if($start2>$start[$j]){$start2=$start[$j]}
					}
				}
			}
		}
		for(my$i=$juncstart;$i<=$juncend;$i++){
			if($junctions{$key}[$i][2] ne $strand{$ke}){next}
			if($junctions{$key}[$i][0]>=$start2-2 && $junctions{$key}[$i][0]<=$end2+2){
				for(my$j=0;$j<=$#start;$j++){
					if($junctions{$key}[$i][1]>=$start[$j]-2 && $junctions{$key}[$i][1]<=$end[$j]+2){
						if($end2<$end[$j]){$end2=$end[$j]}
					}
				}
			}
		}
		my@utr5start;
		my@utr5end;
		my@utr3start;
		my@utr3end;	
		my%exon;
		my$tempstart;
		my$score;
		push@{$exon{start}},$start2;
		for(my$i=$juncstart;$i<=$juncend;$i++){
			if($junctions{$key}[$i][2] ne $strand{$ke}){next}
			if($junctions{$key}[$i][0]<=$start2 || $junctions{$key}[$i][1]>=$end2){next}
			if(!$exon{end}[-1]){
				push@{$exon{end}},$junctions{$key}[$i][0];
				push@{$exon{start}},$junctions{$key}[$i][1];$score=$junctions{$key}[$i][3]
			}
			if($exon{end}[-1] && $junctions{$key}[$i][0]<=$exon{start}[-1]){
				if($junctions{$key}[$i][3]>$score){
					$exon{end}[-1]=$junctions{$key}[$i][0];
					$exon{start}[-1]=$junctions{$key}[$i][1];
					$score=$junctions{$key}[$i][3];
				}
			}
			if($exon{end}[-1] && $junctions{$key}[$i][0]>$exon{start}[-1]){
				push@{$exon{end}},$junctions{$key}[$i][0];
				push@{$exon{start}},$junctions{$key}[$i][1];
				$score=$junctions{$key}[$i][3];
			}
		}
		push@{$exon{end}},$end2;
		for(my$i=0;$i<=$#{$exon{start}};$i++){
			if($exon{start}[$i]>$end2 || $exon{end}[$i]<$start2){next}
			if($exon{end}[$i]<$start){
				push @utr5start,$exon{start}[$i];
				push @utr5end,$exon{end}[$i];
			}
			if($exon{start}[$i]<$start && $exon{end}[$i]>=$start){
				my$arg=$start-1;
				push @utr5start,$exon{start}[$i];
				push @utr5end,$arg;
			}
			if($exon{start}[$i]<=$end && $exon{end}[$i]>$end){
				my$arg=$end+1;
				push @utr3start,$arg;
				push @utr3end,$exon{end}[$i];
			}
			if($exon{start}[$i]>$end){
				push @utr3start,$exon{start}[$i];
				push @utr3end,$exon{end}[$i];
			}
		}
		my$overlap;
		my$shallow;
		my$shallow_loc;
		if($strand{$ke} eq '+' && $utr5start[0]){
			for(my$i=0;$i<=$#sort;$i++){
				my@overlap=&overlap($sort[$i],$coverage_end{$key}{$sort[$i]},$utr5start[0],$utr5end[0]);
				if($overlap[0]>0){
					$overlap=$i;
					last;
				}
			}
			$shallow=$coverage_depth{$key}{$sort[$overlap]};
			$shallow_loc=$utr5start[0];
			for(my$i=$overlap;$i>=0;$i--){
				if($coverage_depth{$key}{$sort[$i]}<$CovDepth){last}
				if($coverage_depth{$key}{$sort[$i]}<$shallow){
					$shallow_loc=$sort[$i];
					$shallow=$coverage_depth{$key}{$sort[$i]};
				}
			}
			$utr5start[0]=$shallow_loc;
			$start2=$shallow_loc;
		}
		if($strand{$ke} eq '-' && $utr3start[0]){
			for(my$i=$#sort;$i>=0;$i--){
				my@overlap=&overlap($sort[$i],$coverage_end{$key}{$sort[$i]},$utr3start[-1],$utr3end[-1]);
				if($overlap[0]>0){
					$overlap=$i;
					last;
				}
			}
			$shallow=$coverage_depth{$key}{$sort[$overlap]};
			$shallow_loc=$utr3end[-1];
			for(my$i=$overlap;$i<=$#sort;$i++){
				if($coverage_depth{$key}{$sort[$i]}<$CovDepth){last}
				if($coverage_depth{$key}{$sort[$i]}<$shallow){
					$shallow_loc=$coverage_end{$key}{$sort[$i]};
					$shallow=$coverage_depth{$key}{$sort[$i]};
				}
			}
			$utr3end[-1]=$shallow_loc;
			$end2=$shallow_loc;
		}
		my$over=0;
		foreach my$k(keys %{$loc{$key}}){
			if($k eq $ke){next}
			my@overlap=&overlap($loc{$key}{$k}[0],$loc{$key}{$k}[1],$start2,$end2);
			if($overlap[0]>0){$over++}
		}
		if($over>0){
			print OUT "$key\tTrans\tmRNA\t$loc{$key}{$ke}[0]\t$loc{$key}{$ke}[1]\t.\t$strand{$ke}\t.\tID=$ke;\n";
			for(my$i=0;$i<=$#{$start{$key}{$ke}};$i++){
				print OUT "$key\tTrans\tCDS\t$start{$key}{$ke}[$i]\t$end{$key}{$ke}[$i]\t.\t$strand{$ke}\t.\tParent=$ke;\n";
			}
			next;
		}
		$loc{$key}{$ke}[0]=$start2;
		$loc{$key}{$ke}[1]=$end2;
		my$arg1;
		my$arg2;
		if($strand{$ke} eq '+'){$arg1='UTR5';$arg2='UTR3'}else{$arg1='UTR3',$arg2='UTR5'}
		print OUT "$key\tGLEAN\tmRNA\t$start2\t$end2\t.\t$strand{$ke}\t.\tID=$ke;\n";
		for(my$i=0;$i<=$#utr5start;$i++){
			print OUT "$key\tTRANS\t$arg1\t$utr5start[$i]\t$utr5end[$i]\t.\t$strand{$ke}\t.\tParent=$ke;\n";
		}
		for(my$i=0;$i<=$#{$start{$key}{$ke}};$i++){
			print OUT "$key\tGLEAN\tCDS\t$start{$key}{$ke}[$i]\t$end{$key}{$ke}[$i]\t.\t$strand{$ke}\t.\tParent=$ke;\n";
		}
		for(my$i=0;$i<=$#utr3start;$i++){
			print OUT "$key\tTRANS\t$arg2\t$utr3start[$i]\t$utr3end[$i]\t.\t$strand{$ke}\t.\tParent=$ke;\n";
		}
	}
}
sub sam{
	my$scaffold=shift;
	my$start=shift;
	my$end=shift;
	my$strand=shift;
	system "/hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/zhangpei/bin/samtools-0.1.18/samtools view $BamFile $scaffold:$start-$end >temp.sam";
	if($Unique){
		system "perl /hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/zhangpei/bin/transcriptome/uniq_filter_sam.pl temp.sam temp.uniq.sam";
		system "mv temp.uniq.sam temp.sam";
	}
	if($StrandSpecific){
		open SUBIN,"temp.sam";
		open SUBOUT, ">temp.filter.sam";
		while(<SUBIN>){
			chomp;
			my@A=split(/\t/);
			my$REVERSE=$A[1] & 16;
			my$FQ=$A[1] & 128;
			if($strand eq '+'){
				if($REVERSE == 16 && $FQ == 0){print SUBOUT "$_\n"}
				if($REVERSE == 0 && $FQ == 128){print SUBOUT "$_\n"}
			}
			if($strand eq '-'){
				if($REVERSE ==16 && $FQ == 128){print SUBOUT "$_\n"}
				if($REVERSE ==0 && $FQ == 0){print SUBOUT "$_\n"}
			}
		}
		close SUBIN;
		close SUBOUT;
		system "mv temp.filter.sam temp.sam";
	}
}
sub coverage{
	my$coverage_depth=shift;
	my$coverage_end=shift;
	system "/hwfssz1/ST_EARTH/Reference/ST_DIVERSITY/PUB/USER/zhangpei/bin/tophat-1.3.1/bin/wiggles temp.sam temp.wiggle";
	open SUBIN,"temp.wiggle";
	while(<SUBIN>){
		chomp;
		if(/track/){next}
		my@A=split(/\t/);
		my$start_wig=$A[1]+1;
		$$coverage_depth{$A[0]}{$start_wig}=$A[3];
		$$coverage_end{$A[0]}{$start_wig}=$A[2];
	}
	close SUBIN;
	system "rm temp.wiggle";
}
sub Junction{
	my$junctions=shift;
	my$junctions_num=shift;
	my%junc;
	open SUBIN,"temp.sam";
	while(<SUBIN>){
		chomp;
		my@A=split(/\t/);
		if($A[5] !~ /N/){next}
		if($A[5] =~ /I/){next}
		if($A[5] =~ /D/){next}
		if($A[5] =~ /S/){next}
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
#       	print "$startChr\t$endChr\n";
			my$strand;
			if($startChr eq 'GT' && $endChr eq 'AG'){
				$strand = "+"
			}elsif($startChr eq 'CT' && $endChr eq 'AC'){
				$strand = "-"
			}else{next}
			if(!$junc{$A[2]}{$start}{$end}{$strand}[0]){$junc{$A[2]}{$start}{$end}{$strand}[0]=$value[$i]}
			if(!$junc{$A[2]}{$start}{$end}{$strand}[1]){$junc{$A[2]}{$start}{$end}{$strand}[1]=$length-$value[$i]}
			if($value[$i]>$junc{$A[2]}{$start}{$end}{$strand}[0]){$junc{$A[2]}{$start}{$end}{$strand}[0]=$value[$i]}
			if($length-$value[$i]>$junc{$A[2]}{$start}{$end}{$strand}[1]){$junc{$A[2]}{$start}{$end}{$strand}[1]=$length-$value[$i]}
			$junc{$A[2]}{$start}{$end}{$strand}[2]++;
		}
	}
	close SUBIN;
	foreach my$chr(keys %junc){
		my@start=sort{$a<=>$b}keys%{$junc{$chr}};
		foreach my$ele(@start){
			my@end=sort{$a<=>$b}keys%{$junc{$chr}{$ele}};
			foreach my$ele2(@end){
				foreach my$strand(keys %{$junc{$chr}{$ele}{$ele2}}){
					if($junc{$chr}{$ele}{$ele2}{$strand}[2]<$JuncDepth){next}
					if(!$$junctions_num{$chr}){$$junctions_num{$chr}=0}
					$$junctions{$chr}[$junctions_num{$chr}][0]=$ele;
					$$junctions{$chr}[$junctions_num{$chr}][1]=$ele2;
					$$junctions{$chr}[$junctions_num{$chr}][2]=$strand;
					$$junctions{$chr}[$junctions_num{$chr}][3]=$junc{$chr}{$ele}{$ele2}{$strand}[2];
					$$junctions_num{$chr}++;
				}
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
#           die "$loc\t$tag[$i-1]\n";
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
#   die "@seq1\n@seq2\n";
	for(my$i=0;$i<=$#seq1;$i++){
		if($seq1[$i] ne $seq2[$i]){
			$mis{$i+1}++;
#           print $i+1 ."\n"
		}
	}
#   die "\n";
	my@tag=split(/([A-Z-])/,$_[0]);
#   print "@tag\n";
	my$arg=1;
	my@loc;
	my$loc=0;
	for(my$i=0;$i<=$#tag;$i++){
		if($tag[$i] eq 'M'){$loc+=$tag[$i-1]}
		if($tag[$i] eq 'D'){$loc+=$tag[$i-1]}
		if($tag[$i] eq 'N'){
			if($tag[$i-2] ne 'M' || $tag[$i+2] ne 'M'){$arg=0; last}
			if($tag[$i-3] <10 || $tag[$i+1] <10){$arg=0; last}
			push @loc,($loc, $loc+$tag[$i-1]+1);
			foreach my$ele($loc-9..$loc){
				if($mis{$ele}){$arg=0;last}
			}
			$loc=$loc+1;
			foreach my$ele($loc..$loc+9){
				if($mis{$ele}){$arg=0;last}
			}
			$loc+=$tag[$i-1]-1;
		}
	}
	if($arg==0){$loc[0]="NA"}
	return(@loc);
}
