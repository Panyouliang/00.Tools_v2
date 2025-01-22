#!/usr/bin/perl -w
use strict;
die "Usage: <genewise result> <query pep>\n" unless @ARGV == 2;

my %pepLen;
open IN, $ARGV[1];
$/ = ">";
<IN>;
while (<IN>) {
	/(.+)\n/;
	my $id = (split /\s+/, $1)[0];
	s/.+\n//;
	s/\s+|>//g;
	my $len = length($_);
	$pepLen{$id} = $len;
}
$/ = "\n";
close IN;

my %result;
open IN, $ARGV[0];
$/ = "//\ngenewise";
while (<IN>) {
	/Query\s+protein:\s+(.+?)\n.+Target\s+Sequence\s+(.+?)\nStrand:\s+(.+?)\n/s;
	my ($query, $target, $strand) = ($1, $2, $3);
	#ENSMMUP00000022738_chr10_33962498_33966016
	#my ($chr, $chr_bg) = (split /_/, $target)[1,2];
	my @chr_tmp = split /##/, $target;
	my ($chr, $chr_bg) = ($chr_tmp[-3], $chr_tmp[-2]);
	if ($strand eq "forward") {
		$strand = "+";
	} elsif ($strand eq "reverse") {
		$strand = "-";
	} else {
		die "strand $strand cannot be parse";
	}

	## get the location of each exon 
	my @info = split /\n/;
	my @exon_location;
	my @pep_align_block;
	for (my $i = 0; $i < @info; $i ++) {
		if ($info[$i] =~ /\s+Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/) {
			my ($exon_bg, $exon_ed, $phase) = ($1, $2, $3);
			$exon_bg = $exon_bg + $chr_bg - 1;
			$exon_ed = $exon_ed + $chr_bg - 1;
			($exon_bg, $exon_ed) = sort {$a <=> $b} ($exon_bg, $exon_ed);
			push @exon_location, [$exon_bg, $exon_ed, $phase];
		} elsif ($info[$i] =~ /\s+Supporting\s+\d+\s+\d+\s+(\d+)\s+(\d+)/) {
			push @pep_align_block, ($1, $2);
			
		}
	}

	## get alignRate
	my ($pep_align_bg, $pep_align_ed) = ($pep_align_block[0], $pep_align_block[-1]);
	my $alignRate = sprintf "%.4f", ($pep_align_ed - $pep_align_bg + 1) / $pepLen{$query};
	
	## sort the exon by bg site
	my %p_to_bg;
	foreach my $p (@exon_location) {
		my $bg = $p->[0];
		$p_to_bg{$p} = $bg;
	}
	@exon_location = sort {$p_to_bg{$a} <=> $p_to_bg{$b}} @exon_location;

	## get the bg and ed site of the gene
	my ($gene_bg, $gene_ed) = ($exon_location[0]->[0], $exon_location[-1]->[1]);

	## output gff format
	my $out;
	#print "$chr\tprotein_coding\tmRNA\t$gene_bg\t$gene_ed\t.\t$strand\t.\tID=$query;\n";
	$out .= "$chr\tGenewise\tmRNA\t$gene_bg\t$gene_ed\t.\t$strand\t.\tID=$query;alignRate=$alignRate;\n";
	foreach my $p (@exon_location) {
		my ($exon_bg, $exon_ed, $phase) = @$p;
		#print "$chr\tprotein_coding\tCDS\t$exon_bg\t$exon_ed\t.\t$strand\t$phase\tParent=$query;\n";
		$out .= "$chr\tGenewise\tCDS\t$exon_bg\t$exon_ed\t.\t$strand\t$phase\tParent=$query;\n";
	}
	$result{$chr}{$gene_bg}{$query} = $out;
}
$/ = "\n";
close IN;

foreach my $chr (sort keys %result) {
	foreach my $bg (sort {$a <=> $b} keys %{$result{$chr}}) {
		foreach my $gene (keys %{$result{$chr}{$bg}}) {
			print $result{$chr}{$bg}{$gene};
		}
	}
}
