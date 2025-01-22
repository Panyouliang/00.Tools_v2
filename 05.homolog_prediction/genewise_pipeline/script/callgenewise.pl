#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
my ($pepFile, $nucleotideFile, $strandList, $outDir, $job_num, $keyword);
my ($Verbose,$Help);
GetOptions(
	"pep:s"=>\$pepFile,
	"nuc:s"=>\$nucleotideFile,
	"list:s"=>\$strandList,
	"out:s"=>\$outDir,
	"key:s"=>\$keyword,
	"num:s"=>\$job_num,
	"help"=>\$Help,
	"verbose"=>\$Verbose,
);
if (!$pepFile || !$nucleotideFile || $Help) {
	print <<"Usage End."; 

Description:
	This program is used to call genewise with a fasta file of pep sequcences and a fasta file of the corresponding nuc sequences.

	Version: 2.0  Date: 2009-11-5
	Author:  liqiye <liqiye\@genomics.org.cn>
			
Usage:
	--pep	       pep sequences, fasta format only
	--nuc          nuc sequences, fasta format only
	--list         strand of each gene, optional. Format: geneID  +/-. Such as
		       ENSMMUP00000000013	-
		       If the strand of a gene is "+", -tfor will be set. If the strand of a gene is "-", -trev will be set.
		       and if --list is not given, -both will be set.
	--key          the prefix of the jobs qsub, default="".
	--num          number of job to qsub, default=10.
	--out          set the output directory, default="./"
	--verbose      output running progress information to screen
	--help         output help information to screen

Example:
	perl 
	
Usage End.
	
	exit;
}
$outDir ||= "./";
mkdir $outDir unless -e $outDir;
$job_num ||= 10;
$keyword ||= "";
#my $program = "/share/raid1/genome/bin/genewise";
#my $program = "/ifs1/pub/genome/bin/genewise";
#my $program = "/nas/RD_09A/zhangpei/bin/wise2.2.0/src/bin/genewise";
#my $program = "/hwfssz1/ST_EARTH/P18Z10200N0107/P17Z10200N0101_Metazoa_RNA_Editing/liqiye/PC_PA_UN/bin/Wise2/wise2.2.0/src/bin/genewise";
my $program = "/hwfssz1/ST_EARTH/P18Z10200N0107/P17Z10200N0101_Metazoa_RNA_Editing/liqiye/PC_PA_UN/bin/Wise2/wise2.2.0/src/bin/genewise";
foreach my $p (\$pepFile, \$nucleotideFile, \$outDir, \$strandList) {
	$$p = abs_path($$p);
}
#################################################################################################################


##############################################################
## prepare sequences
my $pepOutDir = &get_sequence_outDir($outDir, "pep");
&sequence_parser($pepFile, $pepOutDir, "pep");
my $nucleotideOutDir = &get_sequence_outDir($outDir, "nuc");
&sequence_parser($nucleotideFile, $nucleotideOutDir, "nuc");
##############################################################


##############################################################
my %geneStrand;
if ($strandList) {
	&strandList_parser($strandList, \%geneStrand);
}
##############################################################


#######################################################################
## call genewise
my %pepToFile;
opendir PEP, $pepOutDir;
while (my $file = readdir PEP) {
	next unless $file =~ /(.+)\.fa$/;
	my $gene = $1;
	my $in_file = "$pepOutDir/$file";
	$pepToFile{$gene} = $in_file;
}
closedir PEP;

my $genewiseOutDir = &get_sequence_outDir($outDir, "result");
my @shell;
opendir NUC, $nucleotideOutDir;
while (my $file = readdir NUC) {
	next unless $file =~ /(.+)\.fa$/;
	my $gene = $1;
	next unless $pepToFile{$gene};
	my $in_file = "$nucleotideOutDir/$file";
	my $out_file = $genewiseOutDir . "/" . basename($in_file) . ".gw";
	my $compare_direct;
	if ($geneStrand{$gene}) {
		$compare_direct = ($geneStrand{$gene} eq "+") ? "tfor" : "trev";
	} else {
		$compare_direct = "both";
	}
	push @shell, "$program $pepToFile{$gene} $in_file -genesf -$compare_direct -quiet >$out_file\n";
}
closedir NUC;

## split shell
my @sh_files;
my $total_line_num = @shell;
my $each_file_num = int($total_line_num / $job_num + 1);
my ($file_count, $line_count) = (1, 0);
my $sh_file = "$genewiseOutDir/${keyword}_gw$file_count.sh";
open SH, ">$sh_file";
print SH "export WISECONFIGDIR=/hwfssz1/ST_EARTH/P18Z10200N0107/P17Z10200N0101_Metazoa_RNA_Editing/liqiye/PC_PA_UN/bin/Wise2/wise2.2.0/wisecfg\n";
push @sh_files, $sh_file;
foreach (@shell) {
        $line_count ++;
	print SH $_;
	unless ($line_count % $each_file_num) {
		close SH;
		$file_count ++;
		$sh_file = "$genewiseOutDir/${keyword}_gw$file_count.sh";
		open SH, ">$sh_file";
		print SH "export WISECONFIGDIR=/hwfssz1/ST_EARTH/P18Z10200N0107/P17Z10200N0101_Metazoa_RNA_Editing/liqiye/PC_PA_UN/bin/Wise2/wise2.2.0/wisecfg\n";
		push @sh_files, $sh_file;
	}
}
close SH;
## qsub shell file
foreach (@sh_files) {
	#print "$_\n";
	system "qsub -S /bin/sh -cwd -l vf=1G,num_proc=1 -binding linear:1 -q st.q -P P18Z10200N0107 $_";
	#system "sh $_ &";
}
#######################################################################################

## subroutine
#################################################################################################
sub sequence_parser {
my $in_file = shift;
my $out_dir = shift;
my $type = shift;
open IN, $in_file;
open OUT, $in_file;
my $out_file;
while (<IN>) {
	if (/^\s*>/) {
		close OUT;
		my $id;
		if ($type eq "nuc") {
			$id = (split /\s+/)[0];
			$id =~ s/^>//;
			$id = (split /##/, $id)[0];
		} elsif ($type eq "pep") {
			$id = (split /\s+/)[0];
			$id =~ s/^>//;
		}
		$out_file = "$out_dir/$id.fa";
		$out_file =~ s/\|/_/g;
		open OUT, ">$out_file";
		print OUT;
	} else {
		s/^\s+$//;
		s/^\s+|\s+$//g;
		print OUT uc($_), "\n" if $_;
	}
}
close IN;
}
#####################

#####################
sub get_sequence_outDir {
my $out_dir = shift;
my $type = shift;
my $seq_out_dir = "$out_dir/$type";
mkdir $seq_out_dir unless -e $seq_out_dir;
return $seq_out_dir;
}
#####################


#####################
sub strandList_parser {
my $in_file = shift;
my $ref = shift;
open IN, $in_file;
while (<IN>) {
	my ($gene, $strand) = (split /\s+/)[0,1];
	$ref->{$gene} = $strand;
}
close IN;
}
#####################
