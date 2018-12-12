#!/usr/bin/perl -w
use strict;
my $dir = $ARGV[0];
chomp $dir;
my @fastas = `ls $dir/*.fasta`;
my $newdir = "IntermediateFiles";
unless (-d $newdir) {
	mkdir $newdir;
}
unless (-d "RAXML") {
	mkdir "RAXML";
}
unless (-d "Fasttree") {
	mkdir "Fasttree";
}
foreach my $hg (@fastas) {
	chomp $hg;
	$hg =~ s/\.fasta$//g;
	my $HG = $hg;
	if ($HG =~ /$dir\/(.*)/) {
		$HG = $1;
	}
	#if ($HG == 62) {
	#system "mafft --auto --leavegappyregion $dir/$HG.fasta > $newdir/$HG.aln";
	#system "muscle -in $dir/$HG.fasta -out $newdir/$HG.aln";
	#system "muscle -in $newdir/$HG.aln -refine -out $newdir/$HG.refined.aln";
	#system "trimal -in $newdir/$HG.refined.aln -fasta -out $newdir/$HG.trimmed -gappyout";
	#print "trimal -in $newdir/$HG.refined.aln -fasta -out $newdir/$HG.trimmed -gappyout\n\n";
	unless (-f "$HG.trimmed.contree") {
		system "iqtree -s IntermediateFiles/$HG.trimmed -nt AUTO -m MFP+MERGE -bb 1000";
	}
	#}
}
exit;