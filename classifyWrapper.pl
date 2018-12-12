#!/usr/bin/perl -w
use strict;

my $dir = $ARGV[0];
my @files = `ls $dir/*.fasta`;
foreach my $hg (@files) {
	chomp $hg;
	$hg =~ s/\.fasta//g;
	my @hg = split /\//, $hg;
	$hg = $hg[-1];
	print $hg, "\n";
	system "./homeoClassifier.pl --fastafile $dir/$hg\.fasta --blastfile All.blast --interprofile new-interpro.tsv";
	print "\n\n./homeoClassifier.pl --fastafile $dir/$hg\.fasta --blastfile All.blast --interprofile new-interpro.tsv\n\n***Complete***\n\n";
}
exit;