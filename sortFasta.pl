#!/usr/bin/perl -w
use strict;

my $fasta = $ARGV[0];
chomp $fasta;
my $seqs = `grep "" $fasta`;
my @seqs = split />/, $seqs;
shift @seqs;
my %hash = ();
foreach my $seq (@seqs) {
	chomp $seq;
	my ($wholeheader, $seq) = split /\n/, $seq;
	$hash{$seq} = $wholeheader;
	#">$homology"."_$gene"."_$spec"."_$count" - header format
	#>original|$homology|$count
}
my %fasta = ();
while (my ($values, $keys) = each %hash) {
	chomp $keys;
	chomp $values;
	$values =~ s/[\*\_\.]//g;
	#my ($homology,$gene,$spec,$count) = split /_/, $keys;
	my ($original,$homology,$count) = split /\|/, $keys;
	$original =~ s/[\|\-\.]/_/g;
	my $spec;
	if ($original =~ /(^[A-Za-z0-9]{4})_/) {
		$spec = $1;
	}
	my $end = substr $original, 5;
	$original = "$spec|$end";
	$keys = "$original|$homology|$count";
	print ">$keys\n$values\n";
	unless ($keys =~ /mlig_maker/) {
		$fasta{$homology} .= ">$keys\n$values\n";
	}
}
while (my ($keys, $values) = each %fasta) {
	open(FASTA, ">$keys.fasta");
	print FASTA $values;
	close FASTA;
}
exit;