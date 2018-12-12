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
	my ($header, $seq) = split /\n/, $seq;
	$hash{$seq} = $header;
	#">$homology"."_$gene"."_$spec"."_$count" - header format
}
my %fasta = ();
while (my ($values, $keys) = each %hash) {
	chomp $keys;
	chomp $values;
	$values =~ s/[\*\_\.]//g;
	my ($homology,$gene,$spec,$count) = split /_/, $keys;
	print ">$keys\n$values\n";
	$fasta{$homology} .= ">$keys\n$values\n";
}
while (my ($keys, $values) = each %fasta) {
	open(FASTA, ">$keys.fasta");
	print FASTA $values;
	close FASTA;
}
exit;