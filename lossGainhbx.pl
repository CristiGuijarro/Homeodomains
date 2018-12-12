#!/usr/bin/perl -w
use strict;

my $recountfile = $ARGV[0];
my $originsfile = $ARGV[1];
chomp $recountfile;
chomp $originsfile;
my @recounts = `grep -v "0\$" $recountfile`;
my @origins = `grep '[A-Z]' $originsfile`;
my %ancestral = ();

foreach my $originline (@origins){
	chomp $originline;
	my ($class, $family, $group) = split /,/, $originline;
	$ancestral{$group}++;
}
while (my ($keys, $values) = each %ancestral) {
	print "$keys\t$values\n";
}
exit;