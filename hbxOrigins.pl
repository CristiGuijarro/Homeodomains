#!/usr/bin/perl -w
use strict;
use lib '/home/cristig/Desktop/ComparativeGenomics';
use queryMySQL;
use List::UtilsBy qw(max_by);
use List::Util qw(max);
my $hbxcountfile = $ARGV[0];
my @hbxcount = `grep -v -P ',0' $hbxcountfile`;
my %hash = ();
my %homeos = (
	Human => "Homo sapiens",
	Beetle => "Tribolium castaneum",
	Amphioxus => "Branchiostoma floridae",
	Frog => "Xenopus tropicalis",
	Chicken => "Gallus gallus",
	Fruitfly => "Drosophila melanogaster",
	Zebrafish => "Danio rerio",
	Nematode => "Caenorhabditis elegans",
);
foreach my $hbx (@hbxcount) {
	chomp $hbx;
	my ($spec, $group, $class, $family, $count) = split /,/, $hbx;
	if ($homeos{$spec}) {
		$spec = $homeos{$spec};
	}
	$hash{"$class,$family"}.="$spec,";
}
foreach my $keys (sort keys %hash) {
	my $highindex = 0;
	$hash{$keys} =~ s/,$//;
	my @specs = split /,/, $hash{$keys};
	my $specount = scalar(@specs);
	my %logeny = ();
	my @phylos = ("Kingdom", "Ten", "Nine", "Eight", "Seven", "Six", "Five", "Four", "Three");
	foreach my $spec (@specs) {
		chomp $spec;
		foreach (my $i = 0; $i < @phylos; $i++) {
			#print $phylos[$i], "\n";
			my $type = "";
			while ($type eq "") {
				$type = `echo 'SELECT DISTINCT $phylos[$i] FROM phylogeny WHERE Species="$spec"' | mysql -B -uweb ComparativeGenomics`;
			#print "echo 'SELECT DISTINCT $phylos[$i] FROM phylogeny WHERE Species=\"$spec\"' | mysql -B -uweb ComparativeGenomics\n";
			}
			$logeny{$type}++;
			#print $type, "\n";
		}
	}
	my $highest = max values %logeny;
	foreach my $type (sort keys %logeny) {
		if ($logeny{$type} == $highest) {
			chomp $type;
			#print $logeny{$type}, "\n";
			my ($column, $cladename) = split /\n/, $type;
			#print $type, "\n";
			my $index = grep { "$phylos[$_]" eq "$column" } 0..$#phylos;
			if ($highindex < $index) {
				$highindex = $index;
			}
			if ($index => $highindex) {
				$hash{$keys} = $cladename;
			}
		}
	}
}
while (my ($keys, $values) = each %hash) {
	print $keys,",",$values, "\n";
}
exit;