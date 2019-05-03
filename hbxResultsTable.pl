#!/usr/bin/perl -w
use strict;

my $fams = `grep ',' HomeoDBfamilies.csv`;
my %fams = split /[\n,]/, $fams;
my @fams = ();
while (my ($keys, $values) = each %fams) {
	#push @fams, "$keys\_$values";
	push @fams, "$values\t$keys";
#	print "0.$values\t$keys\n";
}
my $species = `grep '[A-Z]' speciesList.tsv`;
my @species = split /\n/, $species;
foreach my $specie (@species) {
	chomp $specie;
	($specie, my $code) = split /\t/, $specie;
	#print "0.$specie\n";
}
my $classlog = $ARGV[0];
chomp $classlog;
my $entries = `grep '_' $classlog`;
my %hash = ();
my @entries = split /\n/, $entries;
foreach my $entry (@entries) {
	chomp $entry;
	my ($spec,$header,$uid,$fasta,$homeodb,$interpro,$classif,$phylo) = split /\t/, $entry;
	my $class = "";
	my $family = "";
	if ($phylo =~ /_/) {
		$classif =~ s/_/\./g;
		my $homeo = "";
		if (($homeo) = grep/$classif/, @fams) {
			chomp $homeo;
			($family,$class) = split /\t/, $homeo;
			#print "1.$class\t$family\n";
		}
	}
	elsif ($classif !~ /Unknown/) {
		($family,$class) = split /\|/, $classif;
		#print "2.$class\t$family\n";
	}
	#else {
		#print $entry, "\n";
	#}
	my $gene = $header;
	if ($header =~ /GN=(.*)PE=/) {
		$gene = $1;
	}
	$hash{"$class\t$family"}{$spec}.="$gene ";
	#print "1.$spec\n";
}
#foreach my $name (sort keys %hash) {
#    foreach my $subject (keys %{ $hash{$name} }) {
#        print "$name, $subject: $hash{$name}{$subject}\n";
#    }
#}
#exit;
print "HBX Class\tHBX Family\t";
print join("\t", @species), "\n";
foreach my $fam (@fams) {
	chomp $fam;
	print "$fam\t";
	foreach my $taxon (@species) {
		chomp $taxon;
		if ($hash{$fam}{$taxon}) {
			print "$hash{$fam}{$taxon}\t";
		}
		else {
			print "-\t";
		}
	}
	print "\n";
}