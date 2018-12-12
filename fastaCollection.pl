#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);
use Getopt::Long;

#USEAGE: ./fastaCollection.pl --terms "GPCR,G-protein coupled receptor,rhodopsin" --occfile "GPCR/GPCR_Occupancy.csv" --fastafile "GPCR/GPCR.fasta"

my $listofterms = "";
my $occfile = "";
my $fastafile = "";
GetOptions ('terms=s' => \$listofterms, 'occfile=s' => \$occfile, 'fastafile=s' => \$fastafile) or die "--listofterms in a comma seperated string list in \"\" e.g. \"GPCR,G-protein,...\", --occfile and --fastafile names also need to be supplied\n\n";
my @listofterms = split /,/, $listofterms;
my @homologygroups = ();
my @fastas = ();
my @occupancy = ();
my %hasta = ();
my %genes = ();
my $occupancytable = "/home/cristig/Desktop/ComparativeGenomicsP1/fullOccupancy.csv";
foreach my $term(@listofterms) {
	chomp $term;
	my $term2 = $term;
	$term2 =~ s/ /_/g;
	#my @homologys = `echo 'SELECT DISTINCT Homology FROM geneontology WHERE Family LIKE \"%$term%\" AND Homology IN (SELECT DISTINCT Homology FROM genescollected WHERE Header LIKE \"%$term2%\") ORDER BY Homology ASC;' | mysql -B -uweb ComparativeGenomics | tail -n+2`;
	my @homologys = `echo 'SELECT DISTINCT Homology FROM geneontology WHERE Family LIKE \"%$term%\" ORDER BY Homology ASC;' | mysql -B -uweb ComparativeGenomics | tail -n+2`;
	foreach my $homology (@homologys) {
		chomp $homology;
		#my $count = `echo "SELECT COUNT(Header) AS Total FROM genescollected WHERE Homology=$homology AND Header LIKE '%$term2%' GROUP BY Homology;" | mysql -B -uweb ComparativeGenomics | tail -n+2`;
		#if ($count > 2) {
		my $genes = `echo 'SELECT DISTINCT Genes FROM geneontology WHERE NOT Genes REGEXP "_" AND Family REGEXP "homeobox family" AND Genes NOT LIKE "v%" AND NOT Genes LIKE "Loc%" AND Homology=$homology LIMIT 1;'| mysql -B -uweb ComparativeGenomics | tail -n+2`; 
		$genes =~ s/^ //;
		if ($genes =~ / /) {
			my @genes = split / /, $genes;
			$genes = $genes[0];
		}
		chomp $genes;
		$genes = ucfirst($genes);
		$genes{$homology} = $genes;
		if ($genes =~ /[A-Za-z]/) {
			push @homologygroups, "$homology";
		}
		#}
	}
}
push @occupancy, `head -n 1 $occupancytable`;
@homologygroups = uniq(@homologygroups);
open(LOG, ">fastaCollection.log");
foreach my $homology (@homologygroups) {
	chomp $homology;
	print "$homology ",$genes{$homology}, "\n";
	#my ($homology,$gene) = split /,/, $homologys;
	my $occupancy = `grep "^$homology," $occupancytable`;
	chop $occupancy;
	push @occupancy, $occupancy;
	#my $gene = `echo "SELECT DISTINCT Genes FROM geneontology WHERE Homology=$homology AND NOT Genes LIKE "%_%" AND NOT Genes LIKE "% %" LIMIT 1;" | mysql -B -uweb ComparativeGenomics | tail -n+2`;
	my @fastas = `echo "SELECT DISTINCT s.Header, s.Sequence FROM genescollected g, sequences s WHERE s.Header=g.Header AND g.Homology=$homology;" | mysql -B -uweb ComparativeGenomics | tail -n+2`;
	my $count = 0;
	foreach my $fasta (@fastas) {
		chomp $fasta;
		$count++;
		my ($header,$sequence) = split /\t/, $fasta;
		my $spec = "";
		if ($header =~ /(^[A-Za-z0-9]{4})\_/) {
			$spec = $1;
		}
		chomp $header; chomp $sequence;
		my $gene = $genes{$homology};
		if ($header =~ /^[A-Za-z0-9]{4}\_(.*)_ENS/) {
			$gene = $1;
		}
		print LOG "$header -> >$homology"."_$gene"."_$spec"."_$count\n";
		$hasta{">$homology"."_$gene"."_$spec"."_$count"} = $sequence;
	}
}
close(LOG);
open(FASTA, ">$fastafile");
print FASTA join("\n", @occupancy), "\n";
close(FASTA);
open(OCC, ">$occfile");
foreach my $keys(sort keys %hasta) {
	print OCC "$keys\n",$hasta{$keys}, "\n";
}
close(OCC);
exit;