#!/usr/bin/perl -w
use strict;
#sed -e '/^atha|Viridiplantae/ s/atha|/atha|Viridiplantae|/ ; s/atha|Viridiplantae/atha|Viridiplantae/' FastasAnnot/*.fasta | grep 'atha'
my $dir = $ARGV[0];
chomp $dir;
my @headers = `grep -h '>' $dir/*`;
my %homeobox = ();
my %spec = ();
foreach my $header (@headers) {
	chomp $header;
	my $spec = "";
	my $kingdom = "";
	my $family = "";
	my $class = "";
	my $HG = "";
	my $uid = "";
	if ($header =~ /\|.*\|.*\|.*\|.*\|/) {
		($spec,$kingdom,$family,$class,$HG,$uid) = split /\|/, $header;
		$spec =~ s/^>//;
		if ($family eq $class) {
			my ($Class,$Family,$cm,$End) = split /,/, `grep -i '$family' HomeoboxClass.csv`;
			if ($Class) {
				#$homeobox{"$class|$family"} = "$Class|$Family";
			}
			else {
				$homeobox{"$class|$family"} = "$class|$family";
			}
		}
	}
	else {
		unless ($header !~ /\|/) {
			($spec,$family,$class,$HG,$uid) = split /\|/, $header;
			$spec =~ s/^>//;
			$spec = lc $spec;
			my $groups = `grep '$spec' phylogenyTable.csv`;
			my @groups = split /,/, $groups;
			$kingdom = $groups[6];
			$kingdom =~ s{\/}{-}g;
			$spec{$spec} = "$spec|$kingdom";
			if ($family eq $class) {
				my ($Class,$Family,$cm,$End) = split /,/, `grep -i '$family' HomeoboxClass.csv`;
				if ($Class) {
					$homeobox{"$class|$family"} = "$Class|$Family";
				}
			}
		}
	}
}
foreach my $orig ( sort keys %homeobox ) {
	my $new = $homeobox{$orig};
	chomp $orig;
	chomp $new;
	print "$orig -> $new\n";
}
foreach my $orig ( sort keys %spec ) {
	my $new = $spec{$orig};
	chomp $orig;
	chomp $new;
	print "$orig -> $new\n";
}
exit;