#!/usr/bin/perl -w
use strict;
use Lingua::EN::Nums2Words;
use Lingua::EN::Words2Nums;
use Array::Utils qw(:all);
my $recountfile = $ARGV[0];
chomp $recountfile;
my @recounts = `grep -v "0\$" $recountfile`;
my %matrix = ();
my %ancestral = ();
my %trueancestral = ();
my %lost = ();
my %novel = ();
my %order = ();
my $cladename = $_[0];
my $clade = $_[1];
foreach my $count (@recounts) {
	chomp $count;
	my ($species, $groups, $geneclass, $genefamily, $count) = split /,/, $count;
	my $phylogeny = `grep '$species' phylogenyTable.csv`;
	$phylogeny =~ s/\t//g;
	my @array1 = split /,/, $phylogeny;
	my @array = splice(@array1,7,10);
	foreach (my $i = 0; $i< @array; $i++) {
			if (exists $order{$array[$i]}) {
				if ((10-$i) > words2nums($order{$array[$i]})) {
					my $value = num2word(10-$i);
					$value = lc($value);
					$value = ucfirst($value);
					$order{$array[$i]} = $value;
				}
			}
			else {
				my $value = num2word(10-$i);
				$value = lc($value);
				$value = ucfirst($value);
				$order{$array[$i]} = $value;
			}
		}
	my $array = join("\t", @array);
	$matrix{"$geneclass $genefamily"} .= "$array\t";
}
my %hash = ();
while (my ($hbx, $phylogenies) = each %matrix) {
	chomp $hbx; chomp $phylogenies;
	$phylogenies =~ s/\t$//;
	$phylogenies =~ s/^\t//;
	my @arrays = split /\t/, $phylogenies;
	foreach my $group (@arrays) {
		chomp $group;
		$group =~ s/^\t//;
		$hash{"$group\t$hbx"}++;
	}
}
while (my ($key, $count) = each %hash) {
	chomp $key; chomp $count;
	my ($cladename, $hbx) = split /\t/, $key;
	push (@{ $ancestral{$cladename} }, $hbx);
}
my %temp = %ancestral;
while (my ($cladename, $hbx) = each %temp) {
	if (exists $order{$cladename}) {
		my @sisterclades = ();
		my @lowerclades = ();
		my @sisterhbxes = ();
		my @ancestralhbxes = ();
		my @ancestraldups = ();
		my @sisterdups = ();
		my %multiples = ();
		my %ancestraldups = ();
		my $clade = $order{$cladename};
		my $sisterclade = sister_clade($cladename, $clade);
		my $lowerclades = lower_clade($cladename, $clade);
		if ($sisterclade !~ /NONE/) {
			@sisterclades = split /\t/, $sisterclade;
			unless (scalar(@sisterclades) < 1) {
				foreach my $sister (@sisterclades) {
					foreach my $i (0 .. $#{ $ancestral{$sister} } ) {
						push @sisterhbxes, $ancestral{$sister}[$i];
					}
				}
			}
		}
		@lowerclades = split /\t/, $lowerclades;
		if (scalar(@lowerclades) > 1) {
			#print "1. $cladename $lowerclades $hbx\n";
			foreach my $lower (@lowerclades) {
				foreach my $i (0 .. $#{ $ancestral{$lower} } ) {
					push @ancestralhbxes, $ancestral{$lower}[$i];
				}
			}
		}
		else {
			#print "2. $cladename $lowerclades $hbx\n";
			foreach my $i (0 .. $#{ $ancestral{$cladename} } ) {
				push @ancestraldups, $ancestral{$cladename}[$i];
			}
		}
		for (@sisterhbxes) {
			push @sisterdups, $_  if ($multiples{$_}++ >= 2)
		}
		for (@ancestralhbxes) {
			push @ancestraldups, $_ if ($ancestraldups{$_}++ == 1)
		}
		my @loss = array_minus(@sisterdups,@{ $ancestral{$cladename} });
		@loss = unique(@loss);
		push @{ $lost{$cladename} }, @loss;
		my @novel = array_minus(@{ $ancestral{$cladename} },@sisterhbxes);
		@novel = unique(@novel);
		push @{ $novel{$cladename} }, @novel;
		my @ancestral = unique(@ancestraldups);
		#my @ancestral = @ancestraldups;
		push @{ $trueancestral{$cladename} }, @ancestral;
	}
}
print "Clade\tAncestral\tNovel\tLost\n";
for my $k (sort keys %ancestral){
	unless (not defined $novel{$k}) {
		my $ancestral = @{ $trueancestral{$k} };
		my $novel = @{ $novel{$k} };
		my $lost = @{ $lost{$k} };
		print "$k\t$ancestral\t$novel\t$lost\n";
	}
}
###########################################################
sub higher_clade {
	my $cladename = $_[0];
	my $clade = $_[1];
	my @hierarchy = qw (
		Superdomain
		Domain
		Subdomain
		Group
		Kingdom
		Ten
		Nine
		Eight
		Seven
		Six
		Five
		Four
		Three
		Two
		One
		Class
		Order
	);
	my $cladeindice = 0;
	foreach (my $i = 0; $i < @hierarchy; $i++) {
		chomp $hierarchy[$i]; chomp $clade;
		if ($hierarchy[$i] eq $clade) {
			$cladeindice = $i;
		}
	}
	my $higherclade = $hierarchy[$cladeindice-1];
	my @cladelist = `echo 'SELECT DISTINCT \`$higherclade\` FROM phylogeny;' | mysql -B -uweb ComparativeGenomics`;
	foreach my $cladelist (@cladelist) {
		if ($cladelist eq $cladename) {
			$higherclade = $hierarchy[$cladeindice-2];
			next;
		}
	}
	my @higherclade = `echo 'SELECT DISTINCT \`$higherclade\` FROM phylogeny WHERE $clade="$cladename" AND Kingdom="Metazoa";' | mysql -B -uweb ComparativeGenomics`;
	return ($higherclade[1]);
}
sub lower_clade {
	my $cladename = $_[0];
	my $clade = $_[1];
	chomp $cladename;
	chomp $clade;
	my @hierarchy = qw (
		Superdomain
		Domain
		Subdomain
		Group
		Kingdom
		Ten
		Nine
		Eight
		Seven
		Six
		Five
		Four
		Three
		Two
		One
		Class
		Order
	);
	my $cladeindice = 0;
	foreach (my $i = 0; $i < @hierarchy; $i++) {
		chomp $hierarchy[$i]; chomp $clade;
		if ($hierarchy[$i] eq $clade) {
			$cladeindice = $i;
		}
	}
	
	my $lowerclade = $hierarchy[$cladeindice+1];

	my $cladelist = `echo 'SELECT DISTINCT \`$lowerclade\` FROM phylogeny WHERE $clade="$cladename" AND Kingdom="Metazoa";' | mysql -B -uweb ComparativeGenomics`;
	my @cladelist = split /\n/, $cladelist;
	shift @cladelist;
	foreach my $cladelist (@cladelist) {
		if ($cladelist eq $cladename) {
			$lowerclade = $hierarchy[$cladeindice+2];
			$cladelist = `echo 'SELECT DISTINCT \`$lowerclade\` FROM phylogeny WHERE $clade="$cladename" AND Kingdom="Metazoa";' | mysql -B -uweb ComparativeGenomics`;
			@cladelist = split /\n/, $cladelist;
			shift @cladelist;
			next;
		}
	}
	my $lowerclades = join("\t", @cladelist);
	return $lowerclades;
}
sub sister_clade {
	my $cladename = $_[0];
	my $clade = $_[1];
	chomp $cladename;
	chomp $clade;
	my $cladelist = `echo 'SELECT DISTINCT \`$clade\` FROM phylogeny WHERE not $clade="$cladename" AND Kingdom="Metazoa";' | mysql -B -uweb ComparativeGenomics`;
	my @cladelist = split /\n/, $cladelist;
	shift @cladelist;
	my $sisterclades = join("\t", @cladelist);
	if (scalar(@cladelist) < 1) {
		$sisterclades = "NONE";
	}
	return $sisterclades;
}
exit;