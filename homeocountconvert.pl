#!/usr/bin/perl -w
use strict;

my $csv = $ARGV[0];
my $hbxclass = `grep 'Homeobox,' $csv`;
my $hbxfam = `grep '^Species,' $csv`;
my $therest = `grep 'Metazoa,' $csv`;
my %hash = ();
my @hbxclass = split /,/, $hbxclass;
@hbxclass = @hbxclass[ 3 .. $#hbxclass ];
my @hbxfam = split /,/, $hbxfam;
@hbxfam = @hbxfam[ 3 .. $#hbxfam ];
my @therest = split /\n/, $therest;
foreach (my $i = 0; $i < @hbxclass; $i++) {
	chomp $hbxclass[$i];
	chomp $hbxfam[$i];
	if ($hbxfam[$i] =~ /^([A-Z][a-z]+).*/) {
		$hbxfam[$i] = $1;
	}
	$hbxclass[$i].=",$hbxfam[$i]";
}
foreach my $counts (@therest) {
	my @parts = split /,/, $counts;
	my $animal = join(',', @parts[0..2]);
	$animal =~ s/,Metazoa//g;
	my @counts = @parts[ 3 .. $#parts ];
	foreach (my $i = 0; $i < @hbxclass; $i++) {
		if (!exists $hash{"$animal,$hbxclass[$i]"}) {
			$hash{"$animal,$hbxclass[$i]"} = $counts[$i];
		}
		else {
			$hash{"$animal,$hbxclass[$i]"} += $counts[$i];
		}
	}
}
while (my ($keys, $values) = each %hash) {
	print "$keys,$values\n";
}