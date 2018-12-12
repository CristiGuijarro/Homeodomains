#!/usr/bin/perl -w
use strict;

my $dir = $ARGV[0];
my @files = `find ./$dir/ -type f -name "*.fasta"`;
my $newdir = "Reduced";
unless (-d $newdir) {
	mkdir $newdir;
}
foreach my $file (@files) {
	my %hash = ();
	chomp $file;
	my $sequences = `grep -A1 '>' $file`;
	(my $nada, my $olddir, $file) = split /\//, $file;
	my $newfile = "$newdir/$file";
	my @individuals = split />/, $sequences;
	foreach my $indiv (@individuals) {
		chomp $indiv;
		if ($indiv =~ /[A-Za-z0-9]/) {
			my ($header, $sequence) = split /\n/, $indiv;
			$hash{">$header"} = "$sequence";
		}
	}
	open (OUT, ">$newfile");
	foreach my $keys (sort keys %hash) {
		chomp $keys;
		chomp $hash{$keys};
		my $values = $hash{$keys};
		print OUT "$keys\n$values\n";
	}
	close(OUT);
}
exit;
