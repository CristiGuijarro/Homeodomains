#!/usr/bin/perl -w
use strict;

my $dir = $ARGV[0];
my @files = `ls $dir/*.trimmed`;
foreach my $file (@files) {
	chomp $file;
	my %hash = ();
	my @newfasta = ();
	my $fastas = `grep [A-Z\-] $file`;
	#print "grep [A-Z] $file\n";
	my $hg;
	if ($file =~ /([0-9]+)\.trimmed/) {
	#	print $1, "\n";
	#	exit;
		$hg = $1;
	}
	my $newfile = "$dir/$hg.sorted.fasta";
	#print "$newfile\n";
	my @fastas = split />/, $fastas;
	foreach my $fasta (@fastas) {
		chomp $fasta;
		#print $fasta, "\n";
		#$fasta =~ s/\n//g;
		my ($header, $sequence) = split / bp/, $fasta;
		$hash{$header} = $sequence;
		
	}
	foreach my $fasta (@fastas) {
		chomp $fasta;
		my ($header, $sequence) = split / bp/, $fasta;
		$hash{$sequence}++;
		if ($header =~ /Unknown/) {
			unless ($hash{$sequence} > 1) {
				#my @matches = grep { $hash{$_} eq $sequence } keys %hash;
				my @matches = grep { /$sequence/ } keys %hash;
				my @number = grep /Unknown/, @matches;
				if (scalar(@number) == scalar(@matches)) {
					unless (length($sequence) > 85) {
						unless (length($sequence) < 35) {
							chomp $header;
							chomp $sequence;
							$sequence =~ s/^\n//;
							push @newfasta, ">$header bp\n$sequence";
						}
					}
				}
			}
		}
		else {
			unless (length($sequence) > 85) {
				unless (length($sequence) < 60) {
					chomp $header;
					chomp $sequence;
					$sequence =~ s/^\n//;
					push @newfasta, ">$header bp\n$sequence";
				}
			}
		}
	}
	open (OUT, ">$newfile");
	print OUT join("\n", @newfasta);
}