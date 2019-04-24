#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw/ uniq /;

my $domains = `cat Homeo*Class*.csv | grep -v 'Family' | sed 's/,End//g'`;
my @domains = split /[\n,]/, $domains;
push @domains, "Homeobox";
push @domains, "Homeodomain";
my @domain = uniq @domains;
my @interpros;
foreach my $domain(@domain) {
	chomp $domain;
	my @interpro = `grep -h '$domain' *interpro.tsv`;
	push @interpros, @interpro;
}
@interpros = uniq @interpros;
open(INTER, ">new-interpro.tsv");
print INTER join("\n", @interpros);
close(INTER);
my %domainhash = ();
my $dir = "FastaDomains";
unless (-d $dir) {
	mkdir $dir;
}
foreach my $result (@interpros) {
	chomp $result;
	my ($query,$rubbish,$length,$analysis,$accession,$description,$start,$stop,$evalue,$status,$date,$tipid,$annotation) = split /\t/, $result;
	chomp $query;
	unless (exists $domainhash{$query}) {
		$domainhash{$query} = "$start\t$stop";
	}
	if (exists $domainhash{$query}) {
		my ($prevstart,$prevstop) = split /\t/, $domainhash{$query};
		my $prevlength = $prevstop-$prevstart;
		my $currentlength = $stop-$start;
		my $prevdiff = abs(60-$prevlength);
		my $currdiff = abs(60-$currentlength);
		if ($currdiff < $prevdiff) {
			$domainhash{$query} = "$start\t$stop";
		}
	}
}
my @selectedresults = ();
while (my($query,$values) = each %domainhash) {
	chomp $query;
	chomp $values;
	push @selectedresults, "$query\t$values";
}
foreach my $result (@selectedresults) {
	chomp $result;
	my ($query,$start,$stop) = split /\t/, $result;
	my $hg;
	chomp $query;
	if ($query =~ /\|([0-9]+)\|[0-9]+$/) {
		$hg = $1;
	}
	#print $query, "\n";
	#my $query2 = quotemeta($query);
	#print $query2, "\n";
	my $fasta = `grep -A1 '$query' Fastas/$hg.fasta`;
	#print $fasta, "\n";
	#print "grep -A1 '$query' Fastas/$hg.fasta\n";
	my ($header,$sequence) = split /\n/, $fasta;
	chomp $header; chomp $sequence; chomp $start; chomp $stop;
	print "$start\t$stop\t";
	my $end = $stop-$start;
	my $newseq = "$sequence";
	if (length($sequence) < 82) {
		$newseq = $sequence;
	}
	else {
	#	print length($sequence),"\n";
		if (($start > 10) && ($end+20 < length($sequence))) {
			if (length($sequence) < $start-10) {
				$newseq = substr $sequence, $start-10, $end+20;
			}
			else {
				$newseq = substr $sequence, $start, $end;
			}
		}
		else {
			$newseq = substr $sequence, 0, $end+10;
		}
	}
	unless ($newseq =~ /[A-Za-z]/) {
		$newseq = substr $sequence, $start, 60;
	}
	print "$header\t".length($sequence)."\t";
	print length($newseq), "\n";
	chomp $newseq;
	unless (-f "$dir/$hg.fasta") {
		system "touch $dir/$hg.fasta";
	}
	unless (!defined $newseq) {
		system "echo '$header\n$newseq' >> $dir/$hg.fasta";
	}
}
exit;