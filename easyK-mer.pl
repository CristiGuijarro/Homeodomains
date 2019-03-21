#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
use Statistics::Basic::Median;
my $dir = $ARGV[0];
chomp $dir;
my $whole = `grep -h [A-Z] $dir/*.sorted.fasta`;
$whole =~ s/^>//; 
my @fastas = split /\n>/, $whole;
my %hash = ();
my %otherhash = ();
my %hashy = ();
my $cons = "RRKRTTFTKEQLQELEKEFQQNKYPSXEEREELAQKLGLTETQVKVWFQNRRAKWKKQQ";
$hash {"Consensus\tConsensus\tConsensus\tConsensus\tConsensus\tConsensus"} = aa_collector($cons);
foreach my $fasta (@fastas) {
	chomp $fasta;
	my ($header, $sequence) = split /\n/, $fasta;
	my ($spec,$group,$family,$class,$HG,$uid) = split /\|/, $header;
	my $phylum = phylum_select("$spec");
	if ($phylum) {
		$hash{"$spec\t$phylum\t$family\t$class\t$HG\t$uid"} = aa_collector($sequence);
	}
}
my (@aacheck) = aa_combinations();
while (my ($keys, $values) = each %hash) {
	chomp $keys;
	chomp $values;
	my ($spec,$phylum,$family,$class,$HG,$uid) = split /\t/, $keys;
	unless ($class eq "Consensus") {
		$hashy{$class}++;
	}
	my @values = split /\t/, $values;
	foreach (my $i = 0; $i < @values; $i++) {
		if ($phylum ne "Consensus") {
			$otherhash{"$phylum\t$family\t$class\t$aacheck[$i]"}.="$values[$i]\t";
		}
		else {
			while (my ($keys, $newvalues) = each %hashy) {
				$otherhash{"Consensus\tConsensus\t$keys\t$aacheck[$i]"}="$values[$i]";
			}
		}
	}
}

while (my ($keys, $values) = each %otherhash) {
	chomp $keys;
	chomp $values;
	my @values = split /\t/, $values;
	my $mean = mean(@values);
	my $median = median(@values);
	if ($keys =~ /Consensus/) {
		$median = $values;
	}
	print "$keys\t$median\n";
}
############################################
sub aa_combinations {
	my @aas = qw {
		A
		R
		N
		D
		C
		E
		Q
		G
		H
		I
		L
		K
		M
		F
		P
		S
		T
		W
		Y
		V
	};
	my %subhash = ();
	my @array = ();
	for (my $i = 0; $i<@aas; $i++) {
		for (my $j = 0; $j<@aas; $j++) {
			chomp $aas[$i];
			chomp $aas[$j];
			$subhash{$aas[$i].$aas[$j]}++;
		}
	}
	@array = ( sort keys %subhash );
	return(@array);
}
sub aa_consensus_window {
	my $consensus = "RRKRTTFTKEQLQELEKEFQQNKYPSXEEREELAQKLGLTETQVKVWFQNRRAKWKKQQ"; #From EMBOSS cons of all Homeodb sequences.
	my @array = ();
	my %subhash = ();
	my @amacs = split //, $consensus;
	foreach (my $i = 0; $i+1 < @amacs; $i++) {
		#print "\"$amacs[$i]$amacs[$i+1]\", ";
		$subhash{$amacs[$i].$amacs[$i+1]}++;
	}
	#print "\n";
	@array = ( sort keys %subhash );
	return(@array);
}
sub aa_collector {
	my $sequence = shift;
	chomp $sequence;
	#my (@aa) = aa_combinations();
	my (@aa) = aa_consensus_window();
	my @counter = ();
	foreach my $aa (@aa) {
		chomp $aa;
		my @c = $sequence =~ /$aa/g;
		my $count = @c;
		push @counter, $count;
	}
	my $string = join("\t", @counter);
	return $string;
}
sub phylum_select {
	my $spec = shift;
	chomp $spec;
	my %keys = (
		drer => "Chordata",
		xtro => "Chordata",
		acar => "Chordata",
		ggal => "Chordata",
		hsap => "Chordata",
		odio => "Chordata",
		cint => "Chordata",
		bsch => "Chordata",
		bflo => "Chordata",
		bbel => "Chordata",
		pmin => "Echinodermata",
		spur => "Echinodermata",
		lvar => "Echinodermata",
		skow => "Hemichordata",
		pfla => "Hemichordata",
		hduj => "Tardigrada",
		rvar => "Tardigrada",
		smim => "Arthropoda",
		agen => "Arthropoda",
		isca => "Arthropoda",
		mmar => "Arthropoda",
		lpol => "Arthropoda",
		dpul => "Arthropoda",
		phaw => "Arthropoda",
		tcas => "Arthropoda",
		dmel => "Arthropoda",
		znev => "Arthropoda",
		smar => "Arthropoda",
		bmal => "Nematoda",
		cele => "Nematoda",
		tspi => "Nematoda",
		rcul => "Nematoda",
		hrob => "Annelida",
		ctel => "Annelida",
		lana => "Brachiopoda",
		cgig => "Mollusca",
		pfuc => "Mollusca",
		obim => "Mollusca",
		lgig => "Mollusca",
		ilin => "Orthonectida",
		hmic => "Platyhelminthes",
		emul => "Platyhelminthes",
		gsal => "Platyhelminthes",
		smed => "Platyhelminthes",
		mlig => "Platyhelminthes",
		sjap => "Platyhelminthes",
		avag => "Rotifera",
		apal => "Cnidaria",
		nvec => "Cnidaria",
		adig => "Cnidaria",
		hmag => "Cnidaria",
		tkit => "Cnidaria",
		pbac => "Ctenophora",
		mlei => "Ctenophora",
		tadh => "Placozoa",
		lcom => "Porifera",
		scil => "Porifera",
		aque => "Porifera",
		ocar => "Porifera",
	);
	my $assoc = $keys{$spec};
	return $assoc;
}
sub mean {
    return sum(@_)/@_;
}
sub median {
	my @array = @_;
	my $count = scalar(@array);
	@array = sort { $a <=> $b } @array;
	my $median = 0;
	if ($count % 2) { 
		$median = $array[($count/2)]; 
	}
	else { 
		$median = ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	}
	return $median;
}
exit;