#!/usr/bin/perl -w
use strict;

my @phylos = qw {
	drer
	xtro
	acar
	ggal
	hsap
	odio
	cint
	bsch
	bflo
	bbel
	pmin
	spur
	lvar
	skow
	pfla
	hduj
	rvar
	smim
	agen
	isca
	mmar
	lpol
	dpul
	phaw
	tcas
	dmel
	znev
	smar
	bmal
	cele
	tspi
	rcul
	hrob
	ctel
	lana
	cgig
	pfuc
	obim
	lgig
	ilin
	hmic
	emul
	gsal
	smed
	mlig
	sjap
	avag
	apal
	nvec
	adig
	hmag
	tkit
	pbac
	mlei
	tadh
	lcom
	scil
	aque
	ocar
};
my $dir = $ARGV[0];
my $newdir = $ARGV[1];
chomp $newdir;
unless (-d $newdir) {
	mkdir $newdir;
}
chomp $dir;
my $whole = `grep -h [A-Z] $dir/*.sorted.fasta`;
$whole =~ s/^>//; 
my @fastas = split /\n>/, $whole;
my $cons = "RRKRTTFTKEQLQELEKEFQQNKYPSXEEREELAQKLGLTETQVKVWFQNRRAKWKKQQ";
my %hash = ();
foreach my $fasta (@fastas) {
	chomp $fasta;
	my ($header,$sequence) = split /\n/, $fasta;
	my ($spec,$phylum,$family,$class,$hg,$ui) = split /\|/, $header;
	$hash{"$phylum"."_$class"} .= ">$header\n$sequence\n";
}
foreach my $keys (sort keys %hash) {
	$hash{$keys} = ">HomeoDB|Consensus|Sequence|Emboss|001\n$cons\n".$hash{$keys};
}
foreach my $keys (sort keys %hash) {
	chomp $keys;
	my $newfile = "$newdir/$keys.fasta";
	open(OUT, ">$newfile");
	print OUT $hash{$keys};
	close(OUT);
	#system "mafft --auto --leavegappyregion $newfile > $newdir/$keys.aln";
	#system "trimal -in $newdir/$keys.aln -fasta -out $newdir/$keys.trimmed -gappyout";
	#system "echo 'trimal -in $newdir/$keys.aln -fasta -out $newdir/$keys.trimmed -gappyout\n' >> trimalRejig.sh";
	#system "iqtree -s $newdir/$keys.trimmed -nt AUTO -m MFP+MERGE -bb 1000";
}
my @alignmentfiles = `ls $newdir/*.trimmed`;
my %hashcounter = ();
my %indicehash = ();
foreach my $aln (@alignmentfiles) {
	chomp $aln;
	my $alignment =  `grep -h [A-Z] $aln`;
	$alignment =~ s/^>//;
	my @alignments = split /\n>/, $alignment;
	my $conseq = shift @alignments;
	foreach my $fasta (@alignments) {
		my ($header,$sequence) = split /bp\n/, $fasta;
		my ($spec,$phylum,$family,$class,$hg,$ui) = split /\|/, $header;
		my ($conheader,$consequence) = split /bp\n/, $conseq;
		$phylum = phylum_select($spec);
		chomp $phylum;
		if ($class =~ /ANTP/) {
			$class = class_correction($family);
			chomp $class;
		}
		if ( grep( /$spec/, @phylos ) ) {
			$consequence =~ s/ //g;
			$sequence =~ s/ //g;
			my $char = "RR";
			my $indice = index($consequence, $char);
			if ($indice > 5) {
				$char = "KR";
				$indice = index($consequence, $char);
			}
			my @chars = split //, $sequence;
			if (scalar(@chars>58)) {
				foreach (my $i = $indice; $i < @chars; $i++) {
					my $position = $i - $indice;
					if ($position < 60) {
						unless ($chars[$i] !~ /[A-Z\-]/) {
							unless ($class !~ /[A-Z]/) {
								$indicehash{"$phylum\t$class\t$position"} .= $chars[$i];
								$hashcounter{"$phylum\t$class\t$position"}++;
							}
						}
					}
				}
			}
		}
	}
}
my %actualfinalhash = ();
while (my ($keys, $values) = each %indicehash) {
	my @array = split //, $values;
	my %finalhash = ();
	foreach my $arr (@array) {
		$finalhash{$arr}++;
	}
	while (my ($aa, $count) = each %finalhash) {
		my $total = $hashcounter{$keys};
		my $multiplier = 100/$total;
		my $perc = $multiplier*$count;
		$actualfinalhash{"$keys\t$aa"} = $perc;
	}
}
while (my ($keys, $values) = each %actualfinalhash) {
	print "$keys\t$values\n";
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
sub class_correction {
	my $family = shift;
	my %keys = (
		"Cdx", "ANTP-HOXL",
		"Evx", "ANTP-HOXL",
		"Gbx", "ANTP-HOXL",
		"Gsx", "ANTP-HOXL",
		"Hox1", "ANTP-HOXL",
		"Hox2", "ANTP-HOXL",
		"Hox3", "ANTP-HOXL",
		"Hox4", "ANTP-HOXL",
		"Hox6-8", "ANTP-HOXL",
		"Hox9-13", "ANTP-HOXL",
		"Meox", "ANTP-HOXL",
		"Mnx", "ANTP-HOXL",
		"Pdx", "ANTP-HOXL",
		"Abox", "ANTP-NKL", 
		"Ankx", "ANTP-NKL", 
		"Barhl", "ANTP-NKL", 
		"Bari", "ANTP-NKL", 
		"Barx", "ANTP-NKL", 
		"Bsx", "ANTP-NKL", 
		"Dbx", "ANTP-NKL", 
		"Dlx", "ANTP-NKL", 
		"Emx", "ANTP-NKL", 
		"En", "ANTP-NKL", 
		"Hhex", "ANTP-NKL", 
		"Hlx", "ANTP-NKL", 
		"Hx", "ANTP-NKL", 
		"Lbx", "ANTP-NKL", 
		"Lex", "ANTP-NKL", 
		"Msx", "ANTP-NKL", 
		"Msxlx", "ANTP-NKL", 
		"Nanog", "ANTP-NKL", 
		"Nedx", "ANTP-NKL", 
		"Nk1", "ANTP-NKL", 
		"Nk2.1", "ANTP-NKL", 
		"Nk2.2", "ANTP-NKL", 
		"Nk3", "ANTP-NKL", 
		"Nk5-Hmx", "ANTP-NKL", 
		"Nk6", "ANTP-NKL", 
		"Nk7", "ANTP-NKL", 
		"Noto", "ANTP-NKL", 
		"Ro", "ANTP-NKL", 
		"Tlx", "ANTP-NKL", 
		"Vax", "ANTP-NKL", 
		"Ventx", "ANTP-NKL",
	);
	my $assoc = $keys{$family};
	return $assoc;
}
exit;