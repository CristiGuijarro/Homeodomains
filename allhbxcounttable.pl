#!/usr/bin/perl -w
use strict;

my $annotfastadir = $ARGV[0];
chomp $annotfastadir;
my @fastas = `grep -h '>' $annotfastadir/*.fasta`;
my %hash = ();
my %hbxs = ();
my %animals = ();
my @homeoanimals = qw (
	Human
	Beetle
	Amphioxus
	Frog
	Chicken 
	Fruitfly
	Zebrafish
	Nematode
);
#my @uniani = ("Leucosolenia complicata","Sycon ciliatum","Amphimedon queenslandica","Oscarella carmela","Patiria miniata","Strongylocentrotus purpuratus","Lytechinus variegatus","Ptychodera flava","Saccoglossus kowalevskii","Branchiostoma floridae","Branchiostoma belcheri","Oikopleura dioica","Ciona intestinalis","Botryllus schlosseri","Danio rerio","Homo sapiens","Anolis carolinensis","Gallus gallus","Xenopus tropicalis","Trichinella spiralis","Brugia malayi","Romanomermis culicivorax","Caenorhabditis elegans","Stegodyphus mimosarum","Acanthoscurria geniculata","Ixodes scapularis","Mesobuthus martensii","Limulus polyphemus","Daphnia pulex","Parhyale hawaiensis","Tribolium castaneum","Drosophila melanogaster","Zootermopsis nevadensis","Strigamia maritima","Hypsibius dujardini","Ramazzottius varieornatus","Helobdella robusta","Capitella teleta","Lingula anatina","Crassostrea gigas","Pinctada fucata","Octopus bimaculoides","Lottia gigantea","Intoshia linei","Hymenolepis microstoma","Echinococcus multilocularis","Gyrodactylus salaris","Schmidtea mediterranea","Macrostomum lignano","Schistosoma japonicum","Adineta vaga","Aiptasia pallida","Nematostella vectensis","Acropora digitifera","Hydra magnipapillata","Thelohanellus kitauei","Pleurobrachia bachei","Mnemiopsis leidyi","Trichoplax adhaerens");
#my @unigroup = ("Calcarea","Demospongiae","Homoscleromorpha","Ambulacraria","Chordata","Ecdysozoa","Lophotrochozoa","Anthozoa","Hydrozoa","Myxozoa","Tentaculata","Trichoplacidae");
#($spec,$kingdom,$family,$class,$HG,$uid) = split /\|/, $header;
foreach my $fasta(@fastas) {
	chomp $fasta;
	$fasta =~ s/^>//;
	$fasta =~ s/ //g;
	#if ($fasta =~ /(^[A-Za-z0-9]{4})_/) {
	#	$fasta = "$1|Unknown|Unknown|Unknown|0|0";
	#}
	my ($species,$group,$family,$class,$Homology,$uid) = "";
	if ($fasta =~ /.*_.*_.*_.*_.*_.*/) {
		($species,$group,$family,$class,$Homology,$uid) = split /_/, $fasta;
		#print "$fasta $genus\n";
	}
	elsif ($fasta =~ /.*\|.*\|.*\|.*\|HomeoDB\|[0-9]+/) {
		#Amphioxus|Cephalochordata|Barx|ANTP|HomeoDB|77
		($species,$group,$family,$class,$Homology,$uid) = split /\|/, $fasta;
		#print "$fasta $species\n";
	}
	else {
		print "$fasta\n";
	}
	my $kingdom;
	my $six;
	my $spec;
	if (grep( /$species/, @homeoanimals )) {
		$six = $group;
		$kingdom = "Metazoa";
		$spec = $species;
	}
	else {
		my $lookup = `grep -i -h ',$species,' phylogenyTable.csv`;
		#print $species, "\n";
		my @parts = split /,/, $lookup;
		$six = $parts[11];
		$kingdom = $parts[6];
		$spec = $parts[1];
		if ($kingdom eq "Metazoa") {
			$six = phylum_select($spec);
			chomp $six;
		}
		$species = $parts[0];
	}
	if ($kingdom eq "Metazoa" || grep( /^$species$/, @homeoanimals ) ) {
		$class = uc($class);
		$family = lc($family);
		$family = ucfirst($family);
		#print "$fasta\n";
		$hash{"$species,$kingdom,$six,$class,$family"}++;
		$animals{"$species,$kingdom,$six"}++;
		if ($kingdom eq "Metazoa") {
			$hbxs{"$class,$family"}++;
		}
	}
}
my @hbx = ();
my @animals = ();
foreach my $keys3 (sort keys %hbxs) {
	push @hbx, $keys3;
	foreach my $keys2 (sort keys %animals) {
		unless (exists($hash{"$keys2,$keys3"})) {
		$hash{"$keys2,$keys3"} = 0;
		}
	}
}
foreach my $keys2 (sort keys %animals) {
	push @animals, $keys2;
}
my @newlines = (@animals);
foreach (my $i = 0; $i < @animals; $i++) {
	foreach my $hbx (@hbx) {
		$newlines[$i].=",".$hash{"$animals[$i],$hbx"};
	}
}
print ",,Homeobox,";
my @class = ();
my @family = ();
foreach my $hbx (@hbx) {
	my ($class, $family) = split /,/, $hbx;
	push @class, $class;
	push @family, $family;
}
print join(",", @class), "\n";
print "Species,Kingdom,Clade,";
print join(",", @family), "\n";
print join("\n", @newlines), "\n";

sub phylum_select {
	my $spec = shift;
	chomp $spec;
	my %keys = (
		drer => "Vertebrata",
		xtro => "Vertebrata",
		acar => "Vertebrata",
		ggal => "Vertebrata",
		hsap => "Vertebrata",
		odio => "Urochordata",
		cint => "Urochordata",
		bsch => "Urochordata",
		bflo => "Cephalochordata",
		bbel => "Cephalochordata",
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
exit;