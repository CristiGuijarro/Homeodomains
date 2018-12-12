#!/usr/bin/perl -w
use strict;

my $annotfastadir = $ARGV[0];
chomp $annotfastadir;
my @fastas = `grep -h '>' $annotfastadir/*.fasta`;
my %hash = ();
my %hbxs = ();
my %animals = ();
#my @uniani = ("Leucosolenia complicata","Sycon ciliatum","Amphimedon queenslandica","Oscarella carmela","Patiria miniata","Strongylocentrotus purpuratus","Lytechinus variegatus","Ptychodera flava","Saccoglossus kowalevskii","Branchiostoma floridae","Branchiostoma belcheri","Oikopleura dioica","Ciona intestinalis","Botryllus schlosseri","Danio rerio","Homo sapiens","Anolis carolinensis","Gallus gallus","Xenopus tropicalis","Trichinella spiralis","Brugia malayi","Romanomermis culicivorax","Caenorhabditis elegans","Stegodyphus mimosarum","Acanthoscurria geniculata","Ixodes scapularis","Mesobuthus martensii","Limulus polyphemus","Daphnia pulex","Parhyale hawaiensis","Tribolium castaneum","Drosophila melanogaster","Zootermopsis nevadensis","Strigamia maritima","Hypsibius dujardini","Ramazzottius varieornatus","Helobdella robusta","Capitella teleta","Lingula anatina","Crassostrea gigas","Pinctada fucata","Octopus bimaculoides","Lottia gigantea","Intoshia linei","Hymenolepis microstoma","Echinococcus multilocularis","Gyrodactylus salaris","Schmidtea mediterranea","Macrostomum lignano","Schistosoma japonicum","Adineta vaga","Aiptasia pallida","Nematostella vectensis","Acropora digitifera","Hydra magnipapillata","Thelohanellus kitauei","Pleurobrachia bachei","Mnemiopsis leidyi","Trichoplax adhaerens");
#my @unigroup = ("Calcarea","Demospongiae","Homoscleromorpha","Ambulacraria","Chordata","Ecdysozoa","Lophotrochozoa","Anthozoa","Hydrozoa","Myxozoa","Tentaculata","Trichoplacidae");
#($spec,$kingdom,$family,$class,$HG,$uid) = split /\|/, $header;
foreach my $fasta(@fastas) {
	chomp $fasta;
	$fasta =~ s/^>//;
	if ($fasta =~ /(^[A-Za-z0-9]{4})_/) {
		$fasta = "$1|Unknown|Unknown|Unknown|0|0";
	}
	if ($fasta =~ /\|.*\|.*\|.*\|.*\|/) {
		#unless ($fasta =~ /[(Unknown)(homeobox)(unassigned)]/) {
			my ($spec,$group,$family,$class,$Homology,$uid) = split /\|/, $fasta;
			my $lookup = `grep "$spec" phylogenyTable.csv`;
			my @parts = split /,/, $lookup;
			my $six = $parts[11];
			my $kingdom = $parts[6];
			if ($kingdom eq "Metazoa") {
				$spec = $parts[0];
				if ($class =~ /ANTP/) {
					$class = "ANTP";
				}
				$class = uc($class);
				$family = lc($family);
				$family = ucfirst($family);
				$hash{"$spec,$kingdom,$six,$class,$family"}++;
				$animals{"$spec,$kingdom,$six"}++;
				if ($kingdom eq "Metazoa") {
					$hbxs{"$class,$family"}++;
				}
			#}	
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
print " , ,Homeobox,";
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
exit;