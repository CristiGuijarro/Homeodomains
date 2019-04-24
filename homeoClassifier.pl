#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::MoreUtils qw/ uniq /;
use Array::Utils qw(:all);
use Cwd;
#./interproscan.sh -appl <pickprogram> -i All.fasta -b all.interpro
my @remove = qw (
	hsap
	ggal
	xtro
	drer
	bflo
	dmel
	tcas
	cele
);
my $fastafile = "";
my $blastfile = "";
my $interprofile = "";
my $raxmlfile = "";
my $fastreefile = "";

GetOptions ('fastafile=s' => \$fastafile, 'blastfile=s' => \$blastfile, 'interprofile=s' => \$interprofile, 'raxmlfile=s' => \$raxmlfile, 'fastreefile=s' => \$fastreefile) or die "--fastafile needs to be specified, --blastfile and --interprofile will need to be specified also";
chomp $fastafile;
chomp $blastfile;
chomp $interprofile;
my $HG;
if ($fastafile =~ /\/([0-9]+)\.fasta$/) {
	$HG = $1;
}
#print "Analaysing top $HG blasts...\n\n";
my @topblast = blast_determination($blastfile,$HG);
#print "Top $HG blasts analysed.\n\nAnalysing top $HG interpros...\n\n";
my @topinterpro = interpro_determination($interprofile,$HG);
#print "Top $HG interpros analysed.\n\nComparing 4 types of analysis results...\n\n";
my %classifieds = ();
my %topblast = ();
my %topinterpro = ();
my %allheaders = ();
my @homeoclasses = `grep "" HomeoDBClass.csv`;
#### When comparing blast and interpro - interpro to override any result for any organism not metazoan.
foreach my $queries (@topblast) {
	chomp $queries;
	my ($query, $match) = split /\t/, $queries;
	$topblast{$query} = $match;
}
foreach my $queries (@topinterpro) {
	chomp $queries;
	my ($query, $match) = split /\t/, $queries;
	$topinterpro{$query} = $match;
}
my %header = ();
my %interpro = ();
my %blasthomeo = ();
my %determined = ();
my @fastaheader = `grep ">" $fastafile`;
foreach my $header (@fastaheader) {
	chomp $header;
	$header =~ s/>//g;
	my $gene = "Unknown";
	my ($species,$origheader,$homology,$count) = split /\|/, $header;
	chomp $count;
	my $uid = "$homology|$count";
	my $specstring = `grep -h ',$species,' phylogenyTable.csv`;
	my $family = "Unknown";
	my $class = "Unknown";
	my $group = "Unknown";
	my @specarray = split /,/, $specstring;
	if ($specstring =~ /Metazoa/) {
		$group = phylum_select($species);
		chomp $group;
	}
	else {
		$group = $specarray[6];
	}
	my $homeoclasscsv = "HomeoDBClass.csv";
	unless ( grep( /^$species$/i, @remove ) ) {
		my $box = "0";
		unless ($origheader =~ /Uncharacterized_protein/) {
			$box = header_match($origheader);
		}
		if ($box ne "0") {
			$box =~ s/[\_\-\.\(\)\\\/]/\./g;
			#print "1)$group\t$box\n";
			#if (my $homeo = `grep -E -h -m1 -i '$box' $homeoclasscsv`) {
			if (my ($homeo) = grep/$box/i, @homeoclasses) {
				($class, $family, my $motifs, $gene) = split /,/, $homeo;
				$header{$uid} = "$box";
			}
			else {
				$header{$uid} = "$box";
				#print "grep -E -h -m1 -i '$box' $homeoclasscsv\n";
			}
		}
		else {
			$header{$uid} = "-";
			$gene = "Unknown";
			if ($specstring !~ /Metazoa/) {
				if (exists $topinterpro{$header}) {
					$interpro{$uid} = $topinterpro{$header};
					my $pattern = header_match($topinterpro{$header});
					if ($pattern ne "0") {
						#if (my $homeo = `grep -E -h -m1 -i "$pattern" $homeoclasscsv`) {
						if (my ($homeo) = grep/$pattern/i, @homeoclasses) {
							($class, $family, my $motifs, $gene) = split /,/, $homeo;
						}
						#print "2)$group\t$pattern\n";
					}
					else {
						$pattern = build_partial($topinterpro{$header},6);
						#print "3)$group\t$pattern\n";
					}
				}
				else {
					$interpro{$uid} = "-";
				}
			}
			elsif (exists $topblast{$header}) {
				$blasthomeo{$uid} = $topblast{$header};
				(my $blastspec, my $subfamily, $family, $class) = split /\|/, $topblast{$header};
				#if (my $homeo = `grep -E -h -m1 -i "$family" $homeoclasscsv`) {
				if (my ($homeo) = grep/$family/i, @homeoclasses) {
					($class, $family, my $desc, my $motifs) = split /,/, $homeo;
					#print "4)$group\t$family\n";
				}
			}
			else {
				$blasthomeo{$uid} = "-";
				if (exists $topinterpro{$header}) {
					my $box = $topinterpro{$header};
					#$box =~ s/ //g;
					$interpro{$uid} = $topinterpro{$header};
					my $pattern = header_match($box);
					if ($pattern ne "0") {
						#if (my $homeo = `grep -E -h -m1 -i "$pattern" $homeoclasscsv`) {
						if (my ($homeo) = grep/$pattern/i, @homeoclasses) {
							($class, $family, my $motifs, $gene) = split /,/, $homeo;
							#print "5)".$topinterpro{$header}."\t$group\t$pattern\n";
						}
					}
				}
				else {
					$interpro{$uid} = "-";
				}
			}
		}
		my $newheader =  "$species|$group|$family|$class|$homology|$count";
		$newheader =~ s{\\}{-}g;
		$classifieds{$header} = $newheader;
		$determined{$uid} = "$family|$class";
	}
	unless (exists $blasthomeo{$uid}) {
		$blasthomeo{$uid} = "-";
	}
	unless (exists $interpro{$uid}) {
		$interpro{$uid} = "-";
	}
	unless (exists $header{$uid}) {
		$header{$uid} = "-";
	}
	if (exists $determined{$uid}) {
		if ($header{$uid} =~ /Unknown/) {
			$header{$uid} = "-";
		}
		print "$uid\t",
		$header{$uid}, "\t",
		$blasthomeo{$uid}, "\t",
		$interpro{$uid}, "\t",
		$determined{$uid}, "\tPLACEHOLDER#\n";
	}
}
#print "Almost done...\n\n";
unless (-d "LOG") {
	mkdir "LOG";
}
unless (-d "FastasAnnot") {
	mkdir "FastasAnnot";
}
open(LOG,">LOG/classifier.$HG.log");
open(FASTA, ">FastasAnnot/$HG.fasta");
while (my ($keys,$values) = each %classifieds) {
	$values =~ s/Holomycota\/Nucletmycea/Holomycota/;
	$values =~ s/Viridiplantae/Plants/;
	my $sequence = `grep -A1 '$keys' $fastafile`;
	my ($header, $seq) = split /\n/, $sequence;
	chomp $seq;
	chomp $values;
	print FASTA ">$values\n$seq\n";
	print LOG "$keys\t$values\n";
}
close LOG;
close FASTA;

#Sub to build partial patterns for pattern matching
sub build_partial {
    my $str = shift;
	my $min = shift;
	$str =~ s/[^a-zA-Z0-9,]/\./g;
	$str =~ s/Homeobox//ig;
    my @re;
    for (0 .. length($str) - $min) {
        my $str = substr $str, $_;
        for ($min .. length $str) {
            push @re, substr $str, 0, $_
        }
    }
    my $re = join '|' => sort {length $a <=> length $b} @re;
    return ("($re)");
}
##Sub to parse Blast file and return just the top result per protein
sub blast_determination {
	my $blastfile = $_[0];
	my $HG = $_[1];
	my @entries = `grep -E "\|$HG\|" $blastfile`;
	my %bestentry = ();
	foreach my $result (@entries) {
		chomp $result;
		my ($query, $match, $id, $length, $mismatches, $gapopen, $qstart, $qend, $mstart, $mend, $evalue, $bitscore) = split /\t/, $result;
		unless (exists $bestentry{"$query"}) {
			$bestentry{"$query"} = ['match', '0', '0', '100', '10', '0', '0'];
		}
		my $overallscore = 0;
		my @queryheader = split /_/, $query;
		my $queryspec = $queryheader[0];
		my ($matchspec, $family, $subclass, $class) = split /|/, $match;
		my $closestspec = organism_score($queryspec,$matchspec);
		$overallscore +=$closestspec;
		if ($id >= $bestentry{$query}[1]) {
			$overallscore+=1;
		}
		if ($length >= $bestentry{$query}[2]) {
			$overallscore+=1;
		}
		if ($mismatches <= $bestentry{$query}[3]) {
			$overallscore+=1;
		}
		if ($evalue <= $bestentry{$query}[4]) {
			$overallscore+=1;
		}
		if ($bitscore >= $bestentry{$query}[5]) {
			$overallscore+=1;
		}
		if ($overallscore >= $bestentry{$query}[6]) {
			$bestentry{"$query"} = ["$match", "$id", "$length", "$mismatches", "$evalue", "$bitscore", "$overallscore"];
		}
	}
	my @bestentries = ();
	foreach my $query ( sort keys %bestentry ) {
		push @bestentries, "$query\t$bestentry{$query}[0]";
	}
	#print "3. Blasts top hits Determined\n\n";
	return @bestentries;
}
###Sub in sub to determine most reasonable result based on organism as well as top hit
sub organism_score {
	my $queryspec = $_[0];
	my $matchspec = $_[1];
	my %group = ("Nematode" => ["Ecdysozoa"], "Beetle" => ["Ecdysozoa"], "Human" => ["Deuterostomia"], "Amphioxus" => ["Cephalochordata"]);
	my $phyloline = `grep -h ",$queryspec," phylogenyTable.csv`;
	my $groupvalue = $group{$matchspec};
	if (exists $group{$matchspec}) {
		if ($phyloline =~ /$groupvalue/) {
		return 1;
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
}
##Sub to parse Interpro file and return top result per protein
sub interpro_determination {
	my $interprofile = $_[0];
	my $hg = $_[1];
	my @entries = `grep -h -E "\|$hg\|" $interprofile`;
	my %bestentry = ();
	foreach my $result (@entries) {
		chomp $result;
		if ($result !~ /\tHomeobox domain\t/) {
			if ($result !~ /\tHomeobox signature\t/) {
				my ($query,$rubbish,$length,$analysis,$accession,$description,$start,$stop,$evalue,$status,$date,$tipid,$annotation) = split /\t/, $result;
				unless(defined $annotation) {
					next;
				}
				unless(defined $evalue) {
					next;
				}
				unless (exists $bestentry{"$query"}) {
					$bestentry{"$query"} = ['annotation', '0', '10', '0'];
				}
				my $overallscore = 0;
				my $matchlength = $stop-$start;
				if ($matchlength >= $bestentry{$query}[1]) {
					$overallscore+=1;
				}
				if ($evalue <= $bestentry{$query}[2]) {
					$overallscore+=2;
				}
				if ($overallscore >= $bestentry{$query}[3]) {
					$bestentry{"$query"} = ["$annotation", "$matchlength", "$evalue", "$overallscore"];
				}
			}
		}
	}
	my @bestentries = ();
	foreach my $query ( sort keys %bestentry ) {
		push @bestentries, "$query\t$bestentry{$query}[0]";
	}
	#print "4. Interproscan top results determined\n\n";
	return @bestentries;
	
}
#Sub to match headers and pattern match for detected homeodomains
sub header_match {
	my $header = shift;
	my $homeoclasscsv = "HomeoDBClass.csv";
	my @homeoclasses = `grep "" $homeoclasscsv`;
	my @boxes = qw (ZHX1 ZHX2 ZHX3 Barhl Barx Bsx Cdx Dbx Dlx Emx En1 En2 Evx Gbx Gsx Hhex Hlx Hox1 Hox2 Hox3 Hox4 Hox5 Hox6/8 Hox9/15 Lbx Meox Mnx Nanog Nedx Nk1 Nk2 Nk2.1 Nk2.2 Nk3 Nk4 Nk5/Hmx Nk6 Nk7 Noto Pdx Tlx Vax Ventx Cux Onecut Satb Hmbox Hnf1 Isl Lhx1/5 Lhx2/9 Lhx3/4 Lhx6/8 Lmx Bix NANOGNB Sia Hdx Pou1 Pou2 Pou3 Pou4 Pou5 Pou6 Pou2 Alx Argfx Arx Dmbx Dprx Drgx Dux Esx Gsc Hesx Hopx Isx Leutx Mix Nobox Otp Pax3/7 Pax4/6 Phox Pitx Prop Prrx Rax Rhox Sebox Shox Tprx Uncx Vsx Prox Six1/2 Six3/6 Six4/5 Irx Meis Mkx Pbx Pknox pknox/meis Tgif Adnp Tshz Zeb Zfhx Zhx/Homez Engrailed Caudal Developing Distal-less Empty-Spiracles Even-skipped hematopoietically Gastrulation Mesenchyme Pancreas Segment Pseudogene Pancreatic T-Cell-Leukemia Pre-B-Cell-Leukemia Ventral Cut-like One-cut Brachyury-inducible Siamois Highly-divergent Aristalless Arginine-fifty Diencephalon/mesencephalon Divergent-paired Dorsal Double Goosecoid Intestine Pituitary Retina Rentinal Stature Visual Prospero Oculis Iroquois Mohawk Knotted TGFB teashirt-zinc Zinc-finger-protein Zinc-finger-E-box Zinc-finger-plus zf ZAG-1 Zinc-finger-homeobox-potein zinc-finger-E-box POU-specific LIM-type rx1/3 Absent-from-Olfactores);
	my $match = 0;
	unless ($header =~ /Uncharacterized_protein/) {
		foreach my $box (@boxes) {
			chomp $box;
			$box =~ s/ //g;
			if ($box =~ /(^[a-z]+)([0-9]+$)/i) {
				my $one = $1;
				my $two = $2;
				if ($header =~ /$one.$two./i) {
					$match = "$box";
				}	
			}
			if ($box =~ /(^[a-z]+)([0-9]+)\/([0-9]+$)/i) {
				my $one = $1;
				my $two = $2;
				my $three = $3;
				for (my $i=$two; $i<($three+1); $i++) {
					if ($header =~ /$one.$i./i) {
						$match = "$box";
					}
				}
			}
			if ($box =~ /(^[a-z]+)\/([a-z]+$)/i) {
				my $one = $1;
				my $two = $2;
				if ($header =~ /$one./i) {
					$match = $box;
				}
				if ($header =~ /$two./i) {
					$match = $box;
				}
			}
			if ($box =~ /(^[a-z]+)-([a-z]+)-([a-z]+)-([a-z]+$)/i) {
				my $one = $1;
				my $two = $2;
				my $three = $3;
				my $four = $4;
				if ($header =~ /$one.$two.$three.$four./i) {
					$match = $box;
				}
			}
			if ($box =~ /(^[a-z]+)\-([a-z]+)\-([a-z]+$)/i) {
				my $one = $1;
				my $two = $2;
				my $three = $3;
				if ($header =~ /$one.$two.$three./i) {
					$match = $box;
				}
			}
			if ($box =~ /(^[a-z]+)\-([a-z]+$)/i) {
				my $one = $1;
				my $two = $2;
				if ($header =~ /$one.$two./i) {
					$match = $box;
				}
			}
		}
		#print "$match-1\t";
		if ($match eq "0") {
			if ($header =~ /^(.*)\_ENS/) {
					$match = $1;
			}
			elsif ($header =~ /\_(.*)\=/) {
				$header = $1;
			}
			else {
				my $gene = "0";
				if (length($header) > 50) {
					if ($header =~ /tr_(.*)\_GN=(.*)\_/) {
						$header = $1;
						$gene = $2;
					}
					elsif ($header =~ /tr_(.*)=/) {
						$header = $1;
					}
				}
				if ($match eq "0") {
					if ($gene =~ /[A-Z]/) {
						#if (my $homeo = `grep -m1 -E -h -i "$gene" $homeoclasscsv`) {
						if (my ($homeo) = grep/$gene/i, @homeoclasses) {
							#print "\t$homeo-1\t";
							(my $class, my $family, my $desc, my $genes) = split /,/, $homeo;
							$match = $family;
						}
					}
					#else {
					#	my $pattern = build_partial($header,6);
					#	#if (my $homeo = `grep -m1 -E -h -i "$pattern" $homeoclasscsv`) {
					#	if (my ($homeo) = grep/$pattern/i, @homeoclasses) {
					#		print "\t$homeo-2\t";
					#		(my $class, my $family, my $desc, my $genes) = split /,/, $homeo;
					#		$match = $family;
					#	}
					#}
				}
			}
		}
	}
	$match =~ s/\-/\./g;
	if ($match !~ /[A-Z]/i) {
		$match = 0;
	}
	return $match;
}
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