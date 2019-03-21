#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::MoreUtils qw/ uniq /;
use Array::Utils qw(:all);
use Cwd;
#./interproscan.sh -appl <pickprogram> -i /home/cristig/Desktop/Homeodomains/All.fasta -b /home/cristig/Desktop/Homeodomains/all.interpro
my @remove = qw (hsap
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
print "Analaysing top $HG blasts...\n\n";
my @topblast = blast_determination($blastfile,$HG);
print "Top $HG blasts analysed.\n\nAnalysing top $HG interpros...\n\n";
my @topinterpro = interpro_determination($interprofile,$HG);
print "Top $HG interpros analysed.\n\nComparing 4 types of analysis results...\n\n";
my %classifieds = ();
my %topblast = ();
my %topinterpro = ();
my %allheaders = ();
my @homeoclasses = `grep -v "Conserved Motifs" HomeoDBClass.csv`;
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
my @fastaheader = `grep ">" $fastafile`;
foreach my $header (@fastaheader) {
	chomp $header;
	$header =~ s/>//g;
	my $gene = "Unknown";
	#my ($homology,$gene2,$species,$count) = split /_/, $header;
	my ($species,$origheader,$homology,$count) = split /|/, $header;
	my $specstring = `grep -h ',$species,' phylogenyTable.csv`;
	my $family = "Unknown";
	my $class = "Unknown";
	my $group = "Unknown";
	my @specarray = split /,/, $specstring;
	if ($specstring =~ /Metazoa/) {
		$group = $specarray[11];
	}
	else {
		$group = $specarray[6];
	}
	my $homeoclasscsv = "HomeoDBClass.csv";
	if ($specstring =~ /Lophotrochozoa/) {
		$homeoclasscsv = "HomeoboxClassLopho.csv";
	}
	#if ($specstring =~ m/(Uniprot|EnsEMBL)/) {
	#	my $orig = `grep '$header' fastaCollection.log`;
	#	if ($orig) {
	#		my @orig = split / /, $orig;
	#		if ($orig[0] =~ /^$species\_(.*)\_ENS/) {
	#			$gene = $1;
	#		}
	#		elsif ($orig[0] =~ /^$species\_tr.*GN=([A-Za-z0-9]+)_.*/) {
	#			$gene = $1;
	#		}
	#	}
	#	else {
	#		$gene = "-";
	#	}
	#	if ($gene =~ /([A-Za-z]+)/) {
	#		my $pattern = build_partial($gene,3);
	#		if (my $homeo = `grep -E -h "$pattern" $homeoclasscsv`) {
	#			($class, $family, my $motifs, my $end) = split /,/, $homeo;
	#			$family = $gene;
	#		}
	#		else {
	#			if ($specstring !~ /Metazoa/) {
	#				if (exists $topinterpro{$header}) {
	#					$topinterpro{$header} =~ s/ domain//g;
	#					$topinterpro{$header} =~ s/Transcription factor//g;
	#					$topinterpro{$header} =~ s/,//g;
	#					$topinterpro{$header} =~ s/superfamily//g;
	#					$topinterpro{$header} =~ s/family//g;
	#					$pattern = build_partial($topinterpro{$header},3);
	#					if (my $homeo = `grep -E -h "$pattern" $homeoclasscsv`) {
	#						($class, $family, my $motifs) = split /,/, $homeo;
	#					}
	#				}	
	#			}
	#			elsif (exists $topblast{$header}) {
	#				(my $blastspec, my $subfamily, $family, $class) = split /\|/, $topblast{$header};
	#				if (my $homeo = `grep -E -h "$family" $homeoclasscsv`) {
	#					($class, $family, my $motifs) = split /,/, $homeo;
	#				}
	#			}
	#			else {
	#				if (exists $topinterpro{$header}) {
	#					$topinterpro{$header} =~ s/ domain//g;
	#					$topinterpro{$header} =~ s/Transcription factor//g;
	#					$topinterpro{$header} =~ s/,//g;
	#					$topinterpro{$header} =~ s/superfamily//g;
	#					$topinterpro{$header} =~ s/family//g;
	#					$pattern = build_partial($topinterpro{$header},3);
	#					$family = $gene;
	#					if (my $homeo = `grep -E -h "$pattern" $homeoclasscsv`) {
	#						($class, my $family2, my $motifs) = split /,/, $homeo;
	#					}
	#				}	
	#			}
	#		}
	#	}
	#	else {
	#		$gene = "Unknown";
	#		$family = "Unknown";
	#	}
	#}
	unless ( grep( /^$species$/, @remove ) ) {
		my $box = header_match($origheader);
		if ($box != 0) {
			my $homeo = `grep -E "$box" $homeoclasscsv`;
			($class, $family, my $motifs, $gene) = split /,/, $homeo;
		}
		else {
			$gene = "Unknown";
			if ($specstring !~ /Metazoa/) {
				if (exists $topinterpro{$header}) {
					$topinterpro{$header} =~ s/ domain//g;
					$topinterpro{$header} =~ s/Transcription factor//g;
					$topinterpro{$header} =~ s/,//g;
					$topinterpro{$header} =~ s/superfamily//g;
					$topinterpro{$header} =~ s/family//g;
					my $pattern = build_partial($topinterpro{$header},3);
					if (my $homeo = `grep -E -h "$pattern" $homeoclasscsv`) {
						($class, $family, my $motifs, $gene) = split /,/, $homeo;
					}
				}	
			}
			elsif (exists $topblast{$header}) {
				(my $blastspec, my $subfamily, $family, $class) = split /\|/, $topblast{$header};
				if (my $homeo = `grep -E -h "$family" $homeoclasscsv`) {
					($class, $family, my $motifs) = split /,/, $homeo;
				}
			}
			else {
				if (exists $topinterpro{$header}) {
					$topinterpro{$header} =~ s/ domain//g;
					$topinterpro{$header} =~ s/Transcription factor//g;
					$topinterpro{$header} =~ s/,//g;
					$topinterpro{$header} =~ s/superfamily//g;
					$topinterpro{$header} =~ s/family//g;
					my $pattern = build_partial($topinterpro{$header},3);
					if (my $homeo = `grep -E -h "$pattern" $homeoclasscsv`) {
						($class, $family, my $motifs) = split /,/, $homeo;
					}
				}	
			}
		}
		my $newheader =  "$species|$group|$family|$class|$homology|$count";
		$newheader =~ s{\\}{-}g;
		$classifieds{$header} = $newheader;
	}
}
print "Almost done...\n\n";
unless (-d "LOG") {
	mkdir "LOG";
}
open(LOG,">LOG/classifier.$HG.log");
while (my ($keys,$values) = each %classifieds) {
	$values =~ s/Holomycota\/Nucletmycea/Holomycota/;
	$keys = quotemeta($keys);
	$values = quotemeta($values);
	system "sed -i 's/$keys/$values/g' FastasAnnot/$HG.fasta";
	if (-f "fastaCollection.log") {
		my $origline = `grep '$keys' fastaCollection.log`;
		chomp $origline;
		print LOG "$values -> $origline\n"; 
	}
	else {
		print LOG "$keys\t$values\n";
	}
}
close LOG;

sub build_partial {
    my ($str, $min) = (@_, 1);
	$str = ",$str";
    my @re;
    for (0 .. length($str) - $min) {
        my $str = substr $str, $_;
        for ($min .. length $str) {
            push @re, quotemeta substr $str, 0, $_
        }
    }
    my $re = join '|' => sort {length $a <=> length $b} @re;
    return ("($re)");
}
##Sub to parse Blast file and return just the top result per protein
sub blast_determination {
	my $blastfile = $_[0];
	my $HG = $_[1];
	my @entries = `grep "^$HG\_" $blastfile`;
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
	print "3. Blasts top hits Determined\n\n";
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
	my @entries = `grep -h "^$hg\_" $interprofile`;
	my %bestentry = ();
	foreach my $result (@entries) {
		chomp $result;	#query\trubbish\tlength\tanalysis\taccession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tipid\tannotation
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
	my @bestentries = ();
	foreach my $query ( sort keys %bestentry ) {
		push @bestentries, "$query\t$bestentry{$query}[0]";
	}
	print "4. Interproscan top results determined\n\n";
	return @bestentries;
	
}
sub header_match {
	my $header = shift;
	my @boxes = qw (ZHX1 ZHX2 ZHX3 Abox Barhl Barx Bsx Cdx Dbx Dlx Emx En Evx Gbx Gsx Hhex Hlx Hox1 Hox2 Hox3 Hox4 Hox5 Hox6-8 Hox9-13(15) Lbx Meox Mnx Nanog Nedx Nk1 Nk2.1 Nk2.2 Nk3 Nk4 Nk5/Hmx Nk6 Nk7 Noto Pdx Tlx Vax Ventx Cux Onecut Satb Hmbox Hnf1 Isl Lhx1/5 Lhx2/9 Lhx3/4 Lhx6/8 Lmx Bix NANOGNB Sia Hdx Pou1 Pou2 Pou3 Pou4 Pou5 Pou6 Pou2 Alx Argfx Arx Dmbx Dprx Drgx Dux Esx Gsc Hesx Hopx Isx Leutx Mix Nobox Otp Pax3/7 Pax4/6 Phox Pitx Prop Prrx Rax Rhox Sebox Shox Tprx Uncx Vsx Prox Six1/2 Six3/6 Six4/5 Irx Meis Mkx Pbx Pknox Tgif Adnp Tshz Zeb Zfhx Zhx/Homez);
	my $match = 0;
	foreach my $box (@boxes) {
		if ($header =~ /$box/i) {
			$match = $box;
		}
	}
	return $match;
}
exit;