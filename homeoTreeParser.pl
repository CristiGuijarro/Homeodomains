#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
use List::Util qw(max);
use List::UtilsBy qw(max_by);

unless (-d "TreeFiles") {
	mkdir "TreeFiles";
}
my $dir = $ARGV[0]; 
my $tsv = $ARGV[1]; #classificationTable.tsv as written in homeoClassifier.
chomp $tsv;
my $classlog = `grep '^[0-9]' $tsv`;
my @classlog = split /\n/, $classlog;
my @files = `ls $dir/*.contree*`;
foreach my $file (@files) {
	chomp $file;
	my $HG = "";
	if ($file =~ /\/([0-9]+)\.*/) {
		$HG = $1;
	}
	open(TREE,">TreeFiles/$HG.newick");
	my @Leaves = ();
	my %Leaves = ();
	my $treeio = Bio::TreeIO->new(-format => 'newick', -file => $file, -internal_node_id => 'bootstrap');
	my %homeohash = ();
	my @homeoclass = ();
	my %reducedhash = ();
	my @reducedclass = ();
	my %replacements = ();
	while( my $tree = $treeio->next_tree ) {
		$tree->move_id_to_bootstrap;
		my $rootnode = $tree->get_root_node;
		my @nodes = grep { $_->id =~ /[A-Z]+/ } $tree->get_leaf_nodes;
		for my $leaf (@nodes) {
			if ($leaf->id =~ /([A-Za-z0-9]+)_([0-9]+)'?$/) {
				my $homology = $1;
				my $count = $2;
				my $uid = "$homology|$count";
				$Leaves{$uid} = $leaf->id;
			}
			elsif ($leaf->id =~ /[a-z]+$/) {
				$Leaves{$leaf} = $leaf->id;
			}
		}
		for my $node (@nodes) {
			my $uid = "";
			if ($node->id =~ /^'?[A-Za-z0-9]+_[A-Za-z]+_Unknown_Unknown_([0-9]+)_([0-9]+)'?$/) {
				my $homology = $1;
				my $count = $2;
				$uid = "$homology|$count";
				$Leaves{$uid} = $node->id;
				my @sequid = ( grep { $classlog[$_] =~ /\A$homology\|$count\t/ } 0..$#classlog );
				($uid, my $header, my $blast, my $interpro, my $determined, my $placeholder) = split /\t/, $classlog[$sequid[0]];
				$classlog[$sequid[0]] = "$uid\t$header\t$blast\t$interpro\t$determined\t$placeholder";
				my $class = "0";
				my %minihash = ();
				my $ancestor = $node->ancestor;
				for my $child ( $ancestor->get_all_Descendents ) {
					if ($child->id =~ /^'?[A-Z0-9]+_[A-Z]+_(.*)_HomeoDB_[0-9]+'?$/i) {
						$class = $1;
						unless ($class eq "0") {
							$minihash{$class}++;
						}
					}
					#elsif ($child->id =~ /^'?[A-Za-z0-9]+_[A-Za-z]+_(.*_[A-Z]+)_[0-9]+_[0-9]+'?$/) {
					#	my $class = $1;
					#	unless ($class eq "0") {
					#		$minihash{$class}++;
					#	}
					#}
				}
				while ($class eq "0") {
					my $counter = 0;
					#while ($counter < 4) {
						$ancestor = $ancestor->ancestor;
						for my $child ( $ancestor->get_all_Descendents ) {
							if ($child->id =~ /^'?[A-Z0-9]+_[A-Z]+_(.*)_HomeoDB_[0-9]+'?$/i) {
								$class = $1;
								$minihash{$class}++;
							}
							#elsif ($child->id =~ /^'?[A-Za-z0-9]+_[A-Za-z]+_(.*_[A-Z]+)_[0-9]+_[0-9]+'?$/) {
							#	$class = $1;
							#	$minihash{$class}++;
							#}
						}
						$counter++;
						if ($counter == 4) {
							#if ($class eq "0") {
								$class = "unassigned_Other";
								$minihash{$class}++;
								$minihash{$class}++;
							#}
						}
					#}
				}
				my $key = max_by { $minihash{$_} } keys %minihash;
				$classlog[$sequid[0]] =~ s/PLACEHOLDER#/$key/e;
				my $seq = $node->id;
				my $replacement = $seq;
				$replacement =~ s/Unknown_Unknown/$key/e;
				$replacements{$seq} = $replacement;
				$Leaves{$uid} = $replacement;
				$uid = "|$uid|";
				#my $logline =  `sed -i '/$uid/a $replacement' LOG/classifier.$HG.log`;
				#print $logline, "\n";
			}
		}
		my $treefile = `grep '' $file`;
		while (my ($keys, $values) = each %replacements) {
			$treefile =~ s/$keys/$values/e;
		}
		print TREE "$treefile";
	}
	close(TREE);
	while (my ($keys, $values) = each %Leaves) {
		chomp $values;
		push @Leaves, $values;
	}
	my $nexus = figtree_nexus("TreeFiles/$HG.newick", @Leaves);
	my $fasta = generate_fasta($HG,@Leaves);
	my $nexfile = "TreeFiles/$HG.nexus";
	open(NEX, ">$nexfile");
	print NEX $nexus;
	close(NEX);
	my $newfasta = "TreeFiles/$HG.CGDB.fasta";
	open(FASTA, ">$newfasta");
	print FASTA $fasta;
}
open(LOG, ">classifiedTable.tsv");
$classlog = join("\n", @classlog);
$classlog =~ s/PLACEHOLDER#/-/g;
print LOG $classlog, "\n";
close(LOG);

sub figtree_nexus {
	my($tree, @leaves) = @_;
	my $phylogeny = `grep '' $tree`;
	chomp $phylogeny;
	while ($phylogeny =~ m/\)([0-9]+):/g) {
		my $bootstrap = $1;
		$phylogeny =~ s/\)$bootstrap\:/\)[&label=$bootstrap]\:/g;#[&label=89]
	}
	my $string = "#NEXUS\nbegin taxa;\ndimensions ntax=".scalar(@leaves).";\ntaxlabels\n";
	foreach my $leaf (@leaves) {
		chomp $leaf;
		if ($leaf =~ /\./) {
			$leaf ="'$leaf'";
		}
		if ($leaf !~ /HomeoDB/) {
			$leaf .= "[&!color=#9824cf]";
		}
		$string .= "\t$leaf\n";
	}
	$string .= ";\nend;\n\nbegin trees;\n\ttree tree_1 = [&R] ";
	$string .= "$phylogeny\nend;\n\n";
	$string .= "begin figtree;\n".
		"\tset appearance.backgroundColorAttribute=\"Default\";\n".
		"\tset appearance.backgroundColour=#ffffff;\n".
		"\tset appearance.branchColorAttribute=\"User selection\";\n".
		"\tset appearance.branchColorGradient=false;\n".
		"\tset appearance.branchLineWidth=1.0;\n".
		"\tset appearance.branchMinLineWidth=0.0;\n".
		"\tset appearance.branchWidthAttribute=\"Fixed\";\n".
		"\tset appearance.foregroundColour=#000000;\n".
		"\tset appearance.hilightingGradient=false;\n".
		"\tset appearance.selectionColour=#2d3680;\n".
		"\tset branchLabels.colorAttribute=\"User selection\";\n".
		"\tset branchLabels.displayAttribute=\"Branch times\";\n".
		"\tset branchLabels.fontName=\"aakar\";\n".
		"\tset branchLabels.fontSize=8;\n".
		"\tset branchLabels.fontStyle=0;\n".
		"\tset branchLabels.isShown=false;\n".
		"\tset branchLabels.significantDigits=4;\n".
		"\tset layout.expansion=0;\n".
		"\tset layout.layoutType=\"RECTILINEAR\";\n".
		"\tset layout.zoom=0;\n".
		"\tset legend.attribute=\"label\";\n".
		"\tset legend.fontSize=10.0;\n".
		"\tset legend.isShown=false;\n".
		"\tset legend.significantDigits=4;\n".
		"\tset nodeBars.barWidth=4.0;\n".
		"\tset nodeBars.displayAttribute=null;\n".
		"\tset nodeBars.isShown=false;\n".
		"\tset nodeLabels.colorAttribute=\"User selection\";\n".
		"\tset nodeLabels.displayAttribute=\"label\";\n".
		"\tset nodeLabels.fontName=\"aakar\";\n".
		"\tset nodeLabels.fontSize=8;\n".
		"\tset nodeLabels.fontStyle=0;\n".
		"\tset nodeLabels.isShown=true;\n".
		"\tset nodeLabels.significantDigits=4;\n".
		"\tset nodeShapeExternal.colourAttribute=\"User selection\";\n".
		"\tset nodeShapeExternal.isShown=false;\n".
		"\tset nodeShapeExternal.minSize=10.0;\n".
		"\tset nodeShapeExternal.scaleType=Width;\n".
		"\tset nodeShapeExternal.shapeType=Circle;\n".
		"\tset nodeShapeExternal.size=4.0;\n".
		"\tset nodeShapeExternal.sizeAttribute=\"Fixed\";\n".
		"\tset nodeShapeInternal.colourAttribute=\"User selection\";\n".
		"\tset nodeShapeInternal.isShown=false;\n".
		"\tset nodeShapeInternal.minSize=10.0;\n".
		"\tset nodeShapeInternal.scaleType=Width;\n".
		"\tset nodeShapeInternal.shapeType=Circle;\n".
		"\tset nodeShapeInternal.size=4.0;\n".
		"\tset nodeShapeInternal.sizeAttribute=\"Fixed\";\n".
		"\tset polarLayout.alignTipLabels=false;\n".
		"\tset polarLayout.angularRange=0;\n".
		"\tset polarLayout.rootAngle=0;\n".
		"\tset polarLayout.rootLength=100;\n".
		"\tset polarLayout.showRoot=true;\n".
		"\tset radialLayout.spread=0.0;\n".
		"\tset rectilinearLayout.alignTipLabels=false;\n".
		"\tset rectilinearLayout.curvature=0;\n".
		"\tset rectilinearLayout.rootLength=100;\n".
		"\tset scale.offsetAge=0.0;\n".
		"\tset scale.rootAge=1.0;\n".
		"\tset scale.scaleFactor=1.0;\n".
		"\tset scale.scaleRoot=false;\n".
		"\tset scaleAxis.automaticScale=true;\n".
		"\tset scaleAxis.fontSize=8.0;\n".
		"\tset scaleAxis.isShown=false;\n".
		"\tset scaleAxis.lineWidth=1.0;\n".
		"\tset scaleAxis.majorTicks=1.0;\n".
		"\tset scaleAxis.minorTicks=0.5;\n".
		"\tset scaleAxis.origin=0.0;\n".
		"\tset scaleAxis.reverseAxis=false;\n".
		"\tset scaleAxis.showGrid=true;\n".
		"\tset scaleBar.automaticScale=true;\n".
		"\tset scaleBar.fontSize=10.0;\n".
		"\tset scaleBar.isShown=false;\n".
		"\tset scaleBar.lineWidth=1.0;\n".
		"\tset scaleBar.scaleRange=2.0;\n".
		"\tset tipLabels.colorAttribute=\"User selection\";\n".
		"\tset tipLabels.displayAttribute=\"Names\";\n".
		"\tset tipLabels.fontName=\"Liberation Mono\";\n".
		"\tset tipLabels.fontSize=10;\n".
		"\tset tipLabels.fontStyle=0;\n".
		"\tset tipLabels.isShown=true;\n".
		"\tset tipLabels.significantDigits=4;\n".
		"\tset trees.order=true;\n".
		"\tset trees.orderType=\"decreasing\";\n".
		"\tset trees.rooting=false;\n".
		"\tset trees.rootingType=\"User Selection\";\n".
		"\tset trees.transform=false;\n".
		"\tset trees.transformType=\"cladogram\";\n".
		"end;\n\n";
		return $string;
}
sub generate_fasta {
	my ($HG,@headers) = @_;
	my $string = "";
	my $fams = `grep ',' HomeoDBfamilies.csv`;
	my %fams = split /[\n,]/, $fams;
	my @fams = ();
	while (my ($keys, $values) = each %fams) {
		push @fams, "$keys\_$values";
	}
	foreach my $new (@headers) {
		if ($new !~ /HomeoDB/) {
			#print $new, "\n";
			if ($new =~ /'?([A-Z0-9]+)_([A-Z]+)_(.*)_($HG)_([0-9]+)'?$/i) {
				my $spec = $1;
				my $group = $2;
				my $classif = $3;
				my $hg = $4;
				my $num = $5;
				my $uid = ".$hg.$num";
				#print $classif, "\n";
				$classif =~ s/_/\./g;
				my $homeo = "";
				if (($homeo) = grep/$classif/i, @fams) {
					chomp $homeo;
					$new = "$spec\_$group\_$homeo\_$hg\_$num";
				}
				my $origwseq = `grep -A1 -P -e '\|$uid' FastaDomains/$HG.fasta`;
				my ($origheader, $sequence) = split /\n/, $origwseq;
				chomp $sequence;
				chomp $new;
				$string .= ">$new\n$sequence\n";
			}
		}
	}
	return $string;
}
exit;
