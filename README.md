# Homeodomains

## Pipeline

Begin by collecting the sequences of interest, in this case all animal homeobox proteins:

`./fastaCollection.pl --terms "Homeobox,Homeodomain,ANTP,PRD,LIM,POU,HNF,SINE,TALE,CUT,PROS,ZF,CERS" --fastafile "<file>.csv" --occfile "<file>.fasta"`

Next the fastas should be sorted and sepearted into individual HG fasta files:

`./sortFasta.pl <file>.fasta`

Next for the classification step:

Download database of homeodb homeoboxes: http://homeodb.zoo.ox.ac.uk/download.get - Selecting the appropriate parameters.

`cat Fastas/*.fasta > All.fasta`

`makeblastdb -dbtype 'prot' -in homeodb.fasta`

`blastp -db 'homeodb' -query 'All.fasta' -evalue 10e-6 -out All.blast`

And run InterProScan:

`./interproscan.sh -appl 'Panther,Pfam' -i </path/to/All.fasta> -b </path/to/results.interpro>`

Run the domain extraction program:

`./domainPuller` - This one is hardcoded so may need adjustments to fit user directories etc.

Finally to actually classify:

`./classifyWrapper.pl` or to individually classify each $hg.fasta file:
`./homeoClassifier.pl --fastafile path/to/$hg.fasta --blastfile All.blast --interprofile new-interpro.tsv`
`./classifyChecker.pl` Quick checks for formatting as well as classifications

Then run the trees:

`treebuilder.pl <FastaDir>`
`SortFastaMore.pl <Intermediate File Directory with alignments>` To remove duplicate domains after trimming ahead of iqtree.

Final classification step:

`./homeoTreeParser.pl <IntermediateFiles> <classificationTable.tsv>` With directory containing all inferred trees in Newick format, and classificationTable.tsv as current classification log as verbose output from classifyWrapper.pl.

## Further graphical analyses

`./allhbxcounttable.pl TreeFiles > hbxCount.csv`
`./homeocountconvert.pl hbxCount.csv > hbxCountMelt.csv`
`./tempHbx.pl hbxCountMelt.csv`
`./tempHbx.pl hbxCountMelt.csv > hbxOrigins.csv`
`./lossGainhbx.pl hbxCountMelt.csv > hbxLossGain.csv`



