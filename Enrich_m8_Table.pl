use warnings;
use strict;
use Getopt::Long;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;

my $input;
my $output = "Tax_Info.m8";
my $acc2taxid_dir;
my $acc2taxid_ext = "uniref100Id2ncbiTaxid.txt";
my $description_file;
my $names_file;
my $nodes_file;


GetOptions(
'input=s' => \$input,
'output=s' => \$output,
'nodes=s' => \$nodes_file,
'names=s' => \$names_file,
'description=s' => \$description_file,
'acc2taxid_dir=s' => \$acc2taxid_dir,
'acc2taxid_extension=s' => \$acc2taxid_ext,
) or die "Missing arguments!\n";

my %acc2taxid;

my @acc2taxidfiles = glob ("$acc2taxid_dir/*$acc2taxid_ext");

foreach my $acc2taxid_file (@acc2taxidfiles) {
	print "Parsing $acc2taxid_file\n";
	open INPUT, "< $acc2taxid_file" or die "$!";
	my $shouter = 1;
	my $counter = 0;
	while (my $line = <INPUT>) {
		$counter++;
		if ($shouter == $counter) {$shouter += $counter; print "Processed $counter entries\n";}
		chomp $line;
		my @values = split /\t/, $line;
		$acc2taxid{$values[0]} = $values[1];
		#last if ($counter >= 10_000_000);
		}
	close INPUT
}


print "Loading Taxonomy database\n";
my $taxDB = Bio::LITE::Taxonomy::NCBI->new (db=>"NCBI",names=> $names_file,nodes=> $nodes_file);

my %id2desc;
if ($description_file) {
	print "Loading sequence descriptions\n";
	open INPUT, "< $description_file" or die "$!";
	my $shouter = 1;
	my $counter = 0;
	while (my $line = <INPUT>) {
		$counter++;
		if ($shouter == $counter) {$shouter += $counter; print "Processed $counter entries\n";}
		chomp $line;
		my @values = split /\s/, $line;
		my $id = $values[0];
		$id2desc{$id} = $line;
		#last if ($counter >= 10_000_000);
		}
	close INPUT
	
}

my $shouter = 1;
my $counter = 0;
my $missed_counter=0;


print "Parsing $input\n";
open INPUT, "< $input" or die "$!";
open OUTPUT, "> $output" or die "$!";

while (my $line = <INPUT>) {
	$counter++;
	if ($shouter == $counter) {$shouter += $counter; print "Processed $counter lines\n";}
	chomp $line;
	my @values = split /\t/, $line;
	my $acc = $values[1];
	my $tax_id = "NA";
	my ($desc,$skingdom,$phylum,$class,$order,$family,$genus,$species) = qw(NA NA NA NA NA NA NA NA);
	if (defined $acc2taxid{$acc}) {
		$tax_id = $acc2taxid{$acc};
		$skingdom = $taxDB->get_term_at_level($tax_id,"superkingdom");
		$skingdom = "NA" if (($skingdom eq "") or ($skingdom eq "undef"));
		$phylum = $taxDB->get_term_at_level($tax_id,"phylum");
		$phylum = "NA" if (($phylum eq "") or ($phylum eq "undef"));
		$class = $taxDB->get_term_at_level($tax_id,"class");
		$class = "NA" if (($class eq "") or ($class eq "undef"));
		$order = $taxDB->get_term_at_level($tax_id,"order");
		$order = "NA" if (($order eq "") or ($order eq "undef"));
		$family = $taxDB->get_term_at_level($tax_id,"family");
		$family = "NA" if (($family eq "") or ($family eq "undef"));
		$genus = $taxDB->get_term_at_level($tax_id,"genus");
		$genus = "NA" if (($genus eq "") or ($genus eq "undef"));
		$species = $taxDB->get_term_at_level($tax_id,"species");
		$species = "NA" if (($species eq "") or ($species eq "undef"));
		} else {
		$missed_counter++;
	}
	$desc = $id2desc{$acc} if (defined $id2desc{$acc});
	my $fullinfo = join("\t",$desc,$tax_id,$skingdom,$phylum,$class,$order,$family,$genus,$species);
	print OUTPUT "$line\t$fullinfo\n";
	}
	
close INPUT;
close OUTPUT;

my $perc = ($missed_counter/$counter)*100;
$perc = int($perc);
print "No taxid found for $missed_counter entries out of $counter ($perc %)\n";