use warnings;
use strict;
use Getopt::Long;


my $protein_mode;
my $clean_name_mode;
my $creep_names;
my $info_file;
my $matrix_file;
my $fields;
my $posit = 1;
my $ignore_version;

GetOptions(
'info_file=s' => \$info_file,
'matrix_file=s' => \$matrix_file,
'fields=s' => \$fields,
'posit=s' => \$posit,
'protein_mode' => \$protein_mode,
'creep_names' => \$creep_names,
'clean_name_mode' => \$clean_name_mode,
'ignore_version' => \$ignore_version,
) or die "Missing arguments!\n";

my %params = ("file" => $info_file);
my $info_ref = read_table(\%params);
my %info = %$info_ref;

open IN, "< $matrix_file" or die "$!";
$matrix_file =~ s/^(.)*\///;

open OUT, "> $matrix_file+$fields.tsv" or die "$!";


my $header = <IN>;
$fields =~ s/\s//g;
my @fields = split /,/, $fields;
foreach my $field (@fields) {
	print OUT "$field\t";
}
print OUT "$header";

$posit--;
my $found_trait_count = 0;
while (my $line = <IN>) {
	chomp $line;
	my @values = split /\t/, $line;
	if ($clean_name_mode) { $values[$posit] =~ s/\W+/_/g; $values[$posit] =~ s/(_)+/_/g;}
	if ($protein_mode) { 
		$values[$posit] =~ s/\|(.)*$// if ($creep_names);
		$values[$posit] =~ s/_(\d)+$//; 
		$values[$posit] =~ s/-cds(\d)+$//;
		}
	if ($ignore_version) {  $values[$posit] =~ s/\.(\d)$//; }
	foreach my $field (@fields) {
		my $trait = "Unknown";
		if (defined $info{$field}{$values[$posit]}) {$trait = $info{$field}{$values[$posit]}; $found_trait_count++; }
		print OUT "$trait\t";
	}
	print OUT "$line\n";
}
close IN;
print "Found $found_trait_count traits and added to output\n";


sub read_table {
#Reads a table as a 2D hash. Columns are the first dimension and rows the second. Receives a reference of a hash as the parameters to run the subroutine  
	my %parameters	= %{$_[0]};
	die "No Input File Defined.\n" unless (defined $parameters{"file"});
	my $file = $parameters{"file"};
	my $sep = "\t";
	my %blistR;
	my %blistC; 
	$sep = $parameters{"separator"} if (defined $parameters{"separator"});	
	%blistR = %{$parameters{"blistR"}} if (defined $parameters{"blistR"});
	%blistC = %{$parameters{"blistC"}} if (defined $parameters{"blistC"});
	
	my %matrix;

	open INPUT, "< $file" or die "$!";
	my $header = <INPUT>;
	chomp $header;
	my @cols = split /$sep/, $header;
	
	my $values_count;
	while (my $line = <INPUT>) {
		#Read each line and assign it to the 2D hash
		chomp $line;
		$line =~ s/\r//g;
		my @rows = split /$sep/, $line;
		if ($clean_name_mode) { $rows[0] =~ s/\W+/_/g; $rows[0] =~ s/(_)+/_/g;}
		next if (defined $blistR{$rows[0]});
		foreach my $posit (1..(@cols - 1)) {
			next if (defined $blistC{$cols[$posit]}); 			
			$matrix{$cols[$posit]}{$rows[0]} = $rows[$posit];
			$values_count++;
		}
	}
	close INPUT;
	print "Obtained $values_count values from $file.\n";
	return \%matrix;
}
