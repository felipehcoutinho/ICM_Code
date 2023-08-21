use warnings;
use strict;
use Getopt::Long;


my $info_file;
my $matrix_file;
my $true_virus;
my $filter_var;
my $min_filter_val = 0;
my $max_filter_val = 999_999_999_999;

GetOptions(
'info_file=s' => \$info_file,
'matrix_file=s' => \$matrix_file,
'min_filter_val=s' => \$min_filter_val,
'filter_var=s' => \$filter_var,
'max_filter_val=s' => \$max_filter_val,
'true_virus' => \$true_virus,
) or die "Missing arguments!\n";


print "Reading in $info_file\n";
my %params = ("file" => $info_file);
my $info_ref = read_table(\%params);
my %info = %$info_ref;

print "Reading in $matrix_file\n";
%params = ("file" => $matrix_file);
my $matrix_ref = read_table(\%params);
my %matrix = %$matrix_ref;
$matrix_file =~ s/^.*\///;

my @categories = keys %info;

foreach my $categ (@categories) {
	my %c_abund;
	if ($categ =~ /Hit_Counts|Perc_Matched_PEGs|CR_Matched_PEGs|CR_AAI|Score|Confidende|Membership/) {
		print "Skipping $categ\n";
		next;
	}
	
	unless ($categ =~ /Recruitment|Lineage|Level|Taxonomy|Host|Source|Baltimore|Domain|Phylum|Class|Order|Family|Genus|Species|Bin|CR_Taxon|Salinity|Sample|VP|Virus_Type|Bacphlip_Classification/i) {
		print "Skipping $categ\n";
		next;
	}
	print "Computing abundances for category $categ\n";
	foreach my $col (keys %matrix) {
		next if ($col =~ /(.*-var$)|contigLen|totalAvgDepth/);
		foreach my $row (keys %{$matrix{$col}}) {
			if (($true_virus) and ($info{"True_Virus"}{$row} ne "TRUE")) { next }
			my $nrow = "Unknown";
			$nrow = $info{$categ}{$row} if ((defined $info{$categ}{$row}) and ($info{$categ}{$row} ne "NA"));
			if (($categ eq "Phylum") and ($nrow eq "Proteobacteria") and ($info{"Class"}{$row} ne "NA") and ($info{"Class"}{$row} ne "")) { $nrow = $info{"Class"}{$row} }
			if (($categ eq "Level.2.CR_Taxon") and ($nrow eq "Proteobacteria") and ($info{"Level.3.CR_Taxon"}{$row} ne "NA") and ($info{"Level.3.CR_Taxon"}{$row} ne "")) { $nrow = $info{"Level.3.CR_Taxon"}{$row} }
			if (($categ eq "Host_Taxonomy_2_Phylum") and ($nrow eq "Proteobacteria") and ($info{"Host_Taxonomy_3_Class"}{$row} ne "NA") and ($info{"Host_Taxonomy_3_Class"}{$row} ne "")) { $nrow = $info{"Host_Taxonomy_3_Class"}{$row} }
			if (($categ eq "Predicted_Host_2_Phylum") and ($nrow eq "Proteobacteria") and ($info{"Predicted_Host_3_Class"}{$row} ne "NA") and ($info{"Predicted_Host_3_Class"}{$row} ne "")) { $nrow = $info{"Predicted_Host_3_Class"}{$row} }
			if (($categ eq "PHIST_Predicted_Host_2_Phylum") and ($nrow eq "Proteobacteria") and ($info{"PHIST_Predicted_Host_3_Class"}{$row} ne "NA") and ($info{"PHIST_Predicted_Host_3_Class"}{$row} ne "")) { $nrow = $info{"PHIST_Predicted_Host_3_Class"}{$row} }
			if (($categ eq "hostPhylum") and ($nrow eq "Proteobacteria") and ($info{"hostClass"}{$row} ne "NA") and ($info{"hostClass"}{$row} ne "")) { $nrow = $info{"hostClass"}{$row} }
			if (($categ eq "phylum") and ($nrow eq "Proteobacteria") and ($info{"class"}{$row} ne "NA") and ($info{"class"}{$row} ne "")) { $nrow = $info{"class"}{$row} }
			if ($filter_var) {
				unless (defined $info{$filter_var}{$row}) {
					next;
				}
				my $val = $info{$filter_var}{$row};
				if ($val eq "NA") { $val = 0}
				if (($val >= $min_filter_val) and ($val <= $max_filter_val)) {
					$c_abund{$col}{$nrow} += $matrix{$col}{$row};
				}
			} else {
				$c_abund{$col}{$nrow} += $matrix{$col}{$row};
			}
			
		}
	}
	$categ =~ s/\W/_/g;
	my $out_file = "$categ-$matrix_file";
	if ($true_virus) {$out_file = "True_Virus-".$out_file}
	if ($filter_var) {$out_file = "Filtered_by_$filter_var-".$out_file}
	($params{"file"},$params{"matrix"}) = ($out_file,\%c_abund);
	print "Printing abundances for category to $out_file\n";
	print_matrix(\%params);
}


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
	$header =~ s/\r//g;
	chomp $header;
	my @cols = split /$sep/, $header;
	
	my $values_count;
	while (my $line = <INPUT>) {
		#Read each line and assign it to the 2D hash
		$line =~ s/\r//g;
		chomp $line;
		my @rows = split /$sep/, $line;
		#next if (defined $blistR{$rows[0]});
		foreach my $posit (1..(@cols - 1)) {
			#next if (defined $blistC{$cols[$posit]}); 			
			$matrix{$cols[$posit]}{$rows[0]} = $rows[$posit];
			$values_count++;
		}
	}
	close INPUT;
	print "Obtained $values_count values from $file.\n";
	return \%matrix;
}


sub print_matrix {
	my %parameters	= %{$_[0]};
		
	die "No Matrix (2D Hash) Defined.\n" unless (defined $parameters{"matrix"});
	my %matrix = %{$parameters{"matrix"}};	
	my $output_file = "Output_Matrix";
	my $sep = "\t";	
	
	$output_file = $parameters{"file"} if (defined $parameters{"file"});
	$sep = $parameters{"separator"} if (defined $parameters{"separator"});		

	
	my $rows_ref = list_level2(\%matrix);
	my %rows = %{$rows_ref};

	open OUTPUT, "> $output_file" or die "$!";	
	print OUTPUT "Sample";
	#print the Header	
	foreach my $col (sort keys %matrix) {
		print OUTPUT "$sep"."$col";
		}
	print OUTPUT "\n";
	#print the row names followed by the values and a line break
	foreach my $row (sort keys %rows) {
		print OUTPUT "$row";
		foreach my $col (sort keys %matrix) {
			if (defined $matrix{$col}{$row}) {print OUTPUT "$sep"."$matrix{$col}{$row}"} else {print OUTPUT "$sep"."0"}
		}
	print OUTPUT "\n";
	}
	close OUTPUT;
}


sub list_level2 {
	#Receives a reference to a 2D hash. returns a reference to a 1D hash in which the keys are all the keys from the 2nd dimesion of the input hash	
	my %hash = %{$_[0]};
	my %list;	
	foreach my $param1 (keys %hash) {
		foreach my $param2 (keys %{$hash{$param1}}) {
			$list{$param2} = 1;
		}
	}
	return \%list;
}



