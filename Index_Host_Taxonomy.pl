use warnings;
use strict;
use Getopt::Long;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;

my $seq_info;
my $host_taxon_var = "Predicted_Host";
my $sep = "\t";
my $output_file = "Host_Info.tsv";
my $names_file = "/mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/names.dmp";
my $nodes_file = "/mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/nodes.dmp";


GetOptions(
'seq_info=s' => \$seq_info,
'sep=s' => \$sep,
'output_file=s' => \$output_file,
'host_taxon_var=s' => \$host_taxon_var,
) or die "Missing arguments!\n";


my $taxDB;
my %host_info;
my %taxonomy;
my %alias;

print "Reading $seq_info\n";
my %params = ("file" => $seq_info, "separator" => $sep);
my $info_ref = read_table(\%params);
my %seq_info = %$info_ref;

load_db();
parse();
%params = ("matrix" => \%seq_info, "file" => $output_file);
print_matrix(\%params);


sub load_db {
	print "Loading database\n";
	$taxDB = Bio::LITE::Taxonomy::NCBI->new (db=>"NCBI",names=> $names_file,nodes=> $nodes_file);;
}

sub parse {
	open OUT, "> Index_Host_Taxonomy.log" or die "$!";
	my $level2_ref = $seq_info{$host_taxon_var};
	my %level2 = %$level2_ref;
	foreach my $scaffold (keys %level2) {
		my $ident;
		my $host_name;
		
		my $host = $seq_info{$host_taxon_var}{$scaffold};
		my $host_tid = $taxDB->get_taxid_from_name($host);

		
		print OUT "Processing $scaffold\n";
		print OUT "host\t$host\n";
		
		if ($host_tid) {
			$ident = $host_tid;
			$host_name = $host;
		} else {
			print OUT "No Taxid found. Sequence was skipped\n";
			next;
		}
		
		print OUT "Host_Name\t$host_name\n";
			#print "Searching Taxonomy for $host\n";
			$ident = "1386" if ($host =~ /Bacillus/);
			$ident = "635" if ($host =~  /Edwardsiella/);
			$ident = "2053" if ($host =~  /Gordonia/);
			$ident = "581" if ($host =~  /Morganella/);
			$ident = "583" if ($host =~  /Proteus/);
			$ident = "265" if ($host =~  /Paracoccus/);
			$ident = "629" if ($host =~  /Yersinia/);
			$ident = "1827" if ($host =~  /Rhodococcus/);
			$ident = "82202" if ($host =~  /Centipeda/);
			$ident = "32207" if ($host =~  /Rothia/);
			$ident = "373984" if ($host =~  /Rivularia/);
			$ident = "43773" if ($host =~  /Syntrophus/);
			$ident = "776" if ($host =~  /Coxiella/);
			$ident = "159191" if ($host =~  /Nodularia/);

			print OUT "TaxID\t$ident\n";
		
			my $skingdom = $taxDB->get_term_at_level($ident,"superkingdom"); 
			my $phylum = $taxDB->get_term_at_level($ident,"phylum");
			my $class = $taxDB->get_term_at_level($ident,"class");
			my $order = $taxDB->get_term_at_level($ident,"order");
			my $family = $taxDB->get_term_at_level($ident,"family");
			my $genus = $taxDB->get_term_at_level($ident,"genus");
			my $species = $taxDB->get_term_at_level($ident,"species");
			
			print OUT "$skingdom,$phylum,$class,$order,$family,$species\n";
			
			$seq_info{$host_taxon_var."_1_Domain"}{$scaffold} = $skingdom if (defined $skingdom);
			$seq_info{$host_taxon_var."_2_Phylum"}{$scaffold} = $phylum if (defined $phylum);
			$seq_info{$host_taxon_var."_3_Class"}{$scaffold} = $class if (defined $class);
			$seq_info{$host_taxon_var."_4_Order"}{$scaffold} = $order if (defined $order);
			$seq_info{$host_taxon_var."_5_Family"}{$scaffold} = $family if (defined $family);
			$seq_info{$host_taxon_var."_6_Genus"}{$scaffold} = $genus if (defined $genus);
			$seq_info{$host_taxon_var."_7_Species"}{$scaffold} = $species if (defined $species);
		
		
		}	
	close OUT;
}




sub print_matrix {
    my %parameters  = %{$_[0]};
         
    die "No Matrix (2D Hash) Defined.\n" unless (defined $parameters{"matrix"});
    my %matrix = %{$parameters{"matrix"}};  
    my $output_file = "Output_Matrix";
    my $sep = "\t"; 
     
    $output_file = $parameters{"file"} if (defined $parameters{"file"});
    $sep = $parameters{"separator"} if (defined $parameters{"separator"});      
 
     
    my $rows_ref = list_level2(\%matrix);
    my %rows = %{$rows_ref};
 
    open OUTPUT, "> $output_file" or die "$!";   
    print OUTPUT "Sequence";
    #print the Header   
    foreach my $col (sort keys %matrix) {
        print OUTPUT "$sep"."$col";
        }
    print OUTPUT "\n";
    #print the row names followed by the values and a line break
    foreach my $row (sort keys %rows) {
		#next unless (defined $seq_info{"Length"}{$row});
        print OUTPUT "$row";
        foreach my $col (sort keys %matrix) {
            if ((defined $matrix{$col}{$row} and $matrix{$col}{$row} ne "undef")) {print OUTPUT "$sep"."$matrix{$col}{$row}"} else {print OUTPUT "$sep"."NA"};
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
		my $level2_ref = $hash{$param1};
		my %level2 = %$level2_ref;
        foreach my $param2 (keys %level2) {
            $list{$param2} = 1;
        }
    }
    return \%list;
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
	chomp $header;
	my @cols = split /$sep/, $header;
	
	my $values_count;
	while (my $line = <INPUT>) {
		#Read each line and assign it to the 2D hash
		chomp $line;
		my @rows = split /$sep/, $line;
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

