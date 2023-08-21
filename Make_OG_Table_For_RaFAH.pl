#! /opt/ohpc/pub/apps/perl/5.28/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Bio::SearchIO;

my $extension = ".faaxDB"; #extension of files to be processed
my $file_prefix = "RaFAH"; #String to be added as prefix of the output files
my $suffix = "";
my $valid_ogs_file = "/mnt/lustre/bio/users/fcoutinho/RaFAH/HP_Ranger_Model_3_Valid_Cols.txt";
my $max_evalue = 0.00001; #Maximum evalue to consider a match between a hmmer model and a protein sequence when parsing the hmmsearch file
my $min_score = 50; 

GetOptions(
'extension=s' => \$extension,
);

my %seq_info;
my %valids_list;
my %hmmer_scores;

calc_scores();


sub calc_scores {
	my $hits_file = $_[0];
	my %seen_queries;
	my %seen_scaffolds;
	

	#Read the list of valid ogs (i.e. those that MUST be in the Genome x OG table) 
	open INPUT, "< $valid_ogs_file" or die "$!";
	while (my $line = <INPUT>) {
		chomp $line;
		my @values = split /\t/, $line;
		$valids_list{$values[0]} = 1;
	}
	close INPUT;
	my $valids_count = keys %valids_list;
	print "Obtained $valids_count ids from $valid_ogs_file\n";

	
	my @hits_files = glob("*$extension");
	foreach my $hits_file (@hits_files) {
		print "Parsing $hits_file\n";
		#Parse the hmmsearch output 
		my $in = Bio::SearchIO->new( -file => $hits_file, -format => 'hmmer' );
		while ( my $result = $in -> next_result ) {
			 while( my $hit = $result->next_hit ) {
				while( my $hsp = $hit->next_hsp ) {            
					#my $hit_desc = $hit->description;
					my $evalue = $hsp->evalue;
					my $score = $hsp->score;
					next if ($evalue > $max_evalue);
					next if ($score < $min_score);

					my $func = $hit->name;#." (".$hit_desc.")";
					my $query = $result->query_name();
					#my $query_desc = $result->query_description();
					

					next unless (defined $valids_list{$query});

					
					my $scaffold = $func;
					$scaffold =~ s/_(\d)+$//;
					#$scaffold =~ s/-cds(\d)+$//;
					
					
					#$seen_queries{$query} = 1;
					#$seen_scaffolds{$scaffold} = 1;
					
					$hmmer_scores{$query}{$scaffold} = $score unless (defined $hmmer_scores{$query}{$scaffold});
					
					if ($score > $hmmer_scores{$query}{$scaffold}) {
						$hmmer_scores{$query}{$scaffold} = $score;
					}
					
					}
				}
			}
		}
	
	my $genomexog_table_file_name = $file_prefix."_"."Genome_to_OG_Score_Min_Score_$min_score-Max_evalue_$max_evalue".$suffix.".tsv";
	my %params = ("matrix" => \%hmmer_scores, "file" => $genomexog_table_file_name, "missing" => "0");
	print_matrix_valids(\%params);
	
}


sub list_level2 {
    #Receives a reference to a 2D hash. returns a reference to a 1D hash in which the keys are all the keys from the 2nd dimesion of the input hash 
    my %hash = %{$_[0]};
    my %list;   
    foreach my $param1 (keys %hash) {
		my %hash2 = %{$hash{$param1}};
        foreach my $param2 (keys %hash2) {
            $list{$param2} = 1;
        }
    }
    return \%list;
}


sub print_matrix_valids {
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
	my @sorted_cols = sort keys %valids_list;
	
    #print the Header   
    foreach my $col (@sorted_cols) {
        print OUTPUT "$sep"."$col";
        }
    print OUTPUT "\n";
    #print the row names followed by the values and a line break
    foreach my $row (sort keys %rows) {
		#next unless (defined $seq_info{"Length"}{$row});
        print OUTPUT "$row";
        foreach my $col (@sorted_cols) {
            if (defined $matrix{$col}{$row}) {
				print OUTPUT "$sep"."$matrix{$col}{$row}";
			}  else {
				print OUTPUT "$sep"."0";
			}
        }
    print OUTPUT "\n";
    }
    close OUTPUT;
}
