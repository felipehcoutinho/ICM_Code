use warnings;
use strict;
use Getopt::Long;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Bio::SeqIO;
use Bio::TreeIO;

my $input_fasta;
my $output = "Tax_Info.m8";
my $acc2taxid_dir = "/mnt/lustre/repos/bio/databases/public/UniProtKB/UniProt_Release_2021-04-07/";
my $acc2taxid_ext = "uniref100Id2ncbiTaxid.txt";
my $description_file = "/mnt/lustre/bio/users/fcoutinho/Databases/UniRef100/UniRef100_Release_2021_04_07_Id_to_Desc.txt";
my $names_file = "/mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/names.dmp";
my $nodes_file = "/mnt/lustre/repos/bio/databases/public/NCBI/taxonomy/taxonomy_2021-04-06/nodes.dmp";
my $get_taxonomy;
my $labels_category;
my $binary_category;
my $binary_shape = 3;
my $color_category;
my $color_palette = "basic";
my $min_trait_count = 0;
my $input_info_file;
my $output_info_file = "Seq_Info.tsv";
my $protein_mode;
my $barplot_category;
my $tree_file;
my $format = "newick";

GetOptions(
'format=s' => \$format,
'tree_file=s' => \$tree_file,
'input_fasta=s' => \$input_fasta,
'output=s' => \$output,
'nodes=s' => \$nodes_file,
'names=s' => \$names_file,
'description=s' => \$description_file,
'acc2taxid_dir=s' => \$acc2taxid_dir,
'acc2taxid_extension=s' => \$acc2taxid_ext,
'get_taxonomy' => \$get_taxonomy,
'color_category=s' => \$color_category,
'color_palette=s' => \$color_palette,
'min_trait_count=s' => \$min_trait_count,
'labels_category=s' => \$labels_category,
'binary_category=s' => \$binary_category,
'binary_shape=s' => \$binary_shape,
'barplot_category=s' => \$barplot_category,
'input_info_file=s' => \$input_info_file,
'output_info_file=s' => \$output_info_file,
'protein_mode' => \$protein_mode,
) or die "Missing arguments!\n";


my %seq_info;
central();

sub central {
	if ($input_info_file) {
		print "Reading in $input_info_file\n";
		my %params = ("file" => $input_info_file);
		my $info_ref = read_table(\%params);
		%seq_info = %$info_ref;
	}
	if ($input_fasta) {
		index_seqs();	
	}
	if (($input_fasta) and ($get_taxonomy)) {
		get_taxonomy();		
	}
	if ($color_category) {
		make_tree_colors();
	}
	if ($labels_category) {
		make_tree_labels();
	}
	if ($binary_category) {
		make_tree_symbols();
	}
	if ($barplot_category) {
		make_tree_barplots();
	}
	print_seq_info();
}

sub count_node_traits {
	my $category = $_[0];
	
	print "Opening tree\n";
	my $tree_parser = new Bio::TreeIO(-file  => $tree_file, -format => $format);
	my $tree_obj = $tree_parser->next_tree;
	my @nodes = $tree_obj->get_nodes;
	
	my %trait_node_counts;
	my %seen_traits;
	my $node_count = -1;
	my $shouter = 1;
	my $root_node;
	my %skip_nodes;
	print "Counting occurrences of $category across tree nodes\n";
	foreach my $node (@nodes) {
		$node_count++;
		#Keep track of the root node
		unless ($root_node) {$root_node = $node};
		if ($node_count == $shouter) { print "Processed $node_count nodes\n"; $shouter+= $node_count}
		my $parent_id = $node->id;
		next if (defined $skip_nodes{$parent_id});
		#Check if node is a leaf
		my $leaf_status = $node->is_Leaf;
		if ($leaf_status) {
			#next;
		}
		my @children = $node->get_all_Descendents;
		my $children_count = 0;
		my @leaves;

		foreach my $child (@children) {
			next unless ($child->is_Leaf);
			$children_count++;
			my $leaf_id = $child->id;
			if ((defined $seq_info{$category}{$leaf_id}) and ($seq_info{$category}{$leaf_id} ne "Unknown")) {
				my $categ_val = $seq_info{$category}{$leaf_id};
				$seen_traits{$categ_val} = 1;
				$trait_node_counts{$parent_id}{$categ_val}++;
			}
		}		
	}
	my $unique_trait_count = keys %seen_traits;
	print("Detected $unique_trait_count unique traits across $node_count nodes\n");
	return(\%trait_node_counts,\%seen_traits);
}

sub make_tree_barplots {		
	my ($ref_trait_node_counts,$ref_seen_vals) = count_node_traits($barplot_category);
	my %trait_node_counts = %$ref_trait_node_counts;
	my %seen_traits = %$ref_seen_vals;
	my @traits_list = sort keys %seen_traits;
	my $traits_count = @traits_list;
	my $traits_string = join("\t",@traits_list);
	my $colors_ref = get_palette($color_palette,$traits_count);
	my @colors = @$colors_ref;
	my $colors_string = join("\t",@colors);
	my $shapes_string = $traits_string;
	$shapes_string =~ s/(\w)+/1/g;
	print "Printing output barplot file\n";
	open OUT, "> Tree_Barplot_$barplot_category.txt" or die "$!";
	print OUT 
"DATASET_MULTIBAR
SEPARATOR TAB
DATASET_LABEL	Barplot_$barplot_category
COLOR	#000000
FIELD_COLORS	$colors_string
FIELD_LABELS	$traits_string
LEGEND_TITLE	$barplot_category
LEGEND_SHAPES	$shapes_string
LEGEND_COLORS	$colors_string
LEGEND_LABELS	$traits_string
DATASET_SCALE	100-100-#000000-3-1-5	500-500-#000000-3-1-5	1000-1000-#000000-3-1-5
WIDTH	250
BORDER_WIDTH	1
BORDER_COLOR	#000000
DATA
";
	foreach my $node_id (keys %trait_node_counts) {
		print OUT "$node_id";
		foreach my $trait (keys %seen_traits) {
			my $count = 0;
			if (defined $trait_node_counts{$node_id}{$trait}) { 
				$count = $trait_node_counts{$node_id}{$trait}
			}
			print OUT "\t$count";
		}
		print OUT "\n",
	}
close OUT;
}

sub make_tree_symbols {
	#Build the file to be used by iTOL
	print "Printing output binary file\n";
	open OUT, "> Tree_Binary_$binary_category.txt" or die "$!";
	print OUT "DATASET_BINARY
SEPARATOR TAB
DATASET_LABEL\t$binary_category
COLOR\t#ff0000
FIELD_SHAPES\t$binary_shape
FIELD_LABELS\t$binary_category
DATA
";
	foreach my $seq_id (keys %{$seq_info{$binary_category}}) {
		my $trait = $seq_info{$binary_category}{$seq_id};
		if ($trait ne "NA") {
			print OUT "$seq_id\t1\n";
		}
	}
	close OUT;	
}


sub print_seq_info {
	print "Printing seq info to $output_info_file\n";
	my %params = ("matrix" => \%seq_info, "file" => $output_info_file);
	print_matrix(\%params);	
}

sub make_tree_labels {
	#Build the file to be used by iTOL
	print "Printing output labels file\n";
	open OUT, "> Tree_Labels_$labels_category.txt" or die "$!";
	print OUT "LABELS
SEPARATOR TAB
DATA
";
	foreach my $seq_id (keys %{$seq_info{$labels_category}}) {
		my $trait = $seq_info{$labels_category}{$seq_id};
		print OUT "$seq_id\t$trait\n";
	}
	close OUT;	
}

sub get_palette {
	print "Using $_[0] palette with $_[1] colours\n";
	my $pal_name = $_[0];
	my $ncol = $_[1] - 1;
	my @colors;
	@colors = qw(#d11141 #00b159 #00aedb #f37735 #ffc425) if ($pal_name =~ /basic/i);
	@colors = qw(#cccccc #a2a39f) if ($pal_name =~ /gray/i);
	@colors = qw(#ffb3ba #ffdfba #ffffba #baffc9 #bae1ff) if ($pal_name =~ /pastel/i);
	@colors = qw(#e6194b #3cb44b #ffe119 #0082c8 #f58231 #911eb4 #46f0f0 #f032e6 #d2f53c #fabebe #008080 #e6beff #aa6e28 #fffac8 #800000 #aaffc3 #808000 #ffd8b1 #000080 #808080) if ($pal_name =~ /sasha/i);
	@colors = ("#E62E2E","#E6CB2E","#62E62E","#2EE697","#2E97E6","#622EE6","#E62ECB") if ($pal_name =~ /rainbow7/i);
	@colors = ("#E62E2E","#E69C2E","#C1E62E","#53E62E","#2EE677","#2EE6E6","#2E77E6","#532EE6","#C12EE6","#E62E9C")  if ($pal_name =~ /rainbow10/i);
	@colors = ("#E62E2E","#E6772E","#E6C12E","#C1E62E","#77E62E","#2EE62E","#2EE677","#2EE6C1","#2EC1E6","#2E77E6","#2E2EE6","#772EE6","#C12EE6","#E62EC1","#E62E77") if ($pal_name =~ /rainbow15/i);
	@colors = ("#E62E2E","#E6652E","#E69C2E","#E6D32E","#C1E62E","#8AE62E","#53E62E","#2EE640","#2EE677","#2EE6AE","#2EE6E5","#2EAEE6","#2E77E6","#2E40E6","#532EE6","#8A2EE6","#C12EE6","#E62ED3","#E62E9C","#E62E65") if ($pal_name =~ /rainbow20/i);
	@colors = qw(#E62E2E #E63D2E #E64B2E #E65A2E #E6692E #E6772E #E6862E #E6952E #E6A32E #E6B22E #E6C12E #E6CF2E #E6DE2E #DEE62E #CFE62E #C1E62E #B2E62E #A3E62E #95E62E #86E62E #77E62E #69E62E #5AE62E #4BE62E #3DE62E #2EE62E #2EE63D #2EE64B #2EE65A #2EE669 #2EE677 #2EE686 #2EE695 #2EE6A3 #2EE6B2 #2EE6C1 #2EE6CF #2EE6DE #2EDEE6 #2ECFE6 #2EC1E6 #2EB2E6 #2EA3E6 #2E95E6 #2E86E6 #2E77E6 #2E69E6 #2E5AE6 #2E4BE6 #2E3DE6 #2E2EE6 #3D2EE6 #4B2EE6 #5A2EE6 #692EE6 #772EE6 #862EE6 #952EE6 #A32EE6 #B22EE6 #C12EE6 #CF2EE6 #DE2EE6 #E62EDE #E62ECF #E62EC1 #E62EB2 #E62EA3 #E62E95 #E62E86 #E62E77 #E62E69 #E62E5A #E62E4B #E62E3DCC) if ($pal_name =~ /rainbow75/i);
	@colors = qw(#52ccd4 #d85331 #75d64d #5672c8 #c8ca4b #98522f #80d796 #d89d50 #3f8342 #797f35)  if ($pal_name =~ /custom1/i);
	@colors = qw(#081D58 #7FCDBB #FFFFD9)  if ($pal_name =~ /zone/i);
	@colors = @colors[0..$ncol];
	
	return(\@colors);
}

sub make_tree_colors {
	#Count the occurrences of each sequence trait. This is necessary to later sort the traits by their number of occurrences
	print "Counting trait occurrences\n";
	my @seq_ids;
	if ($input_fasta) {
		@seq_ids = keys %{$seq_info{"Original_File"}};
	} else {
		@seq_ids = keys %{$seq_info{$color_category}};
	}
	
	my %trait_counts;
	foreach my $seq_id (@seq_ids) {
		my $oid = $seq_id;
		if ($input_fasta) {
			next unless (defined $seq_info{"Original_File"}{$seq_id})
		}
		if ($protein_mode) {
			$oid =~ s/_(\d)+$//;
		}
		#print "SEQ $seq_id $oid TRAIT $seq_info{$color_category}{$oid}\n";
		my $trait = $seq_info{$color_category}{$oid};
		$trait_counts{$trait}++
	}
	
	#Define color palette based on user input
	my $unique_trait_count = keys %trait_counts;
	my $colors_ref = get_palette($color_palette,$unique_trait_count);
	my @colors = @$colors_ref;
	
	#sort the traits by their number of occurrences
	print "Sorting trait occurrences\n";
	my @sorted_traits =  sort { $trait_counts{$a} <=> $trait_counts{$b} } keys(%trait_counts);

	#Assign colors to each trait. Repeat colors if there are more traits than colors in the palette
	my $col_posit = 0;
	my $posit = 0;
	my $col_num = @colors;
	my $max_posit = $col_num--;
	my %trait2col;
	
	#Ignore "undef", "Unknown" and "NA" traits 
	foreach my $trait (@sorted_traits) {
		if (($trait eq "undef") or ($trait eq "Unknown") or ($trait eq "NA") or ($trait eq "")) {next}
		if ($trait_counts{$trait} < $min_trait_count) {next}
		print "$posit\t$trait\t$trait_counts{$trait}\t$colors[$posit]\n";
		$posit++;
		$posit = 0 if ($posit >= $max_posit);
		$trait2col{$trait} = $colors[$posit] 
	}
	
	#Build the file to be used by iTOL
	print "Printing output colors file\n";
	open OUT, "> Tree_Colors_By_Index_Number_$color_category-Palette_$color_palette.txt" or die "$!";
	print OUT "TREE_COLORS
SEPARATOR TAB
DATA
";

	foreach my $seq_id (@seq_ids) {
		my $oid = $seq_id;
		if ($protein_mode) {
			$oid =~ s/_(\d)+$//;
		}
		my $trait = $seq_info{$color_category}{$oid};
		if ($trait2col{$trait}) {
			$seq_info{"Color-$color_palette-$color_category"}{$seq_id} = $trait2col{$trait};
			#print "$seq_id\trange\t$trait2col{$trait}\t$trait\n";
			print OUT "$seq_id\trange\t$trait2col{$trait}\t$trait\n";
		}
	}
	close OUT;
}

sub index_seqs {
	#Iterate over the sequences in the input fasta file. Collect basic info 
	my @seq_files = ($input_fasta);
	my $counter = 0;
	foreach my $file (@seq_files) {
		my $file_name = $file;
		$file_name =~ s/(.)*\///;
		print "\tProcessing $file_name\n";
		my $seq_in = Bio::SeqIO->new(-file => "< $file", -extension => "fasta");
		while (my $seq_obj = $seq_in->next_seq) {
			$counter++;
			my $seq_id = $seq_obj->id;
			my $seq_length = $seq_obj->length;
			my $seq_desc = $seq_obj->desc;
			die "Repeated id $seq_id in $file\n" if (defined $seq_info{"Original_File"}{$seq_id});
			$seq_info{"Original_File"}{$seq_id} = $file_name;
			$seq_info{"Description"}{$seq_id} = $seq_desc;
			$seq_info{"Length"}{$seq_id} = $seq_length;
		} 
	}	
	print "\tProcessed $counter Sequences\n";
}

sub get_taxonomy {
	#Read in the acc2taxid files
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

	#Load the NCBI taxonomy 
	print "Loading Taxonomy database\n";
	my $taxDB = Bio::LITE::Taxonomy::NCBI->new (db=>"NCBI",names=> $names_file,nodes=> $nodes_file);
	my $shouter = 1;
	my $counter = 0;
	my $missed_counter=0;
	#Iterate over the sequences in the input fasta file.
	foreach my $seq_id (keys %{$seq_info{"Description"}}) {
		$counter++;
		if ($shouter == $counter) {$shouter += $counter; print "Processed Taxonomy for $counter sequences\n";}
		my $tax_id = "NA";
		#Assign taxonomy to these sequences for which there is taxid associated with their ids 
		my ($skingdom,$phylum,$class,$order,$family,$genus,$species) = qw(NA NA NA NA NA NA NA NA);
		$seq_info{"Domain"}{$seq_id} = $skingdom;
		$seq_info{"Phylum"}{$seq_id} = $phylum;
		$seq_info{"Class"}{$seq_id} = $class;
		$seq_info{"Order"}{$seq_id} = $order;
		$seq_info{"Family"}{$seq_id} = $family;
		$seq_info{"Genus"}{$seq_id} = $genus;
		$seq_info{"Species"}{$seq_id} = $species;
		
		if (defined $acc2taxid{$seq_id}) {
			$tax_id = $acc2taxid{$seq_id};
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
			
			$seq_info{"Domain"}{$seq_id} = $skingdom;
			$seq_info{"Phylum"}{$seq_id} = $phylum;
			$seq_info{"Class"}{$seq_id} = $class;
			$seq_info{"Order"}{$seq_id} = $order;
			$seq_info{"Family"}{$seq_id} = $family;
			$seq_info{"Genus"}{$seq_id} = $genus;
			$seq_info{"Species"}{$seq_id} = $species;
			
			#Use the class for proteobacteria instead of phylym
			if ($phylum eq "Proteobacteria") {
				$seq_info{"Phylum"}{$seq_id} = $class;	
			}
		} else {
			$missed_counter++;
		}
	}
	my $perc = ($missed_counter/$counter)*100;
	$perc = int($perc);
	print "No taxid found for $missed_counter entries out of $counter ($perc %)\n";
}



sub read_table {
#Reads a table as a 2D hash. Columns are the first dimension and rows the second. Receives a reference of a hash as the parameters to run the subroutine  
	my %parameters	= %{$_[0]};
	die "No Input File Defined.\n" unless (defined $parameters{"file"});
	my $file = $parameters{"file"};
	my $sep = "\t";
	$sep = $parameters{"separator"} if (defined $parameters{"separator"});	
	
	my %matrix;

	open INPUT, "< $file" or die "$!";
	my $header = <INPUT>;
	$header =~ s/\r$//;
	chomp $header;
	my @cols = split /$sep/, $header;
	
	my $values_count;
	while (my $line = <INPUT>) {
		#Read each line and assign it to the 2D hash
		$line =~ s/\r$//;
		chomp $line;
		my @rows = split /$sep/, $line;
		foreach my $posit (1..(@cols - 1)) {
			$matrix{$cols[$posit]}{$rows[0]} = $rows[$posit];
			$values_count++;
		}
	}
	close INPUT;
	print "Obtained $values_count values from $file.\n";
	return \%matrix;
}


sub print_matrix {
    my %parameters  = %{$_[0]};
         
    die "No Matrix (2D Hash) Defined.\n" unless (defined $parameters{"matrix"});
    my %matrix = %{$parameters{"matrix"}};  
    my $output_file = "Output_Matrix";
    my $sep = "\t"; 
    my $missing = "NA";
	
	$missing = $parameters{"missing"} if (defined $parameters{"missing"});
    $output_file = $parameters{"file"} if (defined $parameters{"file"});
    $sep = $parameters{"separator"} if (defined $parameters{"separator"});      
	
     
    my $rows_ref = list_level2(\%matrix);
    my %rows = %{$rows_ref};
 
    open OUTPUT, "> $output_file" or die "$!";   
    print OUTPUT "Variable";
    #print the Header   
    foreach my $col (sort keys %matrix) {
        print OUTPUT "$sep"."$col";
        }
    print OUTPUT "\n";
    #print the row names followed by the values and a line break
    foreach my $row (sort keys %rows) {
        print OUTPUT "$row";
        foreach my $col (sort keys %matrix) {
            if (defined $matrix{$col}{$row}) {print OUTPUT "$sep"."$matrix{$col}{$row}"} else {print OUTPUT "$sep"."$missing"};
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
		my %hash2 = %{$hash{$param1}};
        foreach my $param2 (keys %hash2) {
            $list{$param2} = 1;
        }
    }
    return \%list;
}
