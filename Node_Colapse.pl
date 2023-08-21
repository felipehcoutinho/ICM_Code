use warnings;
use strict;
use Getopt::Long;
use Bio::TreeIO;

my $tree_file;
my $format = "newick";
my $seq_info_file;
my $category;
my $collapse_full;

GetOptions(
'format=s' => \$format,
'tree_file=s' => \$tree_file,
'seq_info_file=s' => \$seq_info_file,
'category=s' => \$category,
'collapse_full' => \$collapse_full,
) or die "Missing arguments!\n";

my %params = ("file" => $seq_info_file);
my $info_ref = read_table(\%params);
my %info = %$info_ref;

print "Opening tree\n";
# parse in tree
my $tree_parser = new Bio::TreeIO(-file   => $tree_file, -format => $format);
my $tree_obj = $tree_parser->next_tree;
my @nodes = $tree_obj->get_nodes;

open NODECOLLAPSE, "> Nodes_To_Colapse.txt" or die "$!";
print NODECOLLAPSE "COLLAPSE\nDATA\n";

my $node_count = -1;
my $shouter = 1;
my $root_node;
my %skip_nodes;

print "Parsing tree\n";

if ($collapse_full) {
	print "Looking for nodes that define $category\n";
	open NODENIDS, "> Node_New_IDS_$category.txt" or die "$!";
	print NODENIDS "DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tNode_New_IDs\nCOLOR\t#000000\nDATA\n";
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
			print "Collapsed single leaf $parent_id assigned to $category $info{$category}{$parent_id}\n";
			print NODECOLLAPSE "$parent_id\n";
			print NODENIDS "$parent_id\t$info{$category}{$parent_id}\t-1\t#000000\tbold\t1\t0\n";
			next;
		}
		my @children = $node->get_all_Descendents;
		my $children_count = 0;
		my @leaves;
		my $class_count = 0;
		my %seen_vals;

		foreach my $child (@children) {
			next unless ($child->is_Leaf);
			$children_count++;
			my $leaf_id = $child->id;
			if ((defined $info{$category}{$leaf_id}) and ($info{$category}{$leaf_id} ne "Unknown")) {
				my $categ_val = $info{$category}{$leaf_id};
				$seen_vals{$categ_val} = 1;
				$class_count++;
			} else {
				print "Missing $category info for $leaf_id\n";
			}
			
		}
		
		my @vals_list = keys %seen_vals;
		my $seen_val_count = @vals_list;
		if ($seen_val_count == 1) {
			print "Collapsed $parent_id with $children_count children and $seen_val_count valid category values: $vals_list[0]\n";
			print NODECOLLAPSE "$parent_id\n";
			print NODENIDS "$parent_id\t$vals_list[0]\t-1\t#000000\tnormal\t1\t0\n";
			foreach my $child (@children) {
				$skip_nodes{$child->id} = 1;
			}
		}		
	}
} else {
	print "Looking for nodes in which all leaves have no defined $category\n";
	foreach my $node (@nodes) {
		$node_count++;
		#Keep track of the root node
		unless ($root_node) {$root_node = $node};
		if ($node_count == $shouter) { print "Processed $node_count nodes\n"; $shouter+= $node_count}
		my $parent_id = $node->id;
		next if (defined $skip_nodes{$parent_id});
		#Check if node is a leaf
		my $leaf_status = $node->is_Leaf;
		if ($leaf_status) { next };
		my @children = $node->get_all_Descendents;
		my $children_count = 0;
		my @leaves;
		my $class_count = 0;

		foreach my $child (@children) {
			next unless ($child->is_Leaf);
			$children_count++;
			my $leaf_id = $child->id;
			$class_count++ if ((defined $info{$category}{$leaf_id}) and ($info{$category}{$leaf_id} ne "Unknown"));
			last if ($class_count > 0);
		}
		
		if ($class_count == 0) {
			print "Collapsed $parent_id with $children_count children and $class_count valid category values\n";
			print NODECOLLAPSE "$parent_id\n";
			foreach my $child (@children) {
				my $node_id = $child->id;
				print "\tAdding $node_id to skip nodes\n";
				$skip_nodes{$node_id} = 1;
			}

		} else {
			print "$parent_id will be kept open\n"
		}
	}	
	
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
	chomp $header;
	my @cols = split /$sep/, $header;
	
	my $values_count;
	while (my $line = <INPUT>) {
		#Read each line and assign it to the 2D hash
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
