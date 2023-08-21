use warnings;
use strict;
use Getopt::Long;
use Bio::TreeIO;


my $in_file;
my $format = "newick";
my $out_file = "Tree_With_node_IDs.newick";

GetOptions(
'format=s' => \$format,
'input_tree=s' => \$in_file,
'output_tree=s' => \$out_file,
) or die "Missing arguments!\n";


print "Opening tree\n";
# parse in newick/new hampshire format
my $tree_parser = new Bio::TreeIO(-file   => $in_file, -format => $format);
my $tree_obj = $tree_parser->next_tree;
my @nodes = $tree_obj->get_nodes;

my $out = new Bio::TreeIO(-file => "> $out_file",
                          -format => 'newick');
						  
my $node_count = 0;
my $shouter = 1;

print "Parsing tree\n";
foreach my $node (@nodes) {
	if ($node_count == $shouter) { print "Processed $node_count nodes\n"; $shouter+= $node_count}
	my $leaf_status = $node->is_Leaf;
	next if ($leaf_status);

	$node->id("Node_$node_count");

	$node_count++;
}

 $out->write_tree($tree_obj);