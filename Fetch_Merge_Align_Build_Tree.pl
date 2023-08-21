use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $list_file;
my $fetch_fasta;
my $merge_fasta;
my $skip_reps;
my $max_length = 999_999_999;
my $min_length = 0;

GetOptions(
'list_file=s' => \$list_file,
'fetch_fasta=s' => \$fetch_fasta,
'merge_fasta=s' => \$merge_fasta,
'min_length=s' => \$min_length,
'max_length=s' => \$max_length,
'skip_reps' => \$skip_reps,
) or die "Missing arguments!\n";

my %list;
open IN, "< $list_file" or die "$!";
while (my $line = <IN>) {
	$line =~ s/\r//g;
	chomp $line;
	$list{$line} = 1;
}
close IN;

my $list_count = keys %list;
print "Obtained $list_count unique IDs from $list_file\n";

my $out_name = $list_file;
$out_name =~ s/^(\.)*\///;
$out_name =~ s/\.(\w)*$//;
$out_name = $out_name."+".$merge_fasta;
#$out_name .= ".fasta";

my $seq_in = Bio::SeqIO->new('-file' => "< $fetch_fasta", '-format' => "Fasta");
my $seq_out = Bio::SeqIO->new('-file' => "> $out_name", '-format' => "Fasta");

my $counter = 0;
my %seen_ids;
my %seen_genomes;
my %passed_list;

while (my $seq_obj = $seq_in->next_seq) {
	my $id = $seq_obj->id();
	my $length = $seq_obj->length();
	if ((defined $seen_ids{$id}) and ($skip_reps)) {print "Skipping repeated ID: $id\n"; next}
	$seen_ids{$id} = 1;
	if ((defined $list{$id}) and ($length <= $max_length) and ($length >= $min_length)) {
		$counter++;
		$seq_out->write_seq($seq_obj);
		$passed_list{$id} = 1;
	}
}
print "$counter sequences added from $fetch_fasta.\n";

$counter = 0;
$seq_in = Bio::SeqIO->new('-file' => "< $merge_fasta", '-format' => "Fasta");
while (my $seq_obj = $seq_in->next_seq) {
	my $length = $seq_obj->length();
	if (($length <= $max_length) and ($length >= $min_length)) {
	$counter++;
	$seq_out->write_seq($seq_obj);
	}
}

print "$counter sequences added from $merge_fasta.\n";

sleep 3;

system("muscle -in $out_name -out Aligned_$out_name");

system("FastTreeMP -nosupport -out Tree_$out_name.newick  Aligned_$out_name");