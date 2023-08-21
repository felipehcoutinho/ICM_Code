use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $fasta;

GetOptions(
'fasta=s' => \$fasta,
) or die "Missing arguments!\n";



my $seq_in = Bio::SeqIO->new('-file' => "< $fasta", '-format' => "Fasta");

my $counter = 0;

while (my $seq_obj = $seq_in->next_seq) {
	my $id = $seq_obj->id();
	my $seq_out = Bio::SeqIO->new('-file' => "> $id.fasta", '-format' => "Fasta");
	$counter++;
	$seq_out->write_seq($seq_obj);
	}
	
print "$counter sequences processed.\n";