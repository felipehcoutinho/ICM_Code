use warnings;
use strict;
use Bio::SeqIO;

my $file_count = 0;
my $seq_count = 0;
my $in_obj = Bio::SeqIO->new(-file => "< $ARGV[0]", -format => 'fasta');
my $out_obj = Bio::SeqIO->new(-file => "> Part_$file_count.fasta", -format => "fasta");
my $size = $ARGV[1];

while (my $seq = $in_obj->next_seq) {
	if ($seq_count == $size) {$seq_count = 0; $file_count++; $out_obj = Bio::SeqIO->new(-file => "> Part_$file_count.fasta", -format => "fasta")}
	$out_obj->write_seq($seq);
	$seq_count++;
	}


