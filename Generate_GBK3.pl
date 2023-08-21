use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SeqFeatureI;

my $nucleotide_file;
my $protein_file;
my $gff_file;
my $annot;
my $max_evalue = 0.001;
my $min_bitscore = 30;
my $min_algn = 0;
my $min_id = 0;
my $explode;
my $skip_unk;

GetOptions(
'annot=s' => \$annot,
'nucleotide_file=s' => \$nucleotide_file,
'protein_file=s' => \$protein_file,
'gff_file=s' => \$gff_file,
'min_bitscore=s' => \$min_bitscore,
'max_evalue=s' => \$max_evalue,
'min_alignment=s'=> \$min_algn, 
'min_identity=s' => \$min_id,
'explode' => \$explode,
'skip_unk' => \$skip_unk,
);

my %anots;
my %seqs2feats;
my %feat_counts;

print "Parsing $protein_file\n";
my $seq_in = Bio::SeqIO->new('-file' => "< $protein_file", '-format' => "Fasta");
while (my $seq_obj = $seq_in->next_seq) {
	my $id = $seq_obj->id();
	my $peg_seq = $seq_obj->seq;
	$peg_seq =~ s/\*$//;
	$anots{"Sequence"}{$id}  = $peg_seq;
}

if ($annot) {	
	print "Parsing $annot\n";
	open INPUT, "< $annot" or die "$!";
	while (my $line = <INPUT>) {
		chomp $line;
		my @values = split /\t/, $line;
		
		my $peg = $values[0];
		next if (defined $anots{"Product"}{$peg});
		
		my $pass = check_cutoffs(@values);
		
		next unless ($pass);
		my $taxon = $values[-1];
		
		my @info = split /(\s)+n=(\d)+(\s)+Tax=/, $values[12];
		my $function = $info[0];
		$function =~ s/n=(\d)+(\s)+Tax=//;
		$function =~ s/^(\w)+_(\w)+(\s)//;
		next if (($skip_unk) and ($function =~ /(Uncharacterized)|(hypothetical)/));
		$anots{"Product"}{$peg} = $function;
		$anots{"Taxon"}{$peg} = $taxon;
	}
}

print "Parsing $gff_file\n";
my $gffio = Bio::Tools::GFF->new(-file => "< $gff_file", -gff_version => 3);
while(my $feature = $gffio->next_feature()) {
	my $type = $feature->primary_tag();
	my $oid = $feature->seq_id();
	$feat_counts{$oid}++;
	
	#my $peg_num = $feature->get_tag_values('ID');
	#$peg_num =~ s/^ID(.)+_//;
	my $peg_num = $feat_counts{$oid};
	my $peg = $oid."_".$peg_num;
	
	#print "$peg ($peg_num) from $oid\n"; #sleep 1;
	
	$feature->add_tag_value('product',$peg);
	$feature->add_tag_value('translation',$anots{"Sequence"}{$peg}) if (defined $anots{"Sequence"}{$peg});
	$feature->add_tag_value('function',$anots{"Product"}{$peg}) if (defined $anots{"Product"}{$peg});
	if (defined $anots{"Taxon"}{$peg}) {$feature->add_tag_value('taxon',$anots{"Taxon"}{$peg})}
	$seqs2feats{$oid}{"$feat_counts{$oid}"} = $feature;
	}
	
print "Parsing $nucleotide_file\n";

if ($explode) {
	$seq_in = Bio::SeqIO->new('-file' => "< $nucleotide_file", '-format' => "Fasta");
	while (my $seq_obj = $seq_in->next_seq) {
		my $oid = $seq_obj->id();
		my $desc = $seq_obj->description();
		$seq_obj->description("Unknown") if ($desc eq "");
		my $seq_out = Bio::SeqIO->new('-file' => "> $oid.gbk", '-format' => "Genbank");
		if (defined $feat_counts{$oid}) {
			foreach my $feat_num (1..$feat_counts{$oid}) {
				my $feature = $seqs2feats{$oid}{$feat_num};
				$seq_obj->add_SeqFeature($feature);
			}
		}
		$seq_out->write_seq($seq_obj);
	}	
} else {
	$seq_in = Bio::SeqIO->new('-file' => "< $nucleotide_file", '-format' => "Fasta");
	my $seq_out = Bio::SeqIO->new('-file' => "> Test.gbk", '-format' => "Genbank");
	while (my $seq_obj = $seq_in->next_seq) {
		my $oid = $seq_obj->id();
		my $desc = $seq_obj->description();
		$seq_obj->description("Unknown") if ($desc eq "");
		if (defined $feat_counts{$oid}) {
			foreach my $feat_num (1..$feat_counts{$oid}) {
				my $feature = $seqs2feats{$oid}{$feat_num};
				$seq_obj->add_SeqFeature($feature);
			}
		}
		$seq_out->write_seq($seq_obj);
	}
}

sub check_cutoffs {
#Receive, as an array, the values of every line from the blast file. Check if te values match the user established criteria for minimum identity, minimum alignment, minimum bitscore, maximum evalue.
	my $identity = $_[2];
	my $algn_length = $_[3];
	my $evalue = $_[10];
	my $bitscore = $_[11];
	return 0 if ($identity < $min_id);
	return 0 if ($algn_length < $min_algn);
	return 0 if ($evalue > $max_evalue);
	return 0 if ($bitscore < $min_bitscore);
	return 1;
}

