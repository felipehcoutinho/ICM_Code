use warnings;
use strict;
use Getopt::Long;

my $abundance_file;
my $seq_info_file;
my $sample_info_file;
my $out_file = "RPKM.tsv";
my $remove_extension;

GetOptions(
'abundance_file=s' => \$abundance_file,
'seq_info_file=s' => \$seq_info_file,
'sample_info_file=s' => \$sample_info_file,
'out_file=s' => \$out_file,
'remove_extension' => \$remove_extension
);

my %params = ("file" => $seq_info_file);
my $seq_info_ref = read_table(\%params);
my %seq_info = %$seq_info_ref;

%params = ("file" => $sample_info_file);
my $sample_info_ref = read_table(\%params);
my %sample_info = %$sample_info_ref;

if ($remove_extension) {
	my %new_sample_info;
	foreach my $var (keys %sample_info) {
		my $level2_ref = $sample_info{$var};
		my %level2 = %$level2_ref;
		foreach my $sample_id (keys %level2) {
			my $sample_nid = $sample_id;
			$sample_nid =~ s/_(\d)\.fastq$//;
			$sample_nid =~ s/_R(\d)\.Paired\.fastq$//;
			#print "Changed $sample_id to $sample_nid\n";
			$new_sample_info{$var}{$sample_nid} = $sample_info{$var}{$sample_id};
		}
	}
	%sample_info = %new_sample_info;
}

%params = ("file" => $abundance_file);
my $abundance_ref = read_table(\%params);
my %abundance = %$abundance_ref;

my %rpkm;

foreach my $sample (keys %abundance) {
	print "Processing sample $sample\n";
	my $level2_ref = $abundance{$sample};
	my %level2 = %$level2_ref;
	foreach my $seq (keys %level2) {
		my $seq_length = $seq_info{"Length"}{$seq};
		my $sample_read_count = $sample_info{"R1_Read_Count"}{$sample};
		#print "Sample: $sample (Reads = $sample_read_count). Sequence: $seq (Length = $seq_length). Abundance = $abundance{$sample}{$seq}\n";
		my $rpkm_val = $abundance{$sample}{$seq} / (($seq_length / 1000) * ($sample_read_count / 1_000_000));
		$rpkm{$sample}{$seq} = $rpkm_val;
	}
}

print "Printing output to $out_file\n";
%params = ("matrix" => \%rpkm, "file" => $out_file);
print_matrix(\%params);

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
            if (defined $matrix{$col}{$row}) {print OUTPUT "$sep"."$matrix{$col}{$row}"} else {print OUTPUT "$sep"."NA"};
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
