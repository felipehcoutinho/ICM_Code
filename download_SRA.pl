use warnings;
use strict;
use Getopt::Long;

my $list_file;


GetOptions(
'list_file=s' => \$list_file,
);

open IN, "< $list_file" or die "$!";

while (my $line = <IN>) {
	chomp $line;
	my $sample_id = $line;
	my $prefetch_command = "/mnt/smart/scratch/vir/felipe/Software/sratoolkit.3.1.1-ubuntu64/bin/prefetch $sample_id --output-directory .";
    system("$prefetch_command");
	system("/mnt/smart/scratch/vir/felipe/Software/sratoolkit.3.1.1-ubuntu64/bin/fasterq-dump --split-files $sample_id");
	system("rm -fr $sample_id");
}
	
