use warnings;
use strict;
use Getopt::Long;

my $chunk_size = 1000;
my $extension = "fa";
my $directory;
my $threads = 48;
my $pplacer_threads = 12;

GetOptions(
'chunk_size=s' => \$chunk_size,
'extension=s' => \$extension,
'directory=s' => \$directory,
'threads=s' => \$threads,
'pplacer_threads=s' => \$pplacer_threads,
);


my $content = `ls $directory | grep $extension`;
my @files = split /\n/, $content;

my $posit = 0;

my $chunk_count = 0;
system("mkdir Chunk_$chunk_count");

my $file_count = 0;

foreach my $file (@files) {
	$file_count++;
	system("cp $directory/$file Chunk_$chunk_count/");
	if (($file_count == $chunk_size) or ($file_count == @files)) {
		system("checkm lineage_wf -x $extension -t $threads --pplacer_threads $pplacer_threads --tab_table --file Chunk_$chunk_count.Checkm_Results.tsv Chunk_$chunk_count/ Checkm_Output_Chunk_$chunk_count/");
		$chunk_count++;
		$file_count = 0;
		system("mkdir Chunk_$chunk_count");
	}
}

system("checkm lineage_wf -x $extension -t $threads --pplacer_threads $pplacer_threads --tab_table --file Chunk_$chunk_count.Checkm_Results.tsv Chunk_$chunk_count/ Checkm_Output_Chunk_$chunk_count/");