use warnings;
use strict;

my @files;

if (@ARGV == 1) {
	#print "Looking for files with extension $ARGV[0] in current directory";
	@files = glob("*$ARGV[0]");
} elsif (@ARGV > 1) {
	@files = @ARGV;
}


print "File\tContigs\tBases\tMax\tN50\tN90\n";

foreach my $file (@files) {
	my $result = `perl /mnt/lustre/bio/users/fcoutinho/Scripts/contig-stats.pl $file`;
	$result .= "NA\n" unless ($result =~ /\n$/);
	print "$file\t$result";
}
