#/usr/bin/perl -w
#source http://seqanswers.com/forums/showthread.php?t=2766  author = BENM
#count and maximum length added by Aaron Liston June 8, 2011
use strict;

my ($len,$total,$contigs)=(0,0,1);
my @x;
while(<>){
	if(/^[\>\@]/){
		if($len>0){
			$total+=$len;
        		$contigs ++;
			push @x,$len;
		}
		$len=0;
	}
	else{
		s/\s//g;
		$len+=length($_);
	}
}
if ($len>0){
	$total+=$len;
	push @x,$len;
}
@x=sort{$b<=>$a} @x; 
my $max_value = $x[0];
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){
	$count+=$x[$j];
	if (($count>=$total/2)&&($half==0)){
		print "$contigs\t";
		print "$total\t";
		print "$max_value\t";
		print "$x[$j]\t";
		$half=$x[$j]
	}elsif ($count>=$total*0.9){
		print "$x[$j]\n";
		exit;
	}
}
