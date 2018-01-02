use warnings;
use strict;


my %cutoffs;
$cutoffs{"RNA18S"} = 3;
$cutoffs{"RNA28S"} = 3;
$cutoffs{"rRNA_extra"} = 3;
$cutoffs{"RNA5-8S"} = 5;

my $fi = "20170228_parsedoutput.ALLdatasets.csv.nopipes.csv";
open(F,$fi);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    if (exists $cutoffs{$tmp[0]}) {
	
	if ($tmp[6] >= $cutoffs{$tmp[0]}) {
	    print "$tmp[1]\n";
	}
    }
}
close(F);
