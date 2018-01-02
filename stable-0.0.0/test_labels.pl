use warnings;
use strict;

my $fi1 = $ARGV[0];
my $fi2 = $ARGV[1];

my %uid2id;
my $match = 0;
my $mismatch = 0;
open(F,$fi1);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    $uid2id{$tmp[0]} = $tmp[1];
}
close(F);

open(G,$fi2);
for my $line (<G>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    if (exists $uid2id{$tmp[0]}) {
	if ($uid2id{$tmp[0]} == $tmp[1]) {
	    $match++;
	} else {
	    $mismatch++;
	}
    }
}
close(G);

print "match $match \n mismatch $mismatch\n";
