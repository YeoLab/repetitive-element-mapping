use warnings;
use strict;

my %hash;
my %sums;

my $fi = "20160508_summary.submitted_only";
open(F,$fi);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    if ($tmp[1] =~ /\_antisense/) {
	next if ($tmp[3] eq "Not_found");
	next if ($tmp[1] eq "Low_complexity_antisense" || $tmp[1] eq "Simple_repeat_antisense");
	$hash{$tmp[1]} = 1;
	$sums{$tmp[0]}{clip} += $tmp[4];
	$sums{$tmp[0]}{input} += $tmp[6];
    }
}
close(F);

for my $k (keys %hash) {
#    print "$k\n";
}


for my $k (keys %sums) {
    print "$k\t".$sums{$k}{clip}."\t".$sums{$k}{input}."\n";
}
