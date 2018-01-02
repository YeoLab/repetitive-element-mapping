use warnings;
use strict;

my $fi = $ARGV[0];
open(F,$fi);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    my $join = join("\t",@tmp);
    print "$join\n";
}
close(F);
