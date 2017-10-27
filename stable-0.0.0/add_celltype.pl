use warnings;
use strict;

my %uid2cellline;
my $fi = $ARGV[0];
open(F,$fi);
for my $line (<F>) {
    chomp($line);
    $line =~ s/\r//g;
    my @tmp = split(/\s+/,$line);
    my $uid = shift(@tmp);
    my $rbp = shift(@tmp);
    my $cellline = shift(@tmp);
    $uid2cellline{$uid} = $cellline;
}
close(F);

my $fi2 = $ARGV[1];	
open(G,$fi2);
for my $line (<G>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $id = $tmp[0];
    my @blah = split(/\_/,$id);
    $tmp[2] = $uid2cellline{$blah[0]};

    print "".join("\t",@tmp)."\n";
}
close(G);
