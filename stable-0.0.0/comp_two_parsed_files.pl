use warnings;
use strict;

my $fi1 = $ARGV[0];
my $fi2 = $ARGV[1];

my %hash;
&parse_file($fi1,"fi1");
&parse_file($fi2,"fi2");

for my $id (keys %hash) {
    
    $hash{$id}{"fi1"} = 0 unless (exists $hash{$id}{"fi1"});
    $hash{$id}{"fi2"} = 0 unless (exists $hash{$id}{"fi2"});
    
    print "$id\t".$hash{$id}{"fi1"}."\t".$hash{$id}{"fi2"}."\n";
}





sub parse_file {
    my $fi = shift;
    my $label = shift;


    open(F,$fi);
    for my $line (<F>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	my $id = $tmp[0]."-".$tmp[1];

	next unless ($tmp[0] eq "TOTAL");
	$hash{$id}{$label} = $tmp[2];
	
    }
    close(F);
}
