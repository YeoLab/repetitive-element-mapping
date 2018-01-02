use warnings;
use strict;

my %element_list;
my %entropy_hash;
my $file = "20170403_ENCODEFINAL_Submittedonly.txt.nopipes.csv";
open(F,$file);
for my $line (<F>) {
    chomp($line);

    my @tmp = split(/\t/,$line);
    next if ($tmp[0] eq "element");
    my ($id,$rep) = split(/\_/,$tmp[2]);

    my $element = $tmp[0];
    my $entropy = $tmp[8];

    $entropy_hash{$id}{$element}{$rep} = $entropy;
    
    $element_list{$element} = 1;

}
close(F);

my @desired_elements;
my @desired_elements2;
my %entropy_geometric_mean;
for my $element (keys %element_list) {
    my %entropy_per_element;
    for my $id (keys %entropy_hash) {
	
	$entropy_geometric_mean{$id}{$element} = log(sqrt( (2 ** $entropy_hash{$id}{$element}{"01"}) * (2 ** $entropy_hash{$id}{$element}{"02"}) ))/log(2);
	$entropy_per_element{$id} = $entropy_geometric_mean{$id}{$element};

    }

    my @sorted = sort {$entropy_per_element{$b} <=> $entropy_per_element{$a}} keys %entropy_per_element;
#    print "$element\tmax_entropy\t".$sorted[0]."\t".$entropy_per_element{$sorted[0]}."\n";
    push @desired_elements2,$element if ($entropy_per_element{$sorted[0]} > 0.1);
}

#eric note - made this manually by sorting the list of eleemnts above > 0.1
my $element_list_fi = "20170403_ENCODEFINAL_Submittedonly.txt.nopipes.csv.formatlab_0.2cutoff.ElementOrderManual.txt";
open(EL,$element_list_fi);
for my $line (<EL>) {
    chomp($line);
    push @desired_elements,$line;

}
close(EL);


my $outfi = $file.".formatlab_0.2cutoff.txt";
open(OUT,">$outfi");
print OUT "expt";
for my $element (@desired_elements) {
    print OUT "\t$element";
}
print OUT "\n";

for my $id (keys %entropy_hash) { 
    print OUT "$id";

    for my $element (@desired_elements) {
	print OUT "\t".$entropy_geometric_mean{$id}{$element};
    }
    print OUT "\n";
}
    
