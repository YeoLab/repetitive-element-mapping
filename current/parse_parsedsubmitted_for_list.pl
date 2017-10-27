use warnings;
use strict;

my %hash;
my $data_file = "20170203_parsedoutput.csv.submitted_only.nopipes.csv";
open(D,$data_file);
for my $line (<D>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    
    $hash{$tmp[0]}{$tmp[3]}{fc} = $tmp[9];
}
close(D);

my @desired_elements = ("RNU1","RNU2","RNU6","RNU11","RNU12","RNU6ATAC");

my %alreadydone;
my $filist_fi = $ARGV[0];
open(FL,$filist_fi);
for my $line (<FL>) {
    chomp($line);
    next unless ($line);
    my ($uid,$rep,$rbp) = split(/\_/,$line);

    my %rep;
    $rep{"1"} = $uid."_01_".$rbp;
    $rep{"2"} = $uid."_02_".$rbp;

    next if (exists $alreadydone{$uid});
    $alreadydone{$uid} = 1;

    for my $replicate (1,2) {
	
	print "".$rep{$replicate}; 
	for my $element (@desired_elements) {
	    print "\t$element\t".$hash{$rep{$replicate}}{$element}{fc};
	}
	print "\n";
    }

}
close(FL);



