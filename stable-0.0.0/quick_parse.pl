use warnings;
use strict;

my %counts;
my $filist = `ls /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_v9_20151209/PolIII_mapping/*.bam.parsed`;
#my $filist = `ls /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_v9_20151209/PolIII_mapping/*APOBEC*.bam.parsed`;
my @filelistt = split(/\s+/,$filist);
for my $file (@filelistt) {
    open(FL,$file);

    my @fisplit = split(/\//,$file);
    my $fi_short = $fisplit[$#fisplit];
    print STDERR "doing $fi_short\n";
    for my $line (<FL>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	if ($tmp[0] eq "TOTAL") {
	    $counts{$tmp[1]}{$file} = $tmp[3];
	}
    }
    close(FL);
}


for my $k (keys %counts) {
    my $outputfi = $k.".alldata_20160127.csv";
    
    open(OUT,">$outputfi");
    my @sorted = sort {$counts{$k}{$b} <=> $counts{$k}{$a}} keys %{$counts{$k}};
    for my $s (@sorted) {
	print OUT "$s\t$counts{$k}{$s}\n";
    }
    close(OUT);
}
