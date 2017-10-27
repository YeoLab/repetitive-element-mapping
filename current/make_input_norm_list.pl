use warnings;
use strict;

my %already_done;
my $current_inputnorm_list = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20160726.txt";
open(CUR,$current_inputnorm_list);
for my $line (<CUR>) {
    chomp($line);
    $line =~ s/\r//g;
    my @tmp = split(/\t/,$line);
    for my $i (0..3) {
#	print "$i\t$tmp[$i]X\n";
    }

    if ($tmp[3] eq "LNG" || $tmp[3] eq "TAG" || $tmp[3] eq "SingleRep" || $tmp[3] eq "OneRep") {
	$already_done{$tmp[1]} = 1;
	$already_done{$tmp[2]} = 1;	
    } elsif ($tmp[0] eq "CORRECT" || $tmp[0] eq "manual") {
	$already_done{$tmp[1]} = 1;
        $already_done{$tmp[2]} = 1;
	$already_done{$tmp[3]} = 1;
#	print "added all 3\n";
    } else {
	print STDERR "didnt' read this line properly $line\n";
    }
}
close(CUR);

my %hash;
my %rephash;
my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_GRCh38_v1.txt";
#my $merging_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/encode_GRCh38_v1.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_master.txt";
open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $fastqs = shift(@tmp);
    my $genome = shift(@tmp);
    my $dataset_label = shift(@tmp);
#    next unless ($dataset_label =~ /^KB/);

    $hash{$dataset_label} = 1;

    if ($dataset_label =~ /\_.+\_.+\_/) {
	print "skipping\t$dataset_label\n";
	next;
    }
    if ($dataset_label =~ /^(\d+)\_(0\d)\_(\w+)$/ || $dataset_label =~ /^(\d+)\_(INPUT)\_(\w+)$/) {
	my ($uid,$rep,$rbp) = ($1,$2,$3);
	$rephash{$uid}{$rep} = $dataset_label;
    }

    if ($dataset_label =~ /^(LNG\d+\-\w)\_(CLIP)\_(\w+)$/ || $dataset_label =~ /^(LNG\d+\-\w)\_(INPUT)\_(\w+)$/) {
	my ($uid,$rep,$rbp) = ($1,$2,$3);
        $rephash{$uid}{$rep} = $dataset_label;
    }

    if ($dataset_label =~ /^(TAG\d+)\_(CLIP)\_(\w+)$/ || $dataset_label =~ /^(TAG\d+)\_(INPUT)\_(\w+)$/) {
	my ($uid,$rep,$rbp) = ($1,$2,$3);
	$rephash{$uid}{$rep} = $dataset_label;
    } 

}
close(M);

for my $uid (keys %rephash) {
    if ($uid =~ /^\d+$/) {
	if (exists $rephash{$uid}{"01"} && exists $rephash{$uid}{"02"} && exists $rephash{$uid}{"INPUT"}) {
	    for my $rep ("01","02","INPUT") {
		$hash{$rephash{$uid}{$rep}} = 2;
	    }
	    next if (exists $already_done{$rephash{$uid}{"01"}} && exists $already_done{$rephash{$uid}{"02"}} && exists $already_done{$rephash{$uid}{"INPUT"}});

	    print "CORRECT\t".$rephash{$uid}{"01"}."\t".$rephash{$uid}{"02"}."\t".$rephash{$uid}{"INPUT"}."\n";
	    

	}
    } elsif ($uid =~ /^LNG/) {
	if (exists $rephash{$uid}{"CLIP"} && exists $rephash{$uid}{"INPUT"}) {
	    for my $rep ("CLIP","INPUT") {
		$hash{$rephash{$uid}{$rep}} = 2;
	    }
	    next if (exists $already_done{$rephash{$uid}{"CLIP"}} && exists $already_done{$rephash{$uid}{"INPUT"}});
	    print "CORRECT\t".$rephash{$uid}{"CLIP"}."\t".$rephash{$uid}{"INPUT"}."\tLNG\n";

	}
    } elsif ($uid =~ /^TAG/) {
	if (exists $rephash{$uid}{"CLIP"} && exists $rephash{$uid}{"INPUT"}) {
            for my $rep ("CLIP","INPUT") {
                $hash{$rephash{$uid}{$rep}} = 2;
            }
	    next if (exists $already_done{$rephash{$uid}{"CLIP"}} && exists $already_done{$rephash{$uid}{"INPUT"}});
            print "CORRECT\t".$rephash{$uid}{"CLIP"}."\t".$rephash{$uid}{"INPUT"}."\tTAG\n";
        }
    } else {
	print "uid doesn't match $uid\n";
    }
}

for my $k (keys %hash) {
    next if ($hash{$k} == 2);
    print "not_found\t$k\n";
}
