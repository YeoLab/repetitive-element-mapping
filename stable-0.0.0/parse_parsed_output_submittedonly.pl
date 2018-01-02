use warnings;
use strict;

my %submitted_hash;
my $submitted_list = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/ALLDATASETS_submittedonly.txt";
#my $submitted_list = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_CLIPperv2_20161120/ALLDATASETS_submittedonly.txt";
#my $submitted_list = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS_submittedonly.txt";
#my $submitted_list = "/home/elvannostrand/data/clip/CLIPseq_analysis/ALLCLIP_v12_20160112/encode_v12_filelist.allencode_20160314.txt.submittedHepG2K562_20160316.txt.fixedgenename_w4000label.txt";
open(S,$submitted_list);
for my $line (<S>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $uid = $tmp[0];
    next if ($uid eq "uID");
    my $rbp = $tmp[1];
    my $orig_rbp = $rbp;
    my $cellline = $tmp[2];
    my $original_rbp = $rbp;
#    for my $rep ("01","02","INPUT") {
    for my $rep ("01","02") {

	my $id = $uid."_".$rep."_".$rbp;
	my $actual_id = $uid."_".$rep."_".$original_rbp;

	if ($uid eq "203" && $rep eq "01") {
	    $id = "271_01_HNRNPC";
	}

	if ($uid eq "353") {
	    $rbp = "KHDRBS1-SAM68";
	    $id = $uid."_".$rep."_".$rbp;
	}

	if ($uid eq "516" || $uid eq "553") {
	    $rbp = "BOP1";
	    $id = $uid."_".$rep."_".$rbp;
	}

	if ($uid eq "366") {
	    $rbp = "TRNC6A";
	    $id = $uid."_".$rep."_".$rbp;
	}

	if ($uid eq "447") {
	    $rbp = "RO60-TROVE2";
	    $id = $uid."_".$rep."_".$rbp;
	}

	if ($uid eq "235x4000") {
	    $id = "235_".$rep."_4000_".$rbp;
	}

	if ($uid eq "390x4000") {
            $id = "390_".$rep."_4000_".$rbp;
	}
	
	if ($uid eq "632x") {
	    $id = "632_".$rep."_bc2rev_".$rbp;
	}
	if ($uid eq "285") {
	    $id = "285_".$rep."_4000_".$rbp;
	}

#	print "id $id\n";
	$submitted_hash{$id}{flag} = 1;
	$submitted_hash{$id}{actual} = $actual_id;
	$submitted_hash{$id}{info} = $orig_rbp."\t".$cellline;
    }
}
close(S);


my %output_hash;
my %output_ids;
my $summary_fi = $ARGV[0];
my $submitted = $summary_fi.".submitted_only";
my $submitted_matlab = $summary_fi.".submitted_only.matlab.csv";
open(OUT,">$submitted");
open(OUTMAT,">$submitted_matlab");
open(S,$summary_fi);
for my $line (<S>) {
    chomp($line);
    
    my @tmp = split(/\t/,$line);
    my $id = $tmp[1];

    if (exists $submitted_hash{$id}) {
	print OUT "".$submitted_hash{$id}{actual}."\t".$submitted_hash{$id}{info}."\t$line\n";
	if ($tmp[6]) {
	    $output_hash{$tmp[0]}{$id} = $tmp[6];
	} else {
	    $output_hash{$tmp[0]}{$id} = 1;
	}
	

	$output_ids{$id} = 1;
	$submitted_hash{$id}{flag} = 2;
    }

}
close(S);

print OUTMAT "element";
for my $id (keys %output_ids) {
    print OUTMAT "\t$id";
}
print OUTMAT "\n";

for my $element (keys %output_hash) {
    print OUTMAT "$element";
    for my $id (keys %output_ids) {
	print OUTMAT "\t".$output_hash{$element}{$id};
    }
    print OUTMAT "\n";
}

for my $id (keys %submitted_hash) {
    if ($submitted_hash{$id}{flag} == 2) {
    } else {
	print STDERR "missing $id\n";
    }
}
