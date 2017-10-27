use warnings;
use strict;

my %uid2info;
my %unique_reads;
my %read_count;
my %read_info;
my %submitted_hash;

my %submitted_uid;
my $repeat_working_directory = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/";

#my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/integrated_paper_figures/DatasetListing_FINAL.txt.ENCODEFINAL.txt";
my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/ALLDATASETS_submittedonly.txt";
&read_manifest_fi($manifest_fi);

my %data;
my %alreadydone_uids;
my $alreadydone_fi = "20170505.ALLENCODEinclnotsubmitted.txt.nopipes.txt.submitted_only.chrM_reparsed.txt";
open(AL,$alreadydone_fi);
for my $line (<AL>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my @ids = split(/\_/,$tmp[0]);
    my $rbp = pop(@ids);
    my $rep = pop(@ids);
    my $uid = join("_",@ids);

    $alreadydone_uids{$uid} = $line;
    my $element = $tmp[3];
    $data{$uid}{$rep}{$element} = $line;
}
close(AL);

for my $uid (keys %submitted_uid) {

    if (exists $alreadydone_uids{$uid}) {
	for my $rep ("01","02") {
	    for my $element ("chrM_reprocess_+strand","chrM_reprocess_-strand") {
		if (exists $data{$uid}{$rep}{$element}) {
		    print "".$data{$uid}{$rep}{$element}."\n";
		} else {
		    print STDERR "missing $uid $rep $element\n";
		}
	    }
	}
	next;
    }

    my $outfi = $repeat_working_directory.$uid.".chrMtmpout.txt";
    if (-e $outfi) {
	open(OUT,$outfi) || die "couldn't open $outfi\n";
	for my $line (<OUT>) {
	    chomp($line);
	    my @tmp = split(/\t/,$line);

	    my @ids = split(/\_/,$tmp[0]);
	    my $rbp = pop(@ids);
	    my $rep = pop(@ids);
	    my $uid = join("_",@ids);
	    $alreadydone_uids{$uid} = $line;
	    my $element= $tmp[3];
	    $data{$uid}{$rep}{$element} = $line;
	}
	close(OUT);

	for my $rep ("01","02") {
            for my $element ("chrM_reprocess_+strand","chrM_reprocess_-strand") {
		if (exists $data{$uid}{$rep}{$element}) {
                    print "".$data{$uid}{$rep}{$element}."\n";
		} else {
                    print STDERR "missing $uid $rep $element\n";
		}
            }
        }
	next;
    }

    print STDERR "still missing data for $uid\n";
}



sub read_manifest_fi {
    my $submitted_list = shift;
    my @rep_listing;
    my $type_flag;
    open(S,$submitted_list);
    for my $line (<S>) {
        chomp($line);
        $line =~ s/\r//g;
        my @tmp = split(/\s+/,$line);
        my $uid = $tmp[0];
        next unless ($uid);
        next if ($uid eq "uID");
        my $rbp = $tmp[1];
        my $cellline = $tmp[2];
        my $original_rbp = $rbp;
        if (scalar(@tmp) == 6) {
            $type_flag = "two_replicate_ENCODEstyle";
        } elsif (scalar(@tmp) == 5) {
            $type_flag = "one_replicate";
        }

        if ($type_flag eq "two_replicate_ENCODEstyle") {
            @rep_listing = ("_01","_02","_INPUT");
        } elsif ($type_flag eq "one_replicate") {
            @rep_listing = ("_01","_INPUT");
        } else {
            print STDERR "TYPE flag is not set properly!!!!\n";
        }

        my $id = $tmp[1]."-".$tmp[2];

        $submitted_uid{$uid} = $id;
        $uid2info{$uid}{rbp} = $rbp;
        $uid2info{$uid}{cellline} = $cellline;
    }
    close(S);

}
