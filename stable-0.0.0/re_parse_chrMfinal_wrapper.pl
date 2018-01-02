use warnings;
use strict;

my %uid2info;
my %unique_reads;
my %read_count;
my %read_info;
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20160906/";
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
my %submitted_hash;

my %submitted_uid;
my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20170304.txt";

my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_GRCh38_v1.txt";

#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
my $repeat_working_directory = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/";
my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/integrated_paper_figures/DatasetListing_FINAL.txt.ENCODEFINAL.txt";

my %clip2input;
my %uid2repids;
open(IN,$inputpairingfi);
for my $line (<IN>) {
    chomp($line);
    $line =~ s/\r//g;

    my @tmp = split(/\s+/,$line);
    
    if ($tmp[3] eq "LNG" || $tmp[3] eq "TAG" || $tmp[3] eq "SingleRep" || $tmp[3] eq "OneRep") {
#       $clip2input{$tmp[1]} = $tmp[2];
        next;
    } else {
        next unless ($tmp[4]);
        $clip2input{$tmp[1]} = $tmp[3];
        $clip2input{$tmp[2]} = $tmp[3] unless ($tmp[2] eq "NoInput");

        my $uid = $tmp[4];
        $uid2repids{$uid}{"01"} = $tmp[1];
        $uid2repids{$uid}{"02"} = $tmp[2];
        $uid2repids{$uid}{"INPUT"} = $tmp[3];
        
    }
}
close(IN);

&read_manifest_fi($manifest_fi);

my $desired_uid = $ARGV[0];


for my $uid (keys %submitted_uid) {

    next unless ($uid eq $desired_uid);
    my %data;
    for my $rep ("01","02","INPUT") {
        my $dataset_label = $uid2repids{$uid}{$rep};

        my $clip_parsedfi = $repeat_working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam.parsed";
        my $clip_samfi = $repeat_working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam.gz";

	my $outline = `perl re_parse_chrMfinal.pl $clip_samfi`;
	print STDERR "outline\n$outline\n";
	my @outlines = split(/\n/,$outline);
	for my $outline_split (@outlines) {
	    my ($id,$count,$total_reads,$rpr) = split(/\t/,$outline_split);
	
	    $data{$id}{count}{$rep} = $count;
	    $data{$id}{total_reads}{$rep} = $total_reads;
	    $data{$id}{rpr}{$rep} = $rpr;
	}
    }

    my $outfi = $repeat_working_directory.$uid.".chrMtmpout.txt";
    open(OUT,">$outfi");
    

    for my $element ("chrM_reprocess_+strand","chrM_reprocess_-strand") {
#    for my $element (keys %data) {

	for my $rep ("01","02") {
	    print STDERR "doing uid $uid rep $rep element $element \n";
	    my $rep_finalid = $uid."_".$rep."_".$uid2info{$uid}{rbp};
	    my $rep_out = $uid2info{$uid}{rbp}."-".$uid2info{$uid}{cellline}."_".$rep;
	    my $rep_dataset = $uid2repids{$uid}{$rep};
	    my $repcount = 0;
	    $repcount = $data{$element}{count}{$rep} if (exists $data{$element}{count}{$rep});
	    my $inputcount = 0;
	    $inputcount = $data{$element}{count}{"INPUT"} if (exists $data{$element}{count}{"INPUT"});
	    my $reprpr;
	    if (exists $data{$element}{rpr}{$rep}) {
		$reprpr = $data{$element}{rpr}{$rep};
	    }
	    my $inputrpr = 0;
	    if (exists $data{$element}{rpr}{"INPUT"}) {
		$inputrpr = $data{$element}{rpr}{"INPUT"};
	    }
	    my $norm_inputrpr = $inputrpr;
	    unless ($inputrpr > 0) {
		$norm_inputrpr = 1 / $data{$element}{total_reads}{"INPUT"};
	    }
	    my $foldenr = $reprpr / $norm_inputrpr;
	    my $entropy = $reprpr * log ($reprpr / $norm_inputrpr) / log(2);
	    
	    print OUT "".$rep_finalid."\t".$uid2info{$uid}{rbp}."\t".$uid2info{$uid}{cellline}."\t".$element."\t".$rep_dataset."\t".$rep_out."\t".$repcount."\t".$reprpr."\t".$inputcount."\t".$inputrpr."\t".$foldenr."\t".$entropy."\n";
	}
    }
    close(OUT);
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
