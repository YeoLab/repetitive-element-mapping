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
    print "uid $uid\n";
    $alreadydone_uids{$uid} = 1;
}
close(AL);


for my $uid (keys %submitted_uid) {
    print STDERR "skipping - alreadydone $uid\n" if (exists $alreadydone_uids{$uid});
    next if (exists $alreadydone_uids{$uid});
    print STDERR "submitting $uid\n";

    my $sh_out = "chrMparse.".$uid.".sh";
    open(SH,">$sh_out");
    print SH "\#\!\/bin\/sh\n";
    print SH "#PBS -N ".$sh_out."\n";
    print SH "#PBS -o ".$sh_out.".out\n";
    print SH "#PBS -e ".$sh_out.".err\n";
    print SH "#PBS -V\n";
    print SH "#PBS -l walltime=8:00:00\n";
    print SH "#PBS -l nodes=1:ppn=4\n";
    print SH "#PBS -A yeo-group\n";
    print SH "#PBS -q condo\n";
#    print SH "#PBS -q home-yeo\n";
    print SH "#PBS -M elvannostrand\@ucsd.edu\n";
    print SH "#PBS -m a\n";
    print SH "cd /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/\n";

    print SH "perl re_parse_chrMfinal_wrapper.pl $uid\n";
    close(SH);

    system("qsub $sh_out");
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
