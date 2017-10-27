use warnings;
use strict;

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

for my $uid (keys %submitted_uid) {
#    next unless ($uid eq "206");
    for my $rep ("01","02","INPUT") {
	my $dataset_label = $uid2repids{$uid}{$rep};

#for my $dataset_label (keys %submitted_hash) {
#    next unless ($dataset_label =~ /437_/);
	my $clip_parsedfi = $repeat_working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam.parsed";
	
	&read_parsed_output($clip_parsedfi,$dataset_label);
	
#	    print "".$uid."\t".$rep."\t".$dataset_label."\t".($read_count{$dataset_label} / $read_info{$dataset_label}{"UsableReads"})."\n";
#	    print "".$dataset_label."\tusable ".$read_info{$dataset_label}{"UsableReads"}."\t unique ".$unique_reads{$dataset_label}."\n";
	print "".$uid."\t".$rep."\t".$dataset_label."\t".$dataset_label."\tinputreads\t".$read_info{$dataset_label}{"AllReads"}."\tusable\t".$read_info{$dataset_label}{"UsableReads"}."\tunique_genome_nonrep\t".$unique_reads{$dataset_label}."\n"; 
    }
}



sub read_parsed_output {
    my $fi = shift;
    my $label = shift;
    my $fipath = $fi;
    print STDERR "$fipath\n";
    open(F,$fipath); # || die "no $fipath\n";
    while (<F>) {
        my $line = $_;
        chomp($line);
        my @tmp = split(/\t/,$line);
        if ($tmp[0] eq "TOTAL") {
	    if ($tmp[1] =~ /unique/) {
		$unique_reads{$label} += $tmp[2];
	    }
	    next if ($tmp[1] =~ /unique/);
	    next if ($tmp[1] =~ /chrM/ || $tmp[1] =~ /H1RNA/ || $tmp[1] =~ /MRP/);
	    $read_count{$label} += $tmp[2];
        } elsif ($tmp[0] eq "#DELETED") {
        } elsif ($tmp[0] eq "#READINFO") {
            $read_info{$label}{$tmp[1]} = $tmp[2];
	}
    }
    close(F);
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
    }
    close(S);

}

