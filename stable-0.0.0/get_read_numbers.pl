use warnings;
use strict;

my %read_info;
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20160906/";
my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
my %submitted_hash;
my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS_submittedonly.txt";
&read_manifest_fi($manifest_fi);

for my $dataset_label (keys %submitted_hash) {
    
    my $clip_parsedfi = $repeat_working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam.parsed";

    &read_parsed_output($clip_parsedfi,$dataset_label);
}

my @types = ("AllReads","UsableReads","GenomicReads","RepFamilyReads");

print "dataset";
for my $type (@types) {
    print "\t$type";
}
print "\n";

for my $dataset_label (keys %submitted_hash) {
    print "$dataset_label";
    for my $type (@types) {
	print "\t".$read_info{$dataset_label}{$type};
    }
    print "\n";
}



sub read_parsed_output {
    my $fi = shift;
    my $label = shift;
    my $fipath = $fi;
#    print STDERR "$fipath\n";
    open(F,$fipath); # || die "no $fipath\n";
    while (<F>) {
        my $line = $_;
        chomp($line);
        my @tmp = split(/\t/,$line);
        if ($tmp[0] eq "TOTAL") {
        } elsif ($tmp[0] eq "#DELETED") {
        } elsif ($tmp[0] eq "#READINFO") {
	    $read_info{$label}{$tmp[1]} = $tmp[2];
	}
    }
    close(F);
}



sub read_manifest_fi {
    my $submitted_list = shift;
    open(S,$submitted_list);
    for my $line (<S>) {
        chomp($line);
        my @tmp = split(/\t/,$line);
        my $uid = $tmp[0];
        next if ($uid eq "uID");
        my $rbp = $tmp[1];
        my $cellline = $tmp[2];
        my $original_rbp = $rbp;
        for my $rep ("01","02","INPUT") {
#       for my $rep ("01","02") {

            my $id = $uid."_".$rep."_".$rbp;
            my $actual_id = $uid."_".$rep."_".$original_rbp;

            if ($uid eq "203" && $rep eq "01") {
                $id = "271_01_HNRNPC";
            }
            if ($uid eq "203" && $rep eq "INPUT") {
                $id = "271_INPUT_HNRNPC";
            }

            if ($uid eq "353") {
                $rbp = "KHDRBS1-SAM68";
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
        

            $submitted_hash{$id}{flag} = 1;
            $submitted_hash{$id}{actual} = $actual_id;
        
        }
    }
    close(S);

}
