use warnings;
use strict;

my $desired_type = "L1||L1HS";
#my $desired_type = "YRNA||HY5_uniquegenomic";
#my $desired_type = "antisense_L1||antisense_L1HS";
#my %desired_genes;
#my $desired_type 


my %pseudocount_val;
my %counts;

my %clip2input;
my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20160726.txt";
#my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepMappingInputPairing_20160409.txt";
open(IN,$inputpairingfi);
for my $line (<IN>) {
    chomp($line);
    $line =~ s/\r//g;

    my @tmp = split(/\t/,$line);
    
    if ($tmp[3] eq "LNG" || $tmp[3] eq "TAG" || $tmp[3] eq "SingleRep" || $tmp[3] eq "OneRep") {
	$clip2input{$tmp[1]} = $tmp[2];
    } else {
	$clip2input{$tmp[1]} = $tmp[3];
        $clip2input{$tmp[2]} = $tmp[3] unless ($tmp[2] eq "NoInput");
    }

}
close(IN);



my %submitted_hash;
my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS_submittedonly.txt";
&read_manifest_fi($manifest_fi);


my $manifestfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS.txt";
#my $manifestfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ALLCLIP_v12_20160112/encode_v12_filelist.allencode_20160314.txt.fixedgenename.txt";

#my $repeat_working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/ALLCLIP_v12_20160112/repeatmapping/";
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_v12_repeatmapping/";
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_v12_20160506repeatmapping/";
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repeatfamilymapping_20160723b/";
my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20160906/";


my %missing_dataset;

my $counter=0;
my %hash;

my $current_file_location_prefix = "/projects/ps-yeolab3/encode/analysis/encode_master/";

my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_GRCh38_v1.txt";


open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $fastqs = shift(@tmp);
    my $genome = shift(@tmp);
    my $dataset_label = shift(@tmp);
#    next unless ($dataset_label =~ /^KB/);
    push @{$hash{$dataset_label}},$fastqs;
#    last if ($counter > 10);
    $counter++;
}
close(M);
my %usable_readnum;
my %alreadydone;
for my $dataset_label (keys %hash) {
    next unless (exists $submitted_hash{$dataset_label});

    for my $fastq (@{$hash{$dataset_label}}) {
	my ($fastq1,$fastq2) = split(/\;/,$fastq);
	
	for ($fastq1,$fastq2) {
	    my @fastq_split = split(/\//,$_);
	    my $fastq_short = $fastq_split[$#fastq_split];
            
	    $fastq_short =~ s/\.fastq\.gz$//;
	    my $orig_fi = $current_file_location_prefix.$fastq_short.".adapterTrim.round2.fastq.gz";  
	    $_ = $fastq_short.".adapterTrim.round2.fastq.gz";
	    my $fastq_path = $repeat_working_directory.$_;
	    
	    my @fastq_fi1_split = split(/\//,$fastq1);
	    my $fastq_fi1_short = $fastq_fi1_split[$#fastq_fi1_split];
	    
	}
    }

    my $merged_bam_file_standardanalysis = $current_file_location_prefix.$dataset_label.".merged.bam";
    my $clip_parsedfi = $repeat_working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam.parsed";
    
    unless (-e $clip_parsedfi) {
	print STDERR "doesn't exists $clip_parsedfi\n";
	$missing_dataset{$dataset_label} = 1;
	next;
    }

    if (exists $alreadydone{$dataset_label}) {
	next;
    }
    $alreadydone{$dataset_label} = 1;

   print STDERR "doing $clip_parsedfi\n";

    &read_parsed_output($clip_parsedfi,$dataset_label);

}

for my $k (keys %counts) {
    for my $fi (keys %clip2input) {
        next unless (exists $submitted_hash{$fi});
        print "$k\t$fi";
	
	$counts{$k}{$fi} = "0\t0" unless (exists $counts{$k}{$fi});
        if (exists $counts{$k}{$fi}) {
            print "\t$counts{$k}{$fi}";

            if (exists $clip2input{$fi}) {
                if (exists $counts{$k}{$clip2input{$fi}}) {
                    my ($clip_n,$clip_rpm) = split(/\t/,$counts{$k}{$fi});
                    my ($input_n,$input_rpm) = split(/\t/,$counts{$k}{$clip2input{$fi}});
                    unless ($input_rpm > 0) {
                        print STDERR "weird error - fi $fi $clip2input{$fi} $clip_rpm $input_rpm\n";
                    }
                    print "\t".$counts{$k}{$clip2input{$fi}}."\t".($clip_rpm / $input_rpm)."\t".&entropy($clip_n,$input_n,$usable_readnum{$fi},$usable_readnum{$clip2input{$fi}});

                } elsif (exists $pseudocount_val{$clip2input{$fi}}) {
                    my ($clip_n,$clip_rpm) = split(/\t/,$counts{$k}{$fi});
                    my ($input_n,$input_rpm) = split(/\t/,$pseudocount_val{$clip2input{$fi}});
                    print "\t".$pseudocount_val{$clip2input{$fi}}."\t".($clip_rpm / $input_rpm)."\t".&entropy($clip_n,$input_n,$usable_readnum{$fi},$usable_readnum{$clip2input{$fi}});
                } elsif (exists $missing_dataset{$clip2input{$fi}}) {
#                   print STDERR "input dataset is missing $fi $clip2input{$fi}\n";                                          
                    print "\tNo_Input\tNo_Input\tNo_Input";

                } else {
                    print STDERR "input present but no pseudocountval? $fi $clip2input{$fi}\n";
                    print "\tNo_Input\tNo_Input\tNo_Input";
                }
            } else {

                print "\tNo_input\tNo_Input\tNo_Input";
            }
        } else {

            print "\tNot_found";
        }
        print "\n";
    }
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
#            $counts{$tmp[1]}{$label} = $tmp[2]."\t".$tmp[3];
	} elsif ($tmp[0] eq "#DELETED") {
	} elsif ($tmp[0] eq "#READINFO") {
            if ($line =~ /UsableReads\t(\d+)\t.+$/) {
#	    if ($line =~ /Union\t(\d+)$/) {
		$pseudocount_val{$label} = "1\t".(1/$1);
		$usable_readnum{$label} = $1;
	    }

	} elsif ($tmp[0] eq "ELEMENT") {
	    my $individ_type = $tmp[4];
	    if ($individ_type eq $desired_type) {
		$counts{$individ_type}{$label} = $tmp[2]."\t".$tmp[3];
#		$counts{$individ_type}{$label}{line} = $line;
	    }
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


sub entropy {
    my $clip_n = shift;
    my $input_n = shift;
    my $clip_usable = shift;
    my $input_usable = shift;

    if ($clip_n == 0) {
        return(0);
    } else {
        my $entropy = (($clip_n / $clip_usable) * log(($clip_n / $clip_usable)/($input_n / $input_usable)) / log(2));
        return($entropy);
    }
}
