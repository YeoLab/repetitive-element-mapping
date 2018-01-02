use warnings;
use strict;

my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20170304.txt";

my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_GRCh38_v1.txt";

#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
my $repeat_working_directory = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/";

#don't think I actually need this anymore here
#my $current_file_location_prefix = "/projects/ps-yeolab3/encode/analysis/encode_master/";


my %usable_readnum;
my %pseudocount_val;
my %counts;

my %clip2input;

#my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_CLIPperv2_20170130/ENCODE_filelist_20170128.txt";


my %submitted_uid;
#my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/ALLDATASETS_submittedonly.txt"; 
my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/integrated_paper_figures/DatasetListing_FINAL.txt.ENCODEFINAL.txt";

my %uid2repids;
open(IN,$inputpairingfi);
for my $line (<IN>) {
    chomp($line);
    $line =~ s/\r//g;

    my @tmp = split(/\s+/,$line);
    
    if ($tmp[3] eq "LNG" || $tmp[3] eq "TAG" || $tmp[3] eq "SingleRep" || $tmp[3] eq "OneRep") {
#	$clip2input{$tmp[1]} = $tmp[2];
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

my %submitted_hash;
my %missing_dataset;
my %hash;


my %file2datasetlabel;
open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $fastqs = shift(@tmp);
    my $genome = shift(@tmp);
    my $dataset_label = shift(@tmp);
#    next unless ($dataset_label =~ /^KB/);
    push @{$hash{$dataset_label}},$fastqs;
#    print STDERR "pushing dataset label $dataset_label\n" if ($dataset_label =~ /^676/);

    my @fastqss = split(/\;/,$fastqs);
    for my $fastq (@fastqss) {
	my @fastqsplit = split(/\//,$fastq);
	my $fastq_short = $fastqsplit[$#fastqsplit];
	$fastq_short =~ s/\.fastq\.gz//;

	my $bam = $fastq_short.".adapterTrim.round2.rmRep.rmDup.sorted.r2.bam";
	$file2datasetlabel{$bam} = $dataset_label;
#	print "bam $bam\n";
    }
}
close(M);

&read_manifest_fi($manifest_fi);




my %alreadydone;
for my $uid (keys %submitted_uid) {
    for my $rep ("01","02","INPUT") {
	my $dataset_label = $uid2repids{$uid}{$rep};
#for my $dataset_label (keys %hash) {
#    next unless (exists $submitted_hash{$dataset_label});
#    next unless ($dataset_label =~ /^676/);

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
}
#exit;
print "element\tfile\tclip_count\tclip_rpr\tinput_count\tinput_rpr\tfold-enrichment\tinformation content\n";
for my $k (keys %counts) {
#    for my $fi (keys %clip2input) {
#	next unless (exists $submitted_hash{$fi});
    for my $uid (keys %submitted_uid) {
	for my $rep ("01","02") {
	    my $fi = $uid2repids{$uid}{$rep};
	    print "$k\t$fi\t".$submitted_uid{$uid}."_".$rep;
	    
	    $counts{$k}{$fi} = "0\t0" unless (exists $counts{$k}{$fi});
	    if (exists $counts{$k}{$fi}) {
		print "\t$counts{$k}{$fi}";
		
#	    print "fi $fi input $clip2input{$fi}\n";
		
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
#		    print STDERR "input dataset is missing $fi $clip2input{$fi}\n";
			print "\tNo_Input\tNo_Input\tNo_Input";
			
		    } else {
			print STDERR "input present but no pseudocountval? $k $fi $clip2input{$fi}\n";
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
            $counts{$tmp[1]}{$label} = $tmp[2]."\t".$tmp[3];
	} elsif ($tmp[0] eq "#DELETED") {
	} elsif ($tmp[0] eq "#READINFO") {
	    if ($line =~ /UsableReads\t(\d+)\t.+$/) {
#		$pseudocount_val{$label} = "1\t".(1000000/$1);
		$pseudocount_val{$label} = "1\t".(1/$1);
#		print STDERR "getting pseudocount val $label\n";
		$usable_readnum{$label} = $1;
#	    } else {
#		print STDERR "couldn't find Union in $fi $label\n";
#		print "reading file $fi $label usable $1\n";
	    }
#        } elsif ($tmp[1] == 1) {
#	    $pseudocount_val{$label} = "1\t".$tmp[2];
#	    print STDERR "label $label $tmp[2]\n";
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

