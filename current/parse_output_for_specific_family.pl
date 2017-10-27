use warnings;
use strict;

my $desired_type = "RNY5";
#my $desired_type = "L1||L1HS";
#my $desired_type = "antisense_L1||antisense_L1HS";

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
	$clip2input{$tmp[2]} = $tmp[3];
    }

}
close(IN);

my $manifestfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS.txt";
#my $manifestfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ALLCLIP_v12_20160112/encode_v12_filelist.allencode_20160314.txt.fixedgenename.txt";

#my $repeat_working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/ALLCLIP_v12_20160112/repeatmapping/";
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_v12_repeatmapping/";
#my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_v12_20160506repeatmapping/";
my $repeat_working_directory = "/home/elvannostrand/scratch/ALLCLIP_repeatfamilymapping_20160723b/";

my %missing_dataset;

my $counter=0;
my %hash;
#my $current_file_location_prefix = "/projects/ps-yeolab2/encode/analysis/encode_master/";
my $current_file_location_prefix = "/projects/ps-yeolab3/encode/analysis/encode_master/";
my $merging_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/encode_GRCh38_v1.txt";
#my $merging_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/test_processing_GRCh38.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_master.txt";
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
my %alreadydone;
for my $dataset_label (keys %hash) {
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
    my $clip_parsedfi = $repeat_working_directory.$dataset_label."_RepMappermerged.rmDup.sam.parsed";
    
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
#    if (exists $clip2input{$dataset_label}) {
#	print STDERR "doing $clip_parsedfi paired with $input_parsedfi\n";	
#	my $input_parsedfi = $repeat_working_directory.$clip2input{$dataset_label}."_RepMappermerged.rmDup.bam.namesort.bam.parsed";
#	&read_parsed_output($input_parsedfi,$clip2input{$dataset_label});
#   }

    ## this is gabe's rmDup final bam file
    if (scalar(@{$hash{$dataset_label}}) == 1) { 
#       if (uc($dataset_label) =~ /INPUT/) {
	
	my $fastq = $hash{$dataset_label}[0];
	my ($fastq1,$fastq2) = split(/\;/,$fastq);
	my @fastq1_split = split(/\//,$fastq1);
	my $fastq1_short = $fastq1_split[$#fastq1_split];
	my $r1f1 = $fastq1_short;
	$r1f1 =~ s/\.fastq.gz//;
	$r1f1 .= ".adapterTrim.round2.rmRep.rmDup.bam";
	$merged_bam_file_standardanalysis = $current_file_location_prefix.$r1f1;
    }
}


for my $fi (keys %clip2input) {
    print "$fi";
    
    if (exists $counts{$fi}) {
	print "\t".join("\|",@{$counts{$fi}{ensts}})."\t".join("\|",@{$counts{$fi}{types}});
	
	if (exists $clip2input{$fi}) {
	    if (exists $counts{$clip2input{$fi}}) {
		my ($clip_n,$clip_rpm) = ($counts{$fi}{count},$counts{$fi}{rpm});
		my ($input_n,$input_rpm) = ($counts{$clip2input{$fi}}{count},$counts{$clip2input{$fi}}{rpm});
		print "\t".$counts{$clip2input{$fi}}{count}."\t".($clip_rpm / $input_rpm);
	    } elsif (exists $pseudocount_val{$clip2input{$fi}}) {
		my ($clip_n,$clip_rpm) = ($counts{$fi}{count},$counts{$fi}{rpm});
		my ($input_n,$input_rpm) = split(/\t/,$pseudocount_val{$clip2input{$fi}});
		print "\t".$pseudocount_val{$clip2input{$fi}}."\t".($clip_rpm / $input_rpm);
	    } elsif (exists $missing_dataset{$clip2input{$fi}}) {
		print STDERR "input dataset is missing $fi $clip2input{$fi}\n";
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
	    if ($line =~ /Union\t(\d+)$/) {
		$pseudocount_val{$label} = "1\t".(1000000/$1);
	    } else {
		print STDERR "couldn't find Union in $fi $label\n";
	    }
#        } elsif ($tmp[1] == 1) {
#	    $pseudocount_val{$label} = "1\t".$tmp[2];
#	    print STDERR "label $label $tmp[2]\n";
	} else {
#	    my $individ_type = $tmp[3];
	    my @individ_types = split(/\|/,$tmp[4]);
	    for my $individ_type (@individ_types) {
		if ($individ_type =~ /^$desired_type/) {
		    $counts{$label}{count} += $tmp[1];
		    $counts{$label}{rpm} += $tmp[2];
		    push @{$counts{$label}{ensts}},$tmp[3];
		    push @{$counts{$label}{types}},$tmp[4];
		}
	    }
	}
    }
    close(F);
}
