use warnings;
use strict;
use POSIX;

my %revstrand;
$revstrand{"+"} = "-";
$revstrand{"-"} = "+";
my $scores = "\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/0123456789:;<=>?\@ABCDEFGHIJ";
my %convert_phred;
my @split_scores = split(//,$scores);
for (my $i=0;$i<@split_scores;$i++) {
    $convert_phred{$split_scores[$i]} = $i;
}

my $hashing_value = 1000;
my %enst2gene;
my %convert_enst2type;
my %peaks;
my $repmask_bed_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/RepeatMask.bed";
&read_peakfi($repmask_bed_fi);
my $path = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/";
#my $filelist_file = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/filelist_POLIII";
#my $filelist_file = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.enst2id.fixed";
my $filelist_file = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat";
&read_in_filelists($filelist_file);

my $filelist_file2 = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/ALLRepBase_elements.id_table.FULL";
&read_in_filelists($filelist_file2);

my ($all_count,$duplicate_count,$unique_count,$unique_genomic_count,$unique_repfamily_count) = (0,0,0,0,0);

# takes in sam file, should be blah.bc.tmp file
my $repfamily_sam = $ARGV[0];
my $gabe_rmRep_sam = $ARGV[1];



my %read_hash;
&read_rep_family($repfamily_sam);
&read_unique_mapped($gabe_rmRep_sam);

my $output_fi = $repfamily_sam.".combined_w_uniquemap.rmDup.sam";
open(OUT,">$output_fi");

&run_pcr_duplicate_removal();

my $count_out = $output_fi.".parsed";
open(COUNT,">$count_out");
print COUNT "#READINFO\tAll reads:\t$all_count\tPCR duplicates removed:\t$duplicate_count\tUsable Remaining:\t$unique_count\tUsable from genomic mapping:\t$unique_genomic_count\tUsable from family mapping:\t$unique_repfamily_count\n";
#print COUNT "#READINFO\tAll reads:\t$all_count\nPCR duplicates removed:\t$duplicate_count\nUsable Remainig:\t$unique_count\nUsable from genomic mapping:\t$unique_genomic_count\nUsable from family mapping:\t$unique_repfamily_count\n";

&count_output();
close(COUNT);

my $out_done = $count_out.".done";
open(DONE,">$out_done");
print DONE "jobs done\n";
close(DONE);


sub count_output {
#    print STDERR "doing unique read counting\n";
    my %counter;
    for my $r1name (keys %read_hash) {
	if (exists $read_hash{$r1name}{file1flag} && $read_hash{$r1name}{file1flag} == 1) {
	    if (exists $read_hash{$r1name}{file2flag} && $read_hash{$r1name}{file2flag} == 1) {
		$counter{"intersect"}++;
	    } else {
		$counter{"bam1only"}++;
	    }
	} elsif (exists $read_hash{$r1name}{file2flag} && $read_hash{$r1name}{file2flag} == 1) {
	    $counter{"bam2only"}++;
	} else {
	    print STDERR "fatal error - should never hit this $r1name\n";
	}
	$counter{"union"}++;
    }

    for my $key ("intersect","union","bam1only","bam2only") {
	$counter{$key} = 0 unless (exists $counter{$key});
    }
    my ($a_only,$b_only,$intersect,$union) = ($counter{"bam1only"},$counter{"bam2only"},$counter{"intersect"},$counter{"union"});
    my $total_unique_mapped_read_num = $unique_count;
    print COUNT "#READINFO2\tRepMappingOnly\t$a_only\tUniqueMappingOnly\t$b_only\tIntersection\t$intersect\tUnion\t$union\n";


#    print STDERR "now doing type counting and final stuff\n";
    my %count;
    my %count_enst;
    for my $read (keys %read_hash) {
	next unless (exists $read_hash{$read}{rep_flag});
	my $ensttype = $read_hash{$read}{ensttype};

	if ($ensttype =~ /Simple\_repeat/) {
	    $ensttype = "Simple_repeat";
	}

	$count{$ensttype}++;

	$count_enst{$ensttype."||".$read_hash{$read}{mult_ensts}}++;
#           $count_enst{$ensttype."|".join("|",@{$read_hash{$read}{mult_ensts}{$ensttype}})}++;                  
    }


    for my $k (keys %count) {
	print COUNT "TOTAL\t$k\t$count{$k}\t".sprintf("%.5f",$count{$k} * 1000000 / $total_unique_mapped_read_num)."\n";
    }

    my @sorted = sort {$count_enst{$b} <=> $count_enst{$a}} keys %count_enst;
    for my $s (@sorted) {
	my ($ensttype,$multensts) = split(/\|\|/,$s);
	my @gids = split(/\|/,$multensts);
	my $type = $ensttype;
	my @genes;
	for my $gid (@gids) {
	    if (exists $enst2gene{$gid}) {
		push @genes,$enst2gene{$gid};
	    } else {
		push @genes,$gid;
	    }
	}
	
	print COUNT "$type\t$count_enst{$s}\t".sprintf("%.5f",$count_enst{$s} * 1000000 / $total_unique_mapped_read_num)."\t$s\t".join("|",@genes)."\n"; 
#    print "$type\t$geneid\t$count_enst{$s}\t$s\n";                      
    }



}

# now do PCR duplicate removal

sub run_pcr_duplicate_removal {
    my %fragment_hash;
    for my $r1name (keys %read_hash) {
	my $r1 = $read_hash{$r1name}{R1};
	my $r2 = $read_hash{$r1name}{R2};
	
	my @tmp_r1 = split(/\t/,$r1);
	my @tmp_r2 = split(/\t/,$r2);
	my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
	my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);
	
	my $r1sam_flag = $tmp_r1[1];
	my $r2sam_flag = $tmp_r2[1];
	unless ($r1sam_flag) {
	    print STDERR "error $r1 $r2\n";
	}
	next if ($r1sam_flag == 77 || $r1sam_flag == 141);
	
	my $frag_strand;
### This section is for only properly paired reads                                                                                                                                                                       
	if ($r1sam_flag == 99 || $r1sam_flag == 355) {
	    $frag_strand = "-";
	} elsif ($r1sam_flag == 83 || $r1sam_flag == 339) {
	    $frag_strand = "+";
	} elsif ($r1sam_flag == 147 || $r1sam_flag == 403) {
	    $frag_strand = "-";
	} elsif ($r1sam_flag == 163 || $r1sam_flag == 419) {
	    $frag_strand = "+";
	}  else {
	    next;
	    print STDERR "R1 strand error $r1sam_flag\n";
	}

	
	my @read_name = split(/\:/,$tmp_r1[0]);
	my $randommer = $read_name[0];

	my $r1_cigar = $tmp_r1[5];
	my $r2_cigar = $tmp_r2[5];

      # 165 = R2 unmapped, R1 rev strand -- frag on fwd strand


	my $r1_chr = $tmp_r1[2];
	my $r1_start = $tmp_r1[3];

	my $r2_chr = $tmp_r2[2];
	my $r2_start = $tmp_r2[3];


	my $mismatch_flags = $tmp_r1[14];

	my %read_regions;
	@{$read_regions{"R1"}} = &parse_cigar_string($r1_start,$r1_cigar,$r1_chr,$frag_strand);
	@{$read_regions{"R2"}} = &parse_cigar_string($r2_start,$r2_cigar,$r2_chr,$frag_strand);

	my $hashing_value;
	if ($frag_strand eq "+") {
	    my $r1_firstregion = $read_regions{"R1"}[scalar(@{$read_regions{"R1"}})-1];
	    my ($r1_firstchr,$r1_firststr,$r1_firstpos) = split(/\:/,$r1_firstregion);
	    my ($r1_firststart,$r1_firststop) = split(/\-/,$r1_firstpos);
	    my $r2_firstregion = $read_regions{"R2"}[0];
	    my ($r2_firstchr,$r2_firststr,$r2_firstpos) = split(/\:/,$r2_firstregion);
	    my ($r2_firststart,$r2_firststop) = split(/\-/,$r2_firstpos);
	    $hashing_value = $r2_firstchr.":".$r2_firststr.":".$r2_firststart."\t".$r1_firstchr.":".$r1_firststr.":".$r1_firststop;
	} elsif ($frag_strand eq "-") {
	    my $r1_firstregion = $read_regions{"R1"}[0];
	    my ($r1_firstchr,$r1_firststr,$r1_firstpos) = split(/\:/,$r1_firstregion);
	    my ($r1_firststart,$r1_firststop) = split(/\-/,$r1_firstpos);
	    my $r2_firstregion = $read_regions{"R2"}[scalar(@{$read_regions{"R2"}})-1];
	    my ($r2_firstchr,$r2_firststr,$r2_firstpos) = split(/\:/,$r2_firstregion);
	    my ($r2_firststart,$r2_firststop) = split(/\-/,$r2_firstpos);
	    $hashing_value = $r1_firstchr.":".$r1_firststr.":".$r1_firststart."\t".$r2_firstchr.":".$r2_firststr.":".$r2_firststop;

	} else {
	    print STDERR "strand error $frag_strand $r1 $r2\n";
	}

	my $full_frag_position = join("|",@{$read_regions{"R1"}})."\t".join("|",@{$read_regions{"R2"}});

	my $debug_flag = 0;


	my %full_fragment;
	$all_count++;
	if (exists $fragment_hash{$r1_chr."|".$frag_strand}{$hashing_value.":".$randommer}) {
	    $duplicate_count++;
	    delete($read_hash{$r1name});
	} else {
	    $fragment_hash{$r1_chr."|".$frag_strand}{$hashing_value.":".$randommer} = 1;

	    if (exists $read_hash{$r1name}{file1flag} && $read_hash{$r1name}{file1flag} == 1) {
		# read is coming from familiy mapping
		print OUT "".$r1."\tRepFamily\t$r1_chr\n".$r2."\tRepFamily\t$r1_chr\n";
		$unique_repfamily_count++;
	    } elsif (exists $read_hash{$r1name}{file2flag} && $read_hash{$r1name}{file2flag} == 1) {
		# read is coming from familiy mapping                                              
                print OUT "".$r1."\tUniqueGenomic\t".$read_hash{$r1name}{ensttype}."||".$read_hash{$r1name}{mult_ensts}."\n".$r2."\tUniqueGenomic\t".$read_hash{$r1name}{ensttype}."||".$read_hash{$r1name}{mult_ensts}."\n";
		$unique_genomic_count++;
	    } else {
		print STDERR "this shouldn't be hit - fi1flag $read_hash{$r1name}{file1flag} fi2flag $read_hash{$r1name}{file2flag} $r1name\n";
	    }


	    $unique_count++;
	}
    }
}








sub read_unique_mapped {
    my %deleted;
    my $fi2_count=0;
    my %read_hash_bam2;
    my $gabe_rmRep_sam_file = shift;
#    print STDERR "opening $gabe_rmRep_sam_file\n";
    if ($gabe_rmRep_sam_file =~ /\.sam$/ || $gabe_rmRep_sam_file =~ /\.tmp$/) {
	open(B,$gabe_rmRep_sam_file)  || die "no $gabe_rmRep_sam_file\n";
    } elsif ($gabe_rmRep_sam_file =~/\.bam$/) {
	open(B,"samtools view -h $gabe_rmRep_sam_file |") || die "no $gabe_rmRep_sam_file\n";
    } else {
	die "couldn't figure out format of $gabe_rmRep_sam_file\n";
    }
    while (<B>) {
	my $r1 = $_;
	next if ($r1 =~ /^\@/);
	$fi2_count++;
	print STDERR "read $fi2_count\n" if ($fi2_count % 100000 == 0);
	my $r2 = <B>;
	chomp($r1);
	chomp($r2);
	
	my @tmp_r1 = split(/\t/,$r1);
	my @tmp_r2 = split(/\t/,$r2);
	my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
	my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);
	
	unless ($tmp_r1[0] eq $tmp_r2[0]) {
	    print STDERR "paired end mismatch error: $tmp_r1[0] $tmp_r2[0]\n";
	}
	
	my $r1sam_flag = $tmp_r1[1];
	my $r2sam_flag = $tmp_r2[1];
	next if ($r1sam_flag == 77 || $r1sam_flag == 141);

        # 77 = R1, unmapped                      
        # 141 = R2, unmapped                     

        # 99 = R1, mapped, fwd strand --- frag on rev strand             
        # 101 = R1 unmapped, R2 mapped rev strand -- frag on rev strand  
        # 73 = R1, mapped, fwd strand --- frag on rev strand             
        # 147 = R2, mapped, rev strand -- frag on rev strand             
        # 153 = R2 mapped (R1 unmapped), rev strand -- frag on rev strand
        # 133 = R2 unmapped, R1 mapped fwd strand -- frag on rev strand  

        # 83 = R1, mapped, rev strand --- frag on fwd strand             
        # 69 = R1 unmapped, R2 mapped fwd strand -- frag on fwd strand   
        # 89 = R1 mapped rev strand, R2 unmapped -- frag on fwd strand   
        # 163 = R2, mapped, fwd strand -- frag on fwd strand             
        # 137 = R2 mapped (R1 unmapped), fwd strand -- frag on fwd strand
        # 165 = R2 unmapped, R1 rev strand -- frag on fwd strand         


	my $frag_strand;
### This section is for only properly paired reads                       
	if ($r1sam_flag == 99) {
	    $frag_strand = "-";
	} elsif ($r1sam_flag == 83) {
	    $frag_strand = "+";
	} elsif ($r1sam_flag == 147) {
	    $frag_strand = "-";
	    ($r1,$r2) = ($r2,$r1);
	    @tmp_r1 = split(/\t/,$r1);
	    @tmp_r2 = split(/\t/,$r2);	    
	    
	} elsif ($r1sam_flag == 163) {
	    ($r1,$r2) =($r2,$r1);
	    @tmp_r1 = split(/\t/,$r1);
	    @tmp_r2 = split(/\t/,$r2);
	    $frag_strand = "+";
	}  else {
#       print STDERR "R1 strand error $r1sam_flag\n";
	    next;
	}
###
	my @read_name = split(/\:/,$tmp_r1[0]);
	my $randommer = $read_name[0];
	
	my $r1_cigar = $tmp_r1[5];
	my $r2_cigar = $tmp_r2[5];
	
	my $r1_chr = $tmp_r1[2];
	my $r1_start = $tmp_r1[3];
	
	my $r2_chr = $tmp_r2[2];
	my $r2_start = $tmp_r2[3];
	
	my @read_regions = &parse_cigar_string($r2_start,$r2_cigar,$r2_chr,$frag_strand);

	my %tmp_hash;
	for my $region (@read_regions) {
	    my ($rchr,$rstr,$rpos) = split(/\:/,$region);
	    my ($rstart,$rstop) = split(/\-/,$rpos);

	    my $verbose_flag = 0;

	    my $rx = int($rstart / $hashing_value);
	    my $ry = int($rstop  / $hashing_value);
	    for my $ri ($rx..$ry) {

		for my $peak (@{$peaks{$rchr}{$rstr}{$ri}}) {
		    my ($pchr,$ppos,$pstr,$ptype) = split(/\:/,$peak);
		    my ($pstart,$pstop) = split(/\-/,$ppos);

#                   next if ($pstart > $rstop || $pstop < $rstart);      
#                   if ($pstop == $rstart) {     
#                   if ($pstart == $rstop) {     
#                   if ($pstart == $rstart) {    
#                       print "region $region $r1 peak $peak\n";         
#                   }    
		    next if ($pstart >= $rstop || $pstop <= $rstart);
		    $tmp_hash{$peak} = 1;

		}
	    }
	}

	my %temp_peak_read_counts;
	my %converted_types;
	for my $peak (keys %tmp_hash) {
	    my ($pchr,$ppos,$pstr,$ptype) = split(/\:/,$peak);
	    $temp_peak_read_counts{$ptype."_uniquegenomic"}++;
	    $converted_types{$convert_enst2type{$ptype}} = 1;
	    print STDERR "peak $peak $ptype\n" unless (exists $convert_enst2type{$ptype});
	}



	my $flags_r1 = join("\t",@tmp_r1[11..$#tmp_r1]);
	my $flags_r2 = join("\t",@tmp_r2[11..$#tmp_r2]);

	my $r1_mismatch;
	my $r2_mismatch;
	if ($flags_r1 =~ /MD\:Z\:(\S+)\s/ || $flags_r1 =~ /MD\:Z\:(\S+?)$/) {                       
	    $r1_mismatch = $1;          
	} 
	if ($flags_r2 =~ /MD\:Z\:(\S+)\s/ || $flags_r2 =~ /MD\:Z\:(\S+?)$/) {        
	    $r2_mismatch = $1;
	} 
          
	my $r1_phred = $tmp_r1[10];                   
	my $r2_phred = $tmp_r2[10];          
	my $r1_seq = $tmp_r1[9];                      
	my $r2_seq = $tmp_r2[9];       
	my $r1_mmscore = &get_alignment_score($r1_mismatch,$r1_cigar,$r1_phred,$r1_seq);               
	my $r2_mmscore = &get_alignment_score($r2_mismatch,$r2_cigar,$r2_phred,$r2_seq);           
	my $total_mmscore = $r1_mmscore+$r2_mmscore;         
	
#check if read already exists from family mapping - 8/31/16 - this runs REGARDLESS of whether unique mapping overlaps a RepBase element or not	
	if (exists $read_hash{$r1name}{file1flag}) {
	    unless (exists $read_hash{$r1name}{rep_score}) {
		print STDERR "err $r1name\n";
	    }
	    if ($total_mmscore > $read_hash{$r1name}{rep_score} + 24) {
            
            # if unique mapping to genome is more than 2 mismatches per read = 2 * 2 * 6 alignment score better than to repeat element, throw out repeat element and use genome mapping
#		print "deleting $r1name from rep element mapping - $total_mmscore $read_hash{$r1name}{rep_score}\n$r1\n$read_hash{$r1name}{R1}\n----\n";
		$deleted{$read_hash{$r1name}{ensttype}}++;
		delete($read_hash{$r1name}) if (exists $read_hash{$r1name});
		
	    } else {
		# repeat mapping exists and is better than unique genome mapping - throw out unique genome mapping
		next;
	    }
	}

	my @sorted_types = sort {$a cmp $b} keys %temp_peak_read_counts;
	my $all_mapped_ensts = join("|",@sorted_types);
	my @sorted_convertedtype = sort {$a cmp $b} keys %converted_types;
	my $all_ensttypes = join("|",@sorted_convertedtype);

	# now throw out reads that uniquely map but don't overlap a RepBase element - removed for now, can add back in later
	unless (scalar(keys %temp_peak_read_counts) >= 1) {
	    $all_ensttypes = "unique_genomic";
	    $all_mapped_ensts = "unique_genomic";
	}
#	next unless (scalar(keys %temp_peak_read_counts) >= 1);



#	print "keeping new mapping $all_mapped_ensts \t $all_ensttypes\n$r1\n";


	$read_hash{$r1name}{rep_flag} = 0;
        $read_hash{$r1name}{rep_score} = $total_mmscore;
        $read_hash{$r1name}{file2flag} = 1;
        $read_hash{$r1name}{R1} = $r1;
        $read_hash{$r1name}{R2} = $r2;
        $read_hash{$r1name}{mult_ensts} = $all_mapped_ensts;
        $read_hash{$r1name}{ensttype} = $all_ensttypes;

    }
    close(B);
}


sub read_rep_family {
    
    my $sam_file_familymapping = shift;
#    print STDERR "opening $sam_file_familymapping\n";
    open(BAM,$sam_file_familymapping) || die "can't open $sam_file_familymapping\n";
#open(BAM,"-|", "samtools view $rmDup_bam");                                                                                                                                    
                                                                 
    while (<BAM>) {
	my $r1 = $_;
	chomp($r1);
	if ($r1 =~ /^\@/) {
	    next;
	}
	
	my $r2 = <BAM>;
	chomp($r2);
	
	my @tmp_r1 = split(/\t/,$r1);
	my @tmp_r2 = split(/\t/,$r2);
	
	my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
	my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);
	
	unless ($r1name eq $r2name) {
	    print STDERR "paired end mismatch error: $sam_file_familymapping r1 $tmp_r1[0] r2 $tmp_r2[0]\n";
	}
	
	my $r1sam_flag = $tmp_r1[1];
	my $r2sam_flag = $tmp_r2[1];
	unless ($r1sam_flag) {
	    print STDERR "error $r1 $r2\n";
	}
	next if ($r1sam_flag == 77 || $r1sam_flag == 141);
	
	my $frag_strand;
### This section is for only properly paired reads   
	if ($r1sam_flag == 99 || $r1sam_flag == 355) {
	    $frag_strand = "-";
	} elsif ($r1sam_flag == 83 || $r1sam_flag == 339) {
	    $frag_strand = "+";
	} elsif ($r1sam_flag == 147 || $r1sam_flag == 403) {
            ($r1,$r2) =($r2,$r1);
	    @tmp_r1 = split(/\t/,$r1);
	    @tmp_r2 = split(/\t/,$r2);
	    $frag_strand = "-";
	} elsif ($r1sam_flag == 163 || $r1sam_flag == 419) {
            ($r1,$r2) =($r2,$r1);
	    @tmp_r1 = split(/\t/,$r1);
	    @tmp_r2 = split(/\t/,$r2);
	    $frag_strand = "+";
	}  else {
	    next;
	    print STDERR "R1 strand error $r1sam_flag\n";
	}
###                                                                                                                                                                                                                                              
	my $flags_r1 = join("\t",@tmp_r1[11..$#tmp_r1]);
	my $flags_r2 = join("\t",@tmp_r2[11..$#tmp_r2]);
	
	my $r1_score;
	if ($flags_r1 =~ /AS\:i\:(\S+)\s/ || $flags_r1 =~ /AS\:i\:(\S+)$/) {
	    $r1_score = $1;
	} else {
	    print STDERR "couldn't find score for r1 $r1\n";
	}
	my $r2_score;
	if ($flags_r2 =~ /AS\:i\:(\S+)\s/ || $flags_r2 =~ /AS\:i\:(\S+)$/) {
	    $r2_score = $1;
	} else {
	    print STDERR "couldn't find score for r2 $r2\n";
	}
	my $total_score = $r1_score + $r2_score;
	
	my $all_mapped_ensts;
	if ($flags_r1 =~ /ZZ\:Z\:(\S+?)\s/ || $flags_r1 =~ /ZZ\:Z\:(\S+?)$/) {
	    $all_mapped_ensts = $1;
	} else {
	    print STDERR "didn't match? $flags_r1\n";
	}
	
        # if read has never been seen before, keep first mapping                                                                                                                                                                                 
	my $all_info = $tmp_r1[2];
	my ($all_ensttypes,$all_primary_enst) = split(/\|\|/,$all_info);
#    my $all_primary_enst = $tmp_r1[2];                                                                                                                                                                                                          
	my @ensttypes = split(/\|/,$all_ensttypes);
	my @mapped_ensts = split(/\|/,$all_primary_enst);
	$read_hash{$r1name}{rep_flag} = 1;
	$read_hash{$r1name}{rep_score} = $total_score;
	$read_hash{$r1name}{file1flag} = 1;
	$read_hash{$r1name}{R1} = $r1;
	$read_hash{$r1name}{R2} = $r2;
#    $read_hash{$r1name}{flags} = $enstpriority;                                                                                                                                                                                                 
	$read_hash{$r1name}{mult_ensts} = $all_mapped_ensts;
#    print "al $all_mapped_ensts\n";                                                                                                                                                                                                             
	$read_hash{$r1name}{ensttype} = $all_ensttypes;
    }
    close(BAM);

}



sub read_peakfi {
    my $fi = shift;
    my $duplicate = 0;
    my $lc=0;
#    print STDERR "reading peak file $fi\n";
    open(F,$fi) || die "no peak file $fi\n";
    while (<F>) {
        my $line = $_;
        chomp($line);
        my @tmp = split(/\t/,$line);
        my $chr = shift(@tmp);
        my $start = shift(@tmp);
        my $stop = shift(@tmp);
        my $gene = shift(@tmp);
        my $pval = shift(@tmp);
        my $strand = shift(@tmp);

        $gene =~ s/\_$//g;
        
        next if ($chr eq "genoName");

        my $x = int($start / $hashing_value);
        my $y = int($stop  / $hashing_value);

        my $peak = $chr.":".$start."-".$stop.":".$strand.":".$gene;
#       if (exists $peak_read_counts{$peak}{allpeaks}) {
#           $duplicate++;
#       } else {
#           $lc++;
#       }
#        $peak_read_counts{$gene}{allpeaks} = 1;
#        $peak_read_counts{$gene."_antisense"}{allpeaks} = 1;
        
        
        for my $i ($x..$y) {
            push @{$peaks{$chr}{$strand}{$i}},$chr.":".$start."-".$stop.":".$strand.":".$gene;
            push @{$peaks{$chr}{$revstrand{$strand}}{$i}},$chr.":".$start."-".$stop.":".$strand.":antisense_".$gene;
        }
    }
    close(F);
#    print STDERR "duplciate $duplicate\nlc $lc\n";
}


sub read_in_filelists {
    my $fi = shift;
    my $priority_n = 0;
    my %corrected_ensts;
    open(F,$fi);
    for my $line (<F>) {
        chomp($line);
        my ($allenst,$allensg,$gid,$type_label,$typefile) = split(/\t/,$line);
#       my ($allensg,$gid,$allenst) = split(/\t/,$line);
        unless ($allenst) {
            print STDERR "error missing enst $line $fi\n";
        }
        my @ensts = split(/\|/,$allenst);
        $gid =~ s/\?$//;
        $gid =~ s/\_$//;
        $type_label =~ s/\?$//;
        for my $enst (@ensts) {
	    if ($enst =~ /\_$/) {
		$corrected_ensts{$enst}++;
		$enst =~ s/\_$//;
	    }
            $enst2gene{$enst} = $gid;
            $enst2gene{"antisense_".$enst} = "antisense_".$gid;
#            $convert_enst2type{$enst} = $type_label.":".$priority_n;
#            $convert_enst2type{$enst."_antisense"} = $type_label."_antisense:".$priority_n;
            $convert_enst2type{$enst} = $type_label;
            $convert_enst2type{"antisense_".$enst} = "antisense_".$type_label;
            $priority_n++;
        }
    }
    close(F);
#    for my $k (keys %corrected_ensts) {
#	print "$k\t$corrected_ensts{$k}\n";
#    }
}

sub min {
    my $x = shift;
    my $y = shift;
    if ($x < $y) {
        return($x);
    } else {
        return($y);
    }
}



sub get_alignment_score {
    my $md_string = shift;
    my $cigar_string = shift;
    my $quality_string = shift;
    my $read_seq = shift;
    my (@insertions) = &parse_cigar_string_foralignmentscore($cigar_string,$quality_string);
    my $gap_score = shift(@insertions);
    my $mm_score = &parse_mismatch_string_foralignmentscore($md_string,$quality_string,\@insertions,$read_seq);
    my $total_score = 0 - $mm_score - $gap_score;
    return($total_score);
}




sub parse_cigar_string {
    my $region_start_pos = shift;
    my $flags = shift;
    my $chr = shift;
    my $strand = shift;

    my $current_pos = $region_start_pos;
    my @regions;

    while ($flags =~ /(\d+)([A-Z])/g) {
       
        if ($2 eq "N") {
            #read has intron of N bases at location

# 1 based, closed ended to 0 based, right side open ended fix            
#            push @regions,$chr.":".$strand.":".$region_start_pos."-".($current_pos-1);
            push @regions,$chr.":".$strand.":".($region_start_pos-1)."-".($current_pos-1);

            $current_pos += $1;
            $region_start_pos = $current_pos;
        } elsif ($2 eq "M") {
            #read and genome match
            $current_pos += $1;
        } elsif ($2 eq "S") {
            #beginning of read is soft-clipped; mapped pos is actually start pos of mapping not start of read
        } elsif ($2 eq "I") {
            #read has insertion relative to genome; doesn't change genome position
        } elsif ($2 eq "D") {
#           push @read_regions,$chr.":".$current_pos."-".($current_pos+=$1);
            $current_pos += $1;
            #read has deletion relative to genome; genome position has to increase
        } else {
            print STDERR "flag $1 $2 $flags\n";

        }
    }

# 1 based, closed ended to 0 based, right side open ended fix            
# $region_start_pos is 1-based, closed ended -> ($region_start_pos-1) is 0-based, closed ended
# ($current_pos-1) is 1-based, closed ended -> ($current_pos-1-1) is 0-based, closed ended -> ($current_pos-1-1+1) is 0-based, open ended

#    push @regions,$chr.":".$strand.":".$region_start_pos."-".($current_pos-1);
    push @regions,$chr.":".$strand.":".($region_start_pos-1)."-".($current_pos-1);

    return(@regions);
}





sub parse_mismatch_string_foralignmentscore {
    my $flags = shift;
    my $phred_scores = shift;
    my $insertionsref = shift;
    my $read_seq = shift;

    my @insertions = @$insertionsref;
    my ($next_insertion_pos,$next_insertion_len);
    my $insertion_flag = 0;
    if (@insertions > 0) {
        if ($insertions[0] eq "none") {
            $insertion_flag = 0;
        } else {
            $insertion_flag = 1;
            my $next_insertion = shift(@insertions);
            ($next_insertion_pos,$next_insertion_len) = split(/\|/,$next_insertion);
        }
    } else {
        $insertion_flag = 0;
    }

    if ($flags =~ /^\d+$/) {
        return(0);
    }
    my @phred_scoress = split(//,$phred_scores);
    my @quality_scores;
    for (my $i=0;$i<@phred_scoress;$i++) {
        $quality_scores[$i] = $convert_phred{$phred_scoress[$i]};
    }

    my $current_pos = 0;

    my $mm_penalty = 0;
    my @regions;
    my $mn = 2;
    my $mx = 6;

    my $current_read_pos = 0;
    while ($flags) {
        if ($flags =~ /^(\d+)/) {
            $current_read_pos += $1;
            if ($insertion_flag == 1) {
                while ($insertion_flag == 1 && $current_read_pos >= $next_insertion_pos) {
                    $current_read_pos += $next_insertion_len;
                    if (@insertions > 0) {
                        $insertion_flag = 1;
                        my $next_insertion = shift(@insertions);
                        ($next_insertion_pos,$next_insertion_len) = split(/\|/,$next_insertion);
                    } else {
                        $insertion_flag = 0;
                    }
                }
            }

            $flags = substr($flags,length($1));
        } elsif ($flags =~ /^([A-Z]+)/) {
            my $len = length($1);

            for my $j (1..$len) {
                my $base_score = $quality_scores[$current_read_pos+$j-1];
                unless ($base_score) {
                    print STDERR "$j $j cur $current_read_pos read $read_seq $phred_scores\n";
                }
                my $base = substr($1,$j-1,1);
                my $base_seq = substr($read_seq,$current_read_pos+$j-1,1);
                my $base_mm_penalty;
                if ($base eq "N" || $base_seq eq "N") {
                    $base_mm_penalty = 1;
                } elsif (($base eq "A" || $base eq "C" || $base eq "G" || $base eq "T") && ($base_seq eq "A" || $base_seq eq "C" || $base_seq eq "G" || $base_seq eq "T")) {
                    $base_mm_penalty = $mn + floor( ($mx-$mn) * (&min($base_score,40)/40) );
                } else {
                    print STDERR "unexpected base - $base\n";
                }
                $mm_penalty += $base_mm_penalty;
            }
            $current_read_pos += $len;
            $flags = substr($flags,length($1));
        } elsif ($flags =~ /^\^([A-Z]+)\d/) {
            my $len = length($1);
            #insertions aren't penalized here; penalized in gap open/close                    
            $flags = substr($flags,$len+1);
        } else {
            print STDERR "this is a flag I'm not expecting $flags\n";
        }
    }
    return($mm_penalty);
}

sub parse_cigar_string_foralignmentscore {
    my $flags = shift;
    my $phred_scores = shift;

    if ($flags =~ /^\d+M$/) {
        return(0);
    }
    my @phred_scoress = split(//,$phred_scores);
    my @quality_scores;
    for (my $i=0;$i<@phred_scoress;$i++) {
        $quality_scores[$i] = $convert_phred{$phred_scoress[$i]};
    }

    my $gap_open = 5;
    my $gap_extend = 3;
    my $current_pos = 0;
    my $mm_penalty = 0;
    my @insertions;
    my $current_read_pos = 0;

    while ($flags =~ /(\d+)([A-Z])/g) {
        if ($2 eq "N") {
            # intron - not sure what to do with this for now? I think skip

            $current_read_pos += $1;
            $current_pos += $1;
        } elsif ($2 eq "M") {
            $current_read_pos += $1;
            $current_pos += $1;
        } elsif ($2 eq "S") {
            $current_read_pos += $1;
            $mm_penalty += $gap_open + $1 * $gap_extend;
        } elsif ($2 eq "I") {
            push @insertions,$current_read_pos."|".$1;
            $current_read_pos += $1;
            $mm_penalty += $gap_open + $1 * $gap_extend;
        } elsif ($2 eq "D") {
            $current_pos += $1;
            $mm_penalty += $gap_open + $1 * $gap_extend;
        } else {
            print STDERR "flag $1 $2 $flags\n";

        }
    }
    return($mm_penalty,@insertions);

}
