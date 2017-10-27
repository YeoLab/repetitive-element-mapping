use warnings;
use strict;

unless ($ARGV[0] && $ARGV[1]) {
    print STDERR "usage: perl wrapper_FULLversion_includemultifamily_submit.pl working_directory manifest_fi fastq_folder_prefix(optional-defaults to encode)\n";
    print STDERR "note - manifest_fi is gabe's pipeline file for merging fastq files\n";
    print STDERR "\n";
    print STDERR "current files used are:\n";
    print STDERR "0 working directory = /home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/\n";
    print STDERR "1 current manifest = /home/gpratt/projects/encode/scripts/encode_GRCh38_v1.txt\n";
    print STDERR "2 current fastq file prefix = /projects/ps-yeolab3/encode/analysis/encode_master/\n";
    print STDERR "3 subset manifest = /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/ALLDATASETS_submittedonly.txt\n";

    exit;
}

my $keep_temp_files_flag = 0;

my $current_file_location_prefix = "/projects/ps-yeolab3/encode/analysis/encode_master/";
if ($ARGV[2]) {
    $current_file_location_prefix = $ARGV[2];
    $current_file_location_prefix =~ s/\/$//;
    $current_file_location_prefix .= "/";
}
#my $current_file_location_prefix = "/projects/ps-yeolab2/encode/analysis/encode_master/";


my $working_directory = $ARGV[0];
# working directory is currently /home/elvannostrand/scratch/ALLCLIP_v12_repeatmapping/

$working_directory =~ s/\/$//;
$working_directory .= "/";

system("mkdir $working_directory") unless (-d $working_directory);
#my $working_directory = "/home/elvannostrand/scratch/ENCODE_v9_20151209/PolIII_mapping/";
#my $working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_v9_20151209/PolIII_mapping/";


#$ARGV[3] should be a subset list
my %subset_hash;
my $subset_flag = 0;
if (exists $ARGV[3] && $ARGV[3]) {
    my $subset_fi = $ARGV[3];
    open(SUB,$subset_fi);
    for my $line (<SUB>) {
	$line =~ s/\r//g;
	chomp($line);
	my @tmp = split(/\s+/,$line);
	my $uid = shift(@tmp);
	my $rbp = shift(@tmp);
	my $celltype = shift(@tmp);
	
	for my $fi (@tmp) {
	    my @splitfi = split(/\//,$fi);
	    my $shortfi = $splitfi[$#splitfi];
#	    if ($shortfi =~ /^(\d+)\_(\d+)\_(.+)\.merged/) {
	    if ($shortfi =~ /^(\S+?)\.merged/ || $shortfi =~ /^(\S+\.unassigned)\.adapterTrim\.round2\.rmRep\.rmDup\.sorted\.r2\.bam/) {
#		$subset_hash{$1."_".$2."_".$3} = 1;
#		$subset_hash{$1."_INPUT_".$3} = 1;
		$subset_hash{$1} = 1;

#		my $input = $1;
#		$input =~ s/_01_/_INPUT_/;
#		$input =~ s/_02_/_INPUT_/;
#		$subset_hash{$input} = 1;


#	    } else {
#		$short_fi =~ s/adapterTrim\.round2\.rmRep\.rmDup\.sorted\.r2\.bam$//;
#		$subset_hash{$1};
	    }
	}
    }
    close(SUB);
    $subset_flag = 1;
    print STDERR "used subset - should be ".scalar(keys %subset_hash)." files\n";
}



my %hash;
my $merging_fi = $ARGV[1];
# current file - /home/gpratt/projects/encode/scripts/encode_GRCh38_v1.txt
#my $merging_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/encode_GRCh38_v1.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_v13.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_v12.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_master.txt";
#my $merging_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/test_input.txt";
open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    next unless ($line);
	
#    next unless ($line =~ /693_CLIP/);

    my @tmp = split(/\t/,$line);
    my $fastqs = shift(@tmp);
    my $genome = shift(@tmp);
    my $dataset_label = shift(@tmp);
    my ($fastq_r1,$fastq_r2) = split(/\;/,$fastqs);
    my @fastq_split = split(/\//,$fastq_r1);
    my $fastq_short = $fastq_split[$#fastq_split];
    
    if ($subset_flag == 1) {
	if (exists $subset_hash{$dataset_label}) {
	    $subset_hash{$dataset_label} = 2;
	    push @{$hash{$dataset_label}},$fastqs;
	} elsif ($fastq_short =~ /^(\S+\.unassigned)\.fastq\.gz/) {
	    my $input_id = $1;
	    if (exists $subset_hash{$input_id}) {
		$subset_hash{$input_id} = 2;
		push @{$hash{$dataset_label}},$fastqs;
	    }
	}
    } else {
	push @{$hash{$dataset_label}},$fastqs;
    }
#    next unless ($dataset_label =~ /^KB/);

#    print STDERR "dataset_label $dataset_label X $line X\n";
}
close(M);
#exit;

my $subset_count=0;
for my $key (keys %subset_hash) {
    $subset_count++ if ($subset_hash{$key} == 2);
}
print STDERR "actually found $subset_count files\n" if ($subset_flag == 1);
unless ($subset_count == scalar(keys %subset_hash)) {
    for my $key (keys %subset_hash) {
	print STDERR "missing file for $key\n" if ($subset_hash{$key} != 2);
    }
}
	
my $count=0;
my $skipped_count = 0;
my %datasets_listed_twice;
for my $dataset_label (keys %hash) {
    $count++;
#    if ($dataset_label =~ /^452/) {
#	print STDERR "dataset $dataset_label\n";
#    }
#    next unless ($dataset_label eq "676_01_RBFOX2");
    my $rmdup_output_merged_sam = $working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam";

    my $parsed_out = $rmdup_output_merged_sam.".parsed";
    if (exists $datasets_listed_twice{$parsed_out}) {
	print STDERR "dataset listed twice $parsed_out\n";
	$skipped_count++;
	next;
    }
    $datasets_listed_twice{$parsed_out} = 1;
    if (-e $parsed_out) {
	print STDERR "skipping $dataset_label , $parsed_out already exists\n";
	$skipped_count++;
	next;
    }
    
    
    my $sh_out = "RM_".$dataset_label.".sh";
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


    my $merged_rmdup = $working_directory.$dataset_label."_merged.rmDup.sam";

    ## this is gabe's rmDup final bam file



    my @unmerged_sams;
    my @unmerged_parsed;

    for my $fastq (@{$hash{$dataset_label}}) {
	my ($fastq1,$fastq2) = split(/\;/,$fastq);
	for ($fastq1,$fastq2) {
	    my @fastq_split = split(/\//,$_);
	    my $fastq_short = $fastq_split[$#fastq_split];
	    
	    $fastq_short =~ s/\.fastq\.gz$//;
	    my $orig_fi = $current_file_location_prefix.$fastq_short.".adapterTrim.round2.fastq.gz";  
	    $_ = $fastq_short.".adapterTrim.round2.fastq.gz";

	    my $fastq_path = $working_directory.$_;
	 
	    print SH "cp $orig_fi $fastq_path\n";
	    print SH "gunzip $fastq_path\n";
	    $_ = $working_directory.$fastq_short.".adapterTrim.round2.fastq";
	}

	my @fastq_fi1_split = split(/\//,$fastq1);
	my $fastq_fi1_short = $fastq_fi1_split[$#fastq_fi1_split];
	my $r1f1 = $fastq_fi1_short;
	$r1f1 =~ s/\.adapterTrim\.round2\.fastq//;
	$r1f1 .= ".adapterTrim.round2.rmRep.bam";
	my $merged_bam_file_standardanalysis = $current_file_location_prefix.$r1f1;
	print SH "perl split_bam_to_subfiles.pl $merged_bam_file_standardanalysis $working_directory\n";

	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat";
#	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed";
	my @bowtie_db_split = split(/\//,$bowtie_db);
	my $bowtie_db_short = $bowtie_db_split[$#bowtie_db_split];
	

	
	my $sam_output = $working_directory.$fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam";
	my $sam_output_short = $fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam";

	print SH "perl parse_bowtie2_output_realtime_includemultifamily.pl $fastq1 $fastq2 $bowtie_db $working_directory $sam_output\n";
	print SH "rm $fastq1\n" unless ($keep_temp_files_flag==1);
	print SH "rm $fastq2\n" unless ($keep_temp_files_flag==1);

	# now split sam file into files separated by first 2 bases of barcode
	print SH "perl split_bam_to_subfiles.pl $sam_output $working_directory\n";
	
	print SH "rm $sam_output\n" unless ($keep_temp_files_flag==1);

	my $predup_output = $working_directory.$sam_output_short.".combined_w_uniquemap.prermDup.sam";

	my $rmdup_output = $working_directory.$sam_output_short.".combined_w_uniquemap.rmDup.sam";
	print SH "rm $predup_output\n";
	print SH "rm $rmdup_output\n";
#	print SH "rm $rmdup_output\n" if (-e $rmdup_output);

	for my $b1 ("A","C","G","T","N") {
	    for my $b2 ("A","C","G","T","N") {
		my $repmapping_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp";
		my $genomemapping_fi = $working_directory.$r1f1.".".$b1.$b2.".tmp";

#		print SH "perl duplicate_removal_inline_paired.pl $repmapping_fi $genomemapping_fi\n";		
		print SH "perl duplicate_removal_inline_paired.count_region_other_reads.pl $repmapping_fi $genomemapping_fi\n";		
		my $tmp_output_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp".".combined_w_uniquemap.rmDup.sam";
		my $tmp_parsed_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp".".combined_w_uniquemap.rmDup.sam.parsed";
		my $tmp_predup_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp".".combined_w_uniquemap.prermDup.sam";

		print SH "cat $tmp_output_fi >> $rmdup_output\n";
		print SH "cat $tmp_predup_fi >> $predup_output\n";

                print SH "rm $repmapping_fi\n" unless ($keep_temp_files_flag==1);
                print SH "rm $genomemapping_fi\n" unless ($keep_temp_files_flag==1);
                print SH "rm $tmp_output_fi\n" unless ($keep_temp_files_flag==1);
		print SH "rm $tmp_predup_fi\n" unless ($keep_temp_files_flag==1);
		push @unmerged_parsed,$tmp_parsed_fi;
	    }
	}
	
	print SH "gzip $predup_output\n";
	push @unmerged_sams,$rmdup_output;
    }

    print SH "perl merge_multiple_parsed_files.pl $parsed_out ".join(" ",@unmerged_parsed)."\n";
    for my $unmerged_parsed_fi (@unmerged_parsed) {
	print SH "rm $unmerged_parsed_fi\n" unless ($keep_temp_files_flag==1);
    }


    print SH "rm $rmdup_output_merged_sam\n";
    my $merge_command = "cp $unmerged_sams[0] $rmdup_output_merged_sam\n";
    for (my $i=1;$i<scalar(@unmerged_sams);$i++) {
	$merge_command .= "cat $unmerged_sams[$i] >> $rmdup_output_merged_sam\n";
    }
    print SH "$merge_command";
    
    for my $rmdup_output_sam (@unmerged_sams) {
	print SH "rm $rmdup_output_sam\n" unless ($keep_temp_files_flag==1);
    }
    print SH "gzip $rmdup_output_merged_sam\n";
## count repbase mapping    

    close(SH);

#    exit;

    system("qsub $sh_out");

    if ($count % 20 == 0) {
#	exit;
	sleep 5;
    }

    if ($count > 10) {
#	exit;
    }



    my $flag = 1;
    while ($flag == 1) {
	my $current = 0;
	my $queue = `qstat -u elvannostrand`;
	my @queuee = split(/\n/,$queue);
	for my $line (@queuee) {
	    my @tmp = split(/\s+/,$line);
	    next unless ($tmp[9]);
	    my $queue_flag = $tmp[9];
#               next if ($tmp[2] eq "home-yeo");
	    next unless ($tmp[2] eq "condo");
	    $current++ if ($queue_flag eq "Q" || $queue_flag eq "R");
	}
	
	print STDERR "current_queued $current\n";
	
	if ($current >= 80) {
#               $count = 0;
	    sleep 100;
	}
	if ($current < 80) {
	    $current++;
	    $flag = 0;
	}
    }
    
}

print STDERR "skipped $skipped_count datasets that were already done\n";
