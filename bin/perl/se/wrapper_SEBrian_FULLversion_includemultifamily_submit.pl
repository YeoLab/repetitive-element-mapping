use warnings;
use strict;

unless ($ARGV[0] && $ARGV[1]) {
    print STDERR "usage: perl wrapper_FULLversion_includemultifamily_submit.pl working_directory manifest_fi fastq_folder_prefix(optional-defaults to encode)\n";
    print STDERR "note - manifest_fi is gabe's pipeline file for merging fastq files\n";
    print STDERR "\n";
    print STDERR "current files used are:\n";
    print STDERR "0 working directory = /home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/\n";
    print STDERR "1 current manifest = /home/elvannostrand/data/clip/CLIPseq_analysis/Collaborations/VSV/infected_manifest.txt";

    exit;
}

my $keep_temp_files_flag = 1;



my $working_directory = $ARGV[0];
# working directory is currently /home/elvannostrand/scratch/ALLCLIP_v12_repeatmapping/

$working_directory =~ s/\/$//;
$working_directory .= "/";

system("mkdir $working_directory") unless (-d $working_directory);
#my $working_directory = "/home/elvannostrand/scratch/ENCODE_v9_20151209/PolIII_mapping/";
#my $working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_v9_20151209/PolIII_mapping/";

my %subset_hash;
my $subset_flag = 0;

my %hash;
my $merging_fi = $ARGV[1];
open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    $line =~ s/\r//g;
    next unless ($line);
    my @tmp = split(/\t/,$line);
    my $uid = shift(@tmp);
    my $rbp = shift(@tmp);
    my $cellline = shift(@tmp);

    for my $bamfi (@tmp) {
	my @splitfi = split(/\//,$bamfi);
	my $shortfi = $splitfi[$#splitfi];

	my $dataset_label = "";
	if ($shortfi =~ /^(.+)\.umi\.r1TrTr\.sorted\.STARUnmapped\.out\.sorted\.STARAligned\.outSo\.rmDupSo\.bam$/) {
	    $dataset_label = $1;
	    $hash{$uid."|".$dataset_label} = $bamfi;
	}
    }
    
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
#dataset_label = "encode4_batch3.SB06_CLIP";
for my $uid_dataset_label (keys %hash) {
    $count++;

    my ($uid,$dataset_label) = split(/\|/,$uid_dataset_label);
    my $bamfi = $hash{$uid_dataset_label};
    my @splitfi = split(/\//,$bamfi);
    pop(@splitfi);
    my $current_file_location_prefix = join("/",@splitfi)."/";


#    next unless ($dataset_label eq "676_01_RBFOX2");
    my $merged_rmdup = $working_directory.$dataset_label."_merged.rmDup.sam";
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
    print SH "cd /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/SEpipeline_Brianfilestructure/\n";

    my @unmerged_sams;
    my @unmerged_parsed;

    my $fastq_short = $dataset_label.".umi.r1TrTr.fq.gz";
    my $fastq_unzipped_short = $dataset_label.".umi.r1TrTr.fq";
    my $fastq_original = $current_file_location_prefix.$fastq_short;
    my $fastq_unzipped_orig = $current_file_location_prefix.$fastq_unzipped_short;
    my $fastq_path = $working_directory.$fastq_short;
    my $fastq_unzipped = $working_directory.$fastq_unzipped_short;

    if (-e $fastq_original) {
	print SH "cp $fastq_original $fastq_path\n";
	print SH "gunzip $fastq_path\n";
    } elsif (-e $fastq_unzipped_orig) {
	print SH "cp $fastq_unzipped_orig $fastq_unzipped\n";
    } else {
	print STDERR "ERROR ERROR - missing both original fastq and unzipped original fastq\n$fastq_original\n$fastq_unzipped_orig\n";
    }



    
    my $merged_bam_file_standardanalysis = $current_file_location_prefix.$dataset_label.".umi.r1TrTr.sorted.STARUnmapped.out.sorted.STARAligned.out.bam";
    print SH "perl split_bam_to_subfiles_SE.pl $merged_bam_file_standardanalysis $working_directory\n";
    
    my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat";
#	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed";
    my @bowtie_db_split = split(/\//,$bowtie_db);
    my $bowtie_db_short = $bowtie_db_split[$#bowtie_db_split];
    
    my $sam_output = $working_directory.$fastq_short.".mapped_vs_".$bowtie_db_short.".sam";
    my $sam_output_short = $fastq_short.".mapped_vs_".$bowtie_db_short.".sam";
    
    print SH "perl parse_bowtie2_output_realtime_includemultifamily_SE.pl $fastq_unzipped $bowtie_db $working_directory $sam_output\n";
#    print SH "perl parse_bowtie2_output_realtime_includemultifamily.pl $fastq1 $fastq2 $bowtie_db $working_directory $sam_output\n";
    print SH "rm $fastq_unzipped\n" unless ($keep_temp_files_flag==1);
    
    # now split sam file into files separated by first 2 bases of barcode
    print SH "perl split_bam_to_subfiles_SE.pl $sam_output $working_directory\n";
    print SH "rm $sam_output\n" unless ($keep_temp_files_flag==1);

    my $predup_output = $working_directory.$sam_output_short.".combined_w_uniquemap.prermDup.sam";
    my $rmdup_output = $working_directory.$sam_output_short.".combined_w_uniquemap.rmDup.sam";
    print SH "rm $predup_output\n";
    print SH "rm $rmdup_output\n";
#	print SH "rm $rmdup_output\n" if (-e $rmdup_output);
    
    for my $b1 ("A","C","G","T","N") {
	for my $b2 ("A","C","G","T","N") {
	    my $repmapping_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp";
	    my $genomemapping_fi = $working_directory.$dataset_label.".umi.r1TrTr.sorted.STARUnmapped.out.sorted.STARAligned.out.bam".".".$b1.$b2.".tmp";
	    
#		print SH "perl duplicate_removal_inline_paired.pl $repmapping_fi $genomemapping_fi\n";		
	    print SH "perl duplicate_removal_inline_paired.count_region_other_reads_SE.pl $repmapping_fi $genomemapping_fi\n";		
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
    my $rmdup_output_merged_sam_gzip = $rmdup_output_merged_sam.".gz";
    print SH "perl reparse_samfile_updatedchrM_fixmultenstsort_SE.pl $rmdup_output_merged_sam_gzip\n";

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
