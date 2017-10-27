use warnings;
use strict;

my $keep_temp_files_flag = 1;

my $current_file_location_prefix = "/projects/ps-yeolab3/encode/analysis/encode_master/";
#my $current_file_location_prefix = "/projects/ps-yeolab2/encode/analysis/encode_master/";

unless ($ARGV[0]) {
    print STDERR "usage perl wrapper_FULLversion_submit.pl working_directory\n";
    exit;
}

my $working_directory = $ARGV[0];
# working directory is currently /home/elvannostrand/scratch/ALLCLIP_v12_repeatmapping/

$working_directory =~ s/\/$//;
$working_directory .= "/";

system("mkdir $working_directory") unless (-d $working_directory);
#my $working_directory = "/home/elvannostrand/scratch/ENCODE_v9_20151209/PolIII_mapping/";
#my $working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_v9_20151209/PolIII_mapping/";

my %hash;
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_v13.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_v12.txt";
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_master.txt";
my $merging_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/test_input.txt";
open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $fastqs = shift(@tmp);
    my $genome = shift(@tmp);
    my $dataset_label = shift(@tmp);
#    next unless ($dataset_label =~ /^KB/);
    push @{$hash{$dataset_label}},$fastqs;
    print STDERR "dataset_label $dataset_label X $line X\n";
}
close(M);
#exit;
	
my $count=0;

for my $dataset_label (keys %hash) {
    $count++;
#    if ($dataset_label =~ /^452/) {
#	print STDERR "dataset $dataset_label\n";
#    }
#    next unless ($dataset_label =~ /^452_01/);
    my $rmdup_output_merged_bam = $working_directory.$dataset_label."_RepMappermerged.rmDup.bam";
    my $rmdup_output_merged_namesort_bam = $rmdup_output_merged_bam.".namesort.bam";
    my $parsed_out = $rmdup_output_merged_namesort_bam.".parsed";
    if (-e $parsed_out) {
#	print STDERR "skipping $dataset_label , $parsed_out already exists\n";
	next;
    }
    
    
    my $sh_out = "RepeatMapper_".$dataset_label.".sh";
    open(SH,">$sh_out");
    print SH "\#\!\/bin\/sh\n";
    print SH "#PBS -N ".$sh_out."\n";
    print SH "#PBS -o ".$sh_out.".out\n";
    print SH "#PBS -e ".$sh_out.".err\n";
    print SH "#PBS -V\n";
    print SH "#PBS -l walltime=24:00:00\n";
    print SH "#PBS -l nodes=1:ppn=16\n";
    print SH "#PBS -A yeo-group\n";
    print SH "#PBS -q home-yeo\n";
    print SH "#PBS -M elvannostrand\@ucsd.edu\n";
    print SH "#PBS -m a\n";
    print SH "cd /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/\n";


    my $merged_rmdup = $working_directory.$dataset_label."_merged.rmDup.sam";

    my @unmerged_bams;
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

	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed";
	my @bowtie_db_split = split(/\//,$bowtie_db);
	my $bowtie_db_short = $bowtie_db_split[$#bowtie_db_split];
	
	my @fastq_fi1_split = split(/\//,$fastq1);
	my $fastq_fi1_short = $fastq_fi1_split[$#fastq_fi1_split];
	
	my $sam_output = $working_directory.$fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam";

	print SH "perl parse_bowtie2_output_realtime.pl $fastq1 $fastq2 $bowtie_db $working_directory $sam_output\n";
	print SH "rm $fastq1\n" unless ($keep_temp_files_flag==1);
	print SH "rm $fastq2\n" unless ($keep_temp_files_flag==1);
	print SH "perl duplicate_remove_inline.pl $sam_output\n";
	my $rmdup_output = $working_directory.$fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam.rmDup.sam";
	my $rmdup_output_bam = $working_directory.$fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam.rmDup.bam";
	print SH "samtools view -b $rmdup_output > $rmdup_output_bam\n";
	push @unmerged_bams,$rmdup_output_bam;
	print SH "rm $sam_output\n" unless ($keep_temp_files_flag==1);
	print SH "rm $rmdup_output\n" unless ($keep_temp_files_flag==1);

    }


    print SH "rm $rmdup_output_merged_bam\n" if (-e $rmdup_output_merged_bam);
    my $merge_command = "samtools merge $rmdup_output_merged_bam ".join(" ",@unmerged_bams);
    print SH "$merge_command\n";
    
    for my $rmdup_output_bam (@unmerged_bams) {
	print SH "rm $rmdup_output_bam\n" unless ($keep_temp_files_flag==1);
    }

    my $temp_resortmerge_prefix = $rmdup_output_merged_bam.".tmp";

    print SH "rm $rmdup_output_merged_namesort_bam\n" if (-e $rmdup_output_merged_namesort_bam);
    print SH "samtools sort -n -O bam -T $temp_resortmerge_prefix -o $rmdup_output_merged_namesort_bam $rmdup_output_merged_bam\n";
    print SH "rm $rmdup_output_merged_bam\n" unless ($keep_temp_files_flag==1);

    my $merged_bam_file_standardanalysis = $current_file_location_prefix.$dataset_label.".merged.bam";
    ## this is gabe's rmDup final bam file
    if (scalar(@{$hash{$dataset_label}}) == 1) { 
# this is for inputs - note that this only happens to work because inputs have 1 file but clips have 2, this is a temporary hack
	my $fastq = $hash{$dataset_label}[0];
	my ($fastq1,$fastq2) = split(/\;/,$fastq);
	my @fastq1_split = split(/\//,$fastq1);
	my $fastq1_short = $fastq1_split[$#fastq1_split];
	my $r1f1 = $fastq1_short;
	$r1f1 =~ s/\.fastq.gz//;
	$r1f1 .= ".adapterTrim.round2.rmRep.rmDup.bam";
	$merged_bam_file_standardanalysis = $current_file_location_prefix.$r1f1;
    }
    
    my $temp_resort_prefix = $working_directory.$dataset_label.".merged.bam.nameresort_tmp";
    my $resorted_fi = $working_directory.$dataset_label.".merged.bam.nameresort.bam";
    print SH "rm $resorted_fi\n" if (-e $resorted_fi);
    print SH "samtools sort -n -O bam -T $temp_resort_prefix -o $resorted_fi $merged_bam_file_standardanalysis\n";

## count repbase mapping    

    print SH "perl count_from_inline_rmDup_fi_andcountbothrmDupfiles.pl $rmdup_output_merged_namesort_bam $resorted_fi > $parsed_out\n";

    print SH "rm $resorted_fi\n" unless ($keep_temp_files_flag==1);
    print SH "rm $rmdup_output_merged_namesort_bam\n" unless ($keep_temp_files_flag==1);

    close(SH);

    system("qsub $sh_out");
    exit;

    if ($count % 20 == 0) {
#	exit;
	sleep 5;
    }
}


