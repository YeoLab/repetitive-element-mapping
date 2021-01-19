use warnings;
use strict;

my $species = "hg38";
## note - this pipeline right now is specific for hg38

my $keep_temp_files_flag = 0;

##
unless ($ARGV[4]) {
    print STDERR "usage: perl RepElement_pipeline_1dataset.pl working_directory fastq1_R1:fastq1_R2::fastq2_R1:fastq2_R2::(etc) genomic_mapping_bam output_prefix se_or_pe_flag\n";
    print STDERR "e.g.\n";
    print STDERR "PE example: perl RepElement_pipeline_1dataset.pl /home/elvannostrand/scratch/20201211_PEtest_204/ /projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/204_RBFOX2_GRCh38/results/204_RBFOX2.rep1_clip.A01.r1.fqTrTr.sorted.fq.gz:/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/204_RBFOX2_GRCh38/results/204_RBFOX2.rep1_clip.A01.r2.fqTrTr.sorted.fq.gz::/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/204_RBFOX2_GRCh38/results/204_RBFOX2.rep1_clip.B06.r1.fqTrTr.sorted.fq.gz:/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/204_RBFOX2_GRCh38/results/204_RBFOX2.rep1_clip.B06.r2.fqTrTr.sorted.fq.gz /projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/204_RBFOX2_GRCh38/results/204_RBFOX2.rep1_clip.A01.r1.fq.genome-mappedSo.bam::/projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/204_RBFOX2_GRCh38/results/204_RBFOX2.rep1_clip.B06.r1.fq.genome-mappedSo.bam 204_01 PE\n";
    print STDERR "SE example: perl RepElement_pipeline_1dataset.pl /home/elvannostrand/scratch/20201211_SEtest/ /projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/INV_B.IP.umi.r1.fqTrTr.sorted.fq.gz /projects/ps-yeolab4/software/eclip/0.7.0/tests/other_examples/INV_B_singleNode/results/INV_B.IP.umi.r1.fq.genome-mappedSoSo.bam RBFOX_seCLIP SE\n";
    exit;
}

my $working_directory = $ARGV[0];
my $fastqs_input = $ARGV[1];
my $standard_analysis_bams_input = $ARGV[2];
my $dataset_label = $ARGV[3];
my $se_or_pe_flag = $ARGV[4];
my @fastq_demux_pairs = split(/\:\:/,$fastqs_input);
my @standard_analysis_bams_split = split(/\:\:/,$standard_analysis_bams_input);

my $rmdup_output_merged_sam = $working_directory.$dataset_label.".combined_w_uniquemap.rmDup.sam";
my $merged_rmdup = $working_directory.$dataset_label."_merged.rmDup.sam";
my $parsed_out = $rmdup_output_merged_sam.".parsed";

if (-e $parsed_out) {
    print STDERR "skipping $dataset_label , $parsed_out already exists\n";
    next;
}

if (-d $working_directory) {
} else {
    system("mkdir $working_directory");
}
$working_directory =~ s/\/$//;
$working_directory .= "/";
    
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
print SH "cd /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/hg38/\n";


my @unmerged_sams;
my @unmerged_parsed;

for my $fastq (@fastq_demux_pairs) {
    my @fastqpairs = split(/\:/,$fastq);
    if (scalar(@fastqpairs) > 2) {
	print STDERR "error - this should be either 1 fastq (single-end eCLIP) or 2 (R1 and R2 for original paired-end eCLIP), but got ".scalar(@fastqpairs)." instead\n$fastq\n";
	exit;
    }

    foreach (@fastqpairs) {
	my @fastq_split = split(/\//,$_);
	my $fastq_short = $fastq_split[$#fastq_split];
	my $newpath = $working_directory.$fastq_short;
	
	print SH "cp $_ $newpath\n";
	print SH "gunzip $newpath\n";

	$fastq_short =~ s/\.gz$//;
	$_ = $working_directory.$fastq_short;
    }
    my $fastq1 = $fastqpairs[0];
    my $fastq2;
    $fastq2 = $fastqpairs[1] if ($se_or_pe_flag eq "PE");

    my @fastq_fi1_split = split(/\//,$fastq1);
    my $fastq_fi1_short = $fastq_fi1_split[$#fastq_fi1_split];
    my $merged_bam_file_standardanalysis = shift(@standard_analysis_bams_split);
    my @merged_bam_file_split = split(/\//,$merged_bam_file_standardanalysis);
    my $merged_bam_file_standardanalysis_short = $merged_bam_file_split[$#merged_bam_file_split];

    print SH "perl split_bam_to_subfiles_SEorPE.pl $merged_bam_file_standardanalysis $working_directory $se_or_pe_flag\n";
    
#	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat";
    #hg19
    my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.20190424";
    if ($species eq "hg38") {
	$bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/hg38/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat";
    }
#	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed";
    my @bowtie_db_split = split(/\//,$bowtie_db);
    my $bowtie_db_short = $bowtie_db_split[$#bowtie_db_split];
    
    my $sam_output = $working_directory.$fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam";
    my $sam_output_short = $fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam";
    
    if ($se_or_pe_flag eq "PE") {
	print SH "perl parse_bowtie2_output_realtime_includemultifamily_PE.pl $fastq1 $fastq2 $bowtie_db $working_directory $sam_output $species\n";
    } else {
	print SH "perl parse_bowtie2_output_realtime_includemultifamily_SE.pl $fastq1 $bowtie_db $working_directory $sam_output $species\n";
    }
    print SH "rm $fastq1\n" unless ($keep_temp_files_flag==1);    
    if ($se_or_pe_flag eq "PE") {
	print SH "rm $fastq2\n" unless ($keep_temp_files_flag==1);
    }

    # now split sam file into files separated by first 2 bases of barcode
    print SH "perl split_bam_to_subfiles_SEorPE.pl $sam_output $working_directory $se_or_pe_flag\n";
    
    print SH "rm $sam_output\n" unless ($keep_temp_files_flag==1);
    
    my $predup_output = $working_directory.$sam_output_short.".combined_w_uniquemap.prermDup.sam";
    
    my $rmdup_output = $working_directory.$sam_output_short.".combined_w_uniquemap.rmDup.sam";
    print SH "rm $predup_output\n";
    print SH "rm $rmdup_output\n";
#	print SH "rm $rmdup_output\n" if (-e $rmdup_output);
    
    for my $b1 ("A","C","G","T","N") {
	for my $b2 ("A","C","G","T","N") {
	    my $repmapping_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp";
	    my $genomemapping_fi = $working_directory.$merged_bam_file_standardanalysis_short.".".$b1.$b2.".tmp";
	    
	    print SH "perl duplicate_removal_inline_paired.count_region_other_reads_masksnRNAs_andreparse_SEandPE_20201210_simple.pl $repmapping_fi $genomemapping_fi $se_or_pe_flag\n";   
	    my $tmp_output_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp".".combined_w_uniquemap.rmDup.sam";
	    my $tmp_parsed_fi = $working_directory.$sam_output_short.".".$b1.$b2.".tmp".".combined_w_uniquemap.rmDup.sam.parsed_v2.20201210.txt";
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

print SH "perl merge_multiple_parsed_files.simplified_20191022.pl $parsed_out ".join(" ",@unmerged_parsed)."\n";
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

system("qsub $sh_out");
