use warnings;
use strict;

my $current_file_location_prefix = "/projects/ps-yeolab2/encode/analysis/encode_master/";

unless ($ARGV[0]) {
    print STDERR "usage perl wrapper_FULLversion_submit.pl working_directory\n";
    exit;
}

my $working_directory = $ARGV[0];
# working directory now is /home/elvannostrand/scratch/ALLCLIP_v12_repeatmapping/

$working_directory =~ s/\/$//;
$working_directory .= "/";

system("mkdir $working_directory") unless (-d $working_directory);
#my $working_directory = "/home/elvannostrand/scratch/ENCODE_v9_20151209/PolIII_mapping/";
#my $working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_v9_20151209/PolIII_mapping/";

my %hash;
#my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_v12.txt";
my $merging_fi = "/home/gpratt/projects/encode/scripts/encode_master.txt";
open(M,$merging_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    my $fastqs = shift(@tmp);
    my $genome = shift(@tmp);
    my $dataset_label = shift(@tmp);
#    next unless ($dataset_label =~ /^KB/);
    push @{$hash{$dataset_label}},$fastqs;
}
close(M);

my $removed_file_count=0;

for my $dataset_label (keys %hash) {

    my $rmdup_output_merged_bam = $working_directory.$dataset_label."_RepMappermerged.rmDup.bam";
    my $rmdup_output_merged_namesort_bam = $rmdup_output_merged_bam.".namesort.bam";
    my $parsed_out = $rmdup_output_merged_namesort_bam.".parsed";
    
    
    my $sh_out = "RepeatMapper_".$dataset_label.".sh";
    my $sh_err = $sh_out.".err";
    my $merged_rmdup = $working_directory.$dataset_label."_merged.rmDup.sam";
    
    my $remove_flag = 0;
    unless (-e $sh_err) {
#	print STDERR "no $sh_err file - skipping\n";
	next;
    }
    open(F,$sh_err);
    my $firstline = (<F>);
    unless ($firstline) {
#        print STDERR "no .err file output $sh_err\n";
        next;
    }
    my @linesplit = split(/\s+/,$firstline);
    
    if ($linesplit[4] && ($linesplit[4] eq "Killed")) {
	print "sh_err\n$firstline\n";
	
	my $rmcommand2 = "rm $sh_err";
	print STDERR "$rmcommand2\n";
	my $rmcommand = "rm $parsed_out";
	print STDERR "$rmcommand\n";
	system($rmcommand);
	system($rmcommand2);
	$removed_file_count++;
	$remove_flag = 1;
    } else {

        if ($firstline =~ /Killed/) {
            print STDERR "weird firstline $firstline $sh_err\n";
        }
        for my $line (<F>) {
            chomp($line);

            if ($line =~ /Killed/) {
                print STDERR "weird $line $sh_err\n";
            }
        }
    }


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
	    $_ = $working_directory.$fastq_short.".adapterTrim.round2.fastq";
	    
	}

	my $bowtie_db = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.fa.fixed";
	my @bowtie_db_split = split(/\//,$bowtie_db);
	my $bowtie_db_short = $bowtie_db_split[$#bowtie_db_split];
	
	my @fastq_fi1_split = split(/\//,$fastq1);
	my $fastq_fi1_short = $fastq_fi1_split[$#fastq_fi1_split];
	
	my $sam_output = $working_directory.$fastq_fi1_short.".mapped_vs_".$bowtie_db_short.".sam";

	my $bowtie_out = $sam_output.".bowtieout";
	
	if ($remove_flag == 1) {
	    my $rmcommand = "rm $bowtie_out";
	    system($rmcommand);
	    print STDERR "$rmcommand\n";
	}
    }

}


print STDERR "removed $removed_file_count files\n";
