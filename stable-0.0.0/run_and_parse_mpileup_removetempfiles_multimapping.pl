use warnings;
use strict;

#my $desired_element = "rRNA_extra||NR_046235.1";                                                                                                    
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/rRNA_extra.fa";  

#

my $samfi = $ARGV[0];
my $desired_element = $ARGV[1];
my $desired_fasta = $ARGV[2];
my ($desfamily,$desenst) = split(/\|\|/,$desired_element);

my $keep_files_flag = 0;
$keep_files_flag = $ARGV[3] if ($ARGV[3]);

my $outsam = $samfi.".tmp.".$desfamily.".sam";
my $outbam = $samfi.".tmp.".$desfamily.".bam";
my $sorted = $outbam.".sorted.bam";
my $mpileup_out = $sorted.".mpileup";
my $mpileup_shortout = $sorted.".mpileup.short";

if (-e $mpileup_shortout) {
    print STDERR "skipping $samfi - already done\n";
    exit;
}

unless (-e $outsam) {
    open(O,">$outsam");
    if ($samfi =~ /\.sam$/) {
	open(S,$samfi);
    } elsif ($samfi =~ /\.sam\.gz$/) {
	open(S,"gunzip -c $samfi |");
    } else {
	print STDERR "error - $samfi not either .sam or .sam.gz\n";
	exit;
    }
    while (<S>) {
        my $line = $_;
#    for my $line (<S>) {                                                                                                                                                                                            
        chomp($line);
        my @tmp = split(/\t/,$line);
        next unless ($tmp[2] eq $desired_element);

        $tmp[2] = $desenst;

        my $del = pop(@tmp);
        my $del2 = pop(@tmp);

        print O "".join("\t",@tmp)."\n";
    }
    close(S);
    close(O);
}


my $desired_fai = $desired_fasta.".fai";

unless (-e $desired_fai) {
    system("samtools faidx $desired_fasta");
}
unless (-e $outbam) {
    system("samtools view -bt $desired_fai $outsam > $outbam");
}
system("rm $outsam") unless ($keep_files_flag == 1);
my $tmp_sorting = $outbam.".tmp";
unless (-e $sorted) {
    system("samtools sort $outbam -T $tmp_sorting -o $sorted");
}
system("rm $outbam") unless ($keep_files_flag == 1);


unless (-e $mpileup_out) {
    system("samtools mpileup -d 100000000 -f $desired_fasta $sorted > $mpileup_out");
    system("cut -f1-4 $mpileup_out > $mpileup_shortout");
}
system("rm $sorted") unless ($keep_files_flag == 1);
system("rm $mpileup_out") unless ($keep_files_flag == 1);

