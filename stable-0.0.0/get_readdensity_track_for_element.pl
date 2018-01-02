use warnings;
use strict;

#my $desired_element = "L1||L1HS";
my $desired_element = $ARGV[1];
#my $desired_element = "RNA18S||ENST00000606783.1";
my ($desfamily,$desenst) = split(/\|\|/,$desired_element);
#print STDERR "desenst $desenst\n";
my $samfi = $ARGV[0];
my $outsam = $samfi.".tmp.".$desfamily.".sam";
my $outbam = $samfi.".tmp.".$desfamily.".bam";
my $sorted = $outbam.".sorted.bam";
my $mpileup_out = $sorted.".mpileup";
if (-e $mpileup_out) {
    print STDERR "skipping $samfi - already done\n";
    exit;
}

unless (-e $outsam) {
    open(O,">$outsam");
    open(S,$samfi);
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


#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/L1HS_only.fa";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/18S_only.fa";
my $desired_fasta = $ARGV[2];
my $desired_fai = $desired_fasta.".fai";



unless (-e $desired_fai) {
    system("samtools faidx $desired_fasta");
}



unless (-e $outbam) {
    system("samtools view -bt $desired_fai $outsam > $outbam");
}
my $tmp_sorting = $outbam.".tmp";
unless (-e $sorted) {
    system("samtools sort $outbam -T $tmp_sorting -o $sorted");
}

unless (-e $mpileup_out) {
    system("samtools mpileup -d 100000000 -f $desired_fasta $sorted > $mpileup_out");
}

system("rm $outsam");
system("rm $outbam");
