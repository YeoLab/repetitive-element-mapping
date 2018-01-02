use warnings;
use strict;

my $filist_fi = $ARGV[0];
my $desired_element = "RNA18S||ENST00000606783.1";
my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/18S_only.fa";  
#my $desired_element = "RNA28S||ENST00000607521.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/28S_only.fa";    
#my $desired_element = "rRNA_extra||NR_046235.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/rRNA_extra.fa";

#my $desired_element = "RNA5-8S||ENST00000474885.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/5.8S_only.fa";

#my $working_dir = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170114/";
#my $working_dir = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
#my $working_dir = "/home/elvannostrand/scratch/Satellite/";
my $working_dir = "/home/elvannostrand/scratch/Enching_RepElement/";
print STDERR "doing $desired_element\n";

my $inputpairingfi = "Enching_list.txt";
#my $inputpairingfi = "Satellite_inputpairing_20170419.txt";
#my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20170227.txt";
my %clip2input;

open(IN,$inputpairingfi);
for my $line (<IN>) {
    chomp($line);
    $line =~ s/\r//g;

    my @tmp = split(/\s+/,$line);
    
    if ($tmp[3] eq "LNG" || $tmp[3] eq "TAG" || $tmp[3] eq "SingleRep" || $tmp[3] eq "OneRep") {
        $clip2input{$tmp[1]} = $tmp[2];
    } else {
        $clip2input{$tmp[1]} = $tmp[3];
        $clip2input{$tmp[2]} = $tmp[3] unless ($tmp[2] eq "NoInput");
    }

}
close(IN);

my %submitted_hash;



open(FL,$filist_fi);
for my $line (<FL>) {
    chomp($line);
    next unless ($line);
    my ($uid,$rep,$rbp) = split(/\_/,$line);

#    if ($rbp eq "NOLC1") {
#	$rbp = "BOP1";
#    }
#    $line = $uid."_".$rep."_".$rbp;
#    my $input = $uid."_INPUT_".$rbp;

    my $input = $clip2input{$line};

    my $rbp_fi = $working_dir.$line.".combined_w_uniquemap.rmDup.sam";
    my $input_fi = $working_dir.$input.".combined_w_uniquemap.rmDup.sam";

    print STDERR "now doign $rbp_fi\n";
    system("perl","get_readdensity_track_for_element.pl","$rbp_fi","$desired_element","$desired_fasta");
    print STDERR "now doign $input_fi\n";
    system("perl","get_readdensity_track_for_element.pl","$input_fi","$desired_element","$desired_fasta");



}
close(FL);
