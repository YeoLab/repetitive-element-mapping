use warnings;
use strict;

my $desired_element = "RNA18S||ENST00000606783.1";
my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/18S_only.fa";  
#my $desired_element = "RNA28S||ENST00000607521.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/28S_only.fa";    
#my $desired_element = "rRNA_extra||NR_046235.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/rRNA_extra.fa";
#my $desired_element = "RNA5-8S||ENST00000474885.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/5.8S_only.fa";

my ($desfamily,$desenst) = split(/\|\|/,$desired_element);
my %seq;
&read_seq($desired_fasta,$desenst);

my $filist_fi = "rRNA_filelist_20170302_inclnonsubmitted.csv";

my $output_file = $filist_fi.".".$desfamily.".csv";

unless (exists $seq{$desenst}) {
    print STDERR "error $desenst not found in $desired_fasta\n";
    exit;
}
my @seqq = split(//,$seq{$desenst});

my $output2 = $output_file.".wseq.csv";
open(OUT,">$output2");
print OUT "ID\tUID";
for (my $i=0;$i<@seqq;$i++) {
    print OUT "\t".($i+1);
}
print OUT "\n";
print OUT "ID\tUID\t".join("\t",@seqq)."\n";
close(OUT);

system("cat $output_file >> $output2");

sub read_seq {
    my $file = shift;
    my $id;
    open(F,$file);
    for my $line (<F>) {
	chomp($line);
	if ($line =~ /^\>(.+)$/) {
	    $id = $1;
	    
	} else {
	    $seq{$id} .= $line; 
	}
    }
    close(F);

}
