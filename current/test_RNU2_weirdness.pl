use warnings;
use strict;

my $fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat";
open(F,$fi);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    
    if ($tmp[0] eq "ENST00000618664.1") {
	for (my $i=0;$i<@tmp;$i++) {
	    print "X$i\t".$tmp[$i]."X\n";
	}
    }
}
close(F);
