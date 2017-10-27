use warnings;
use strict;

my $removed_fi=0;
my $done_fi_list = `ls /home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/*.done`;
my @done_fi_listt = split(/\s+/,$done_fi_list);
for my $done_fi (@done_fi_listt) {
    my $full_text = `cat $done_fi`;
    chomp($full_text);
    if ($full_text eq "jobs done") {
	system("rm $done_fi");
	$removed_fi++;
	print STDERR "removed $removed_fi\r";
    }
}
