use warnings;
use strict;
use Statistics::Basic qw(:all);
use Statistics::Distributions;
use Statistics::R;
my $R = Statistics::R->new() ;

#my $fi_ip = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/539_01_NOL12.combined_w_uniquemap.rmDup.sam.parsed.reparse_primary.txt";
my $fi_ip = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/539_02_NOL12.combined_w_uniquemap.rmDup.sam.parsed.reparse_primary.txt";
my $fi_input = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/539_INPUT_NOL12.combined_w_uniquemap.rmDup.sam.parsed.reparse_primary.txt";

my %mapped_read_count;
my %hash;
&readfi($fi_ip,"IP");
&readfi($fi_input,"INPUT");



my %precalculated_fisher;

for my $element (keys %hash) {
    unless (exists $hash{$element}{"INPUT"}{read_n}) {
	$hash{$element}{"INPUT"}{read_n} = 1;
	$hash{$element}{"INPUT"}{read_rpr}  = 0;
    }
    next unless (exists $hash{$element}{"IP"}{read_n} && $hash{$element}{"IP"}{read_n} > 0);
    my $l2fc = log(($hash{$element}{"IP"}{read_n}/$mapped_read_count{"IP"}{"UsableReads"}) / ($hash{$element}{"INPUT"}{read_n}/$mapped_read_count{"INPUT"}{"UsableReads"})) / log(2);

    my ($chipval,$chival,$chitype,$chienrdepl) = &fisher_or_chisq($hash{$element}{"IP"}{read_n},$mapped_read_count{"IP"}{"UsableReads"}-$hash{$element}{"IP"}{read_n},$hash{$element}{"INPUT"}{read_n},$mapped_read_count{"INPUT"}{"UsableReads"}-$hash{$element}{"INPUT"}{read_n});
    my $log10pval = $chipval > 0 ? -1 * log($chipval)/log(10) : 400 ;

    print "$element\t".$hash{$element}{"IP"}{read_rpr}."\t".$hash{$element}{"INPUT"}{read_rpr}."\t$l2fc\t$chipval\t$chival\t$chitype\t$chienrdepl\t$log10pval\n";
}




sub readfi {
    my $fi = shift;
    my $label = shift;

    open(F,$fi);
    for my $line (<F>) {
	chomp($line);
	my @tmp = split(/\t/,$line);

	if ($tmp[0] eq "#READINFO") {
	    $mapped_read_count{$label}{$tmp[1]} = $tmp[2];
	} elsif ($tmp[0] eq "TOTAL") {
	} elsif ($tmp[0] eq "PRIMARYELEMENT") {
	    $hash{$tmp[4]}{$label}{read_n} = $tmp[2];
	    $hash{$tmp[4]}{$label}{read_rpr} = $tmp[3];
	    $hash{$tmp[4]}{family} = $tmp[1];
	}
    }
    close(F);


}




sub abs {
    my $x = shift;
    if ($x > 0) {
        return($x);
    } else {
        return(-1*$x);
    }
}

sub square {
    my $x = shift;
    return($x * $x);
}



sub fisher_or_chisq {
    my ($a,$b,$c,$d) = @_;
    unless ($a && $b && $c && $d) {
#        return("1","NA","NA");
    }

    my $tot = $a + $b + $c + $d;
    my $expa = ($a+$c)*($a+$b)/$tot;
    my $expb = ($b+$d)*($a+$b)/$tot;
    my $expc = ($a+$c)*($c+$d)/$tot;
    my $expd = ($b+$d)*($c+$d)/$tot;

    my $direction = "enriched";
    if ($a<$expa) {
        $direction = "depleted";
        return(1,"DEPL","N",$direction);
    }



    if ($expa < 5 || $expb < 5 || $expc < 5 || $expd < 5 || $a < 5 || $b < 5 || $c < 5 || $d < 5) {
        if (exists $precalculated_fisher{$a."|".$c}) {
            return($precalculated_fisher{$a."|".$c}{p},$precalculated_fisher{$a."|".$c}{v},"F",$direction);
        } else {
            my ($pval,$val) = &fisher_exact($a,$b,$c,$d);
            $precalculated_fisher{$a."|".$c}{p} = $pval;
            $precalculated_fisher{$a."|".$c}{v} = $val;
            return($pval,$val,"F",$direction);
        }
    } else {
        my ($pval,$val) = &chi_square($a,$b,$c,$d);
        return($pval,$val,"C",$direction);
    }
}
sub chi_square {
    my ($a,$b,$c,$d) = @_;
    #    print "$a\t$b\t$c\t$d\t";
    return(0) unless ($a && $b && $c && $d);
    #    $b = $b-$a;
    #    $c = $c-$a;
    #    $d = $d-$c-$b-$a;
    #    $d = $d - $c;

#        print "$a\t$b\t$c\t$d\t";
#    if ($a >= 5 && $b >= 5 && $c >= 5 && $d >= 5 ){
    my $tot = $a + $b + $c + $d;
    my $expa = ($a+$c)*($a+$b)/$tot;
    my $expb = ($b+$d)*($a+$b)/$tot;
    my $expc = ($a+$c)*($c+$d)/$tot;
    my $expd = ($b+$d)*($c+$d)/$tot;
    
    if ($expa >= 5 || $expb >= 5 || $expc >= 5 || $expd >= 5) {
        my $chival = &square(&abs($a-$expa)-0.5)/$expa +  &square(&abs($b-$expb)-0.5)/$expb + &square(&abs($c-$expc)-0.5)/$expc +  &square(&abs($d-$expd)-0.5)/$expd;
        
        my $pval = Statistics::Distributions::chisqrprob(1,&abs($chival));

        if ($a<$expa) {
            $chival = $chival * -1;
        }
        return ($pval,$chival);
    } else {
         #       print "\n";
        print STDERR "shouldn't get to this - should have been shunted into fisher exact test\n";
        return(1);
    }
}
sub fisher_exact {
    my ($x1,$x2,$y1,$y2) = @_;
    #Run fisher exact test in R                                                                                                                            

    $R->run("rm(list = ls())");
    $R->run("blah <- matrix(c(".$x1.",".$x2.",".$y1.",".$y2."),nrow=2)");
    $R->run("foo <- fisher.test(blah)");
    my $p_value_vs_bgd = $R->get('foo$p.value');
    my $val = "F";
    return($p_value_vs_bgd,$val);
}
