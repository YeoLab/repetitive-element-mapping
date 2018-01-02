use warnings;
use strict;

my $filist_fi = $ARGV[0];
#my $desired_element = "RNA18S||ENST00000606783.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/18S_only.fa";  
#my $desired_element = "RNA28S||ENST00000607521.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/28S_only.fa";    

my $desired_element = "rRNA_extra||NR_046235.1";
my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/rRNA_extra.fa";
#my $desired_element = "RNA5-8S||ENST00000474885.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/5.8S_only.fa";
my ($desfamily,$desenst) = split(/\|\|/,$desired_element);

my $working_dir = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20170304.txt";

#my $inputpairingfi = "Satellite_inputpairing_20170419.txt";
#my $working_dir = "/home/elvannostrand/scratch/Satellite/";

#my $inputpairingfi = "Enching_list.txt";
#my $working_dir = "/home/elvannostrand/scratch/Enching_RepElement/";

print STDERR "doing $desired_element\n";

my $output_file = $filist_fi.".".$desfamily.".csv";
open(OUT,">$output_file");
my $output_entfile = $filist_fi.".".$desfamily.".entropy.csv";
open(OUTENT,">$output_entfile");

my %clip2input;
my %clip2finalid;
open(IN,$inputpairingfi);
for my $line (<IN>) {
    chomp($line);
    $line =~ s/\r//g;

    my @tmp = split(/\s+/,$line);

    if ($tmp[3] eq "LNG" || $tmp[3] eq "TAG" || $tmp[3] eq "SingleRep" || $tmp[3] eq "OneRep") {
        $clip2input{$tmp[1]} = $tmp[2];


	if ($tmp[4] && $tmp[4] ne "x") {
	    my $final_id = $tmp[4]."-".$tmp[5]."-".$tmp[6];
	    $clip2finalid{$tmp[1]} = $final_id;
	    $clip2finalid{$tmp[2]} = $final_id;
	} else {
	    $clip2finalid{$tmp[1]} = $tmp[1];
#	    $clip2finalid{$tmp[1]} = $tmp[1]."-BRCA2-HepG2";
	}

    } else {
        $clip2input{$tmp[1]} = $tmp[3];
        $clip2input{$tmp[2]} = $tmp[3] unless ($tmp[2] eq "NoInput");

	if ($tmp[4] && $tmp[4] ne "x") {
	    my $final_id = $tmp[4]."-".$tmp[5]."-".$tmp[6];
	    $clip2finalid{$tmp[1]} = $final_id;
	    $clip2finalid{$tmp[2]} = $final_id;
	}

    }


}
close(IN);

my %uid2orig;
my @values;
my %fold_enr;
my %entropy;
open(FL,$filist_fi);
for my $line (<FL>) {
    chomp($line);
    next unless ($line);
    my $origline = $line;

    my @ids = split(/\_/,$line);
    my $rbp = $ids[2];
    if ($rbp eq "4000") {
	$rbp = $ids[3];
    }
#    my ($uid,$rep,$rbp) = split(/\_/,$line);
#    if ($rbp eq "NOLC1") {
#        $rbp = "BOP1";
#    }
#    $line = $uid."_".$rep."_".$rbp;
    
    $uid2orig{$line} = $line;

    my $input = $clip2input{$line};
#my $input = $uid."_INPUT_".$rbp;
    print STDERR "doing $line\n";
    my %file;
    $file{"CLIP"} = $working_dir.$line.".combined_w_uniquemap.rmDup.sam";
    $file{"INPUT"} = $working_dir.$input.".combined_w_uniquemap.rmDup.sam";

    my %usable_num;
    for my $type ("CLIP","INPUT") {
	my $outbam = $file{$type}.".tmp.".$desfamily.".bam";
	my $sorted = $outbam.".sorted.bam";
	my $tmp_sorting = $outbam.".tmp";
	my $mpileup_out = $sorted.".mpileup";

	my $readnum_fi = $file{$type}.".parsed";
	$usable_num{$type} = &read_readnumfi($readnum_fi);
	
	&parse_mpileup($mpileup_out,$line,$type);
	
    }
    
    for (my $i=1;$i<scalar(@values);$i++) {
	$values[$i]{$line}{"INPUT"}++;
	$values[$i]{$line}{"CLIP"} = 1 unless (exists $values[$i]{$line}{"CLIP"} && $values[$i]{$line}{"CLIP"} > 0);

#	$fold_enr{$line}[$i] = $values[$i]{$line}{"CLIP"};
	$fold_enr{$line}[$i] = sprintf("%.5f", ($values[$i]{$line}{"CLIP"} / $usable_num{"CLIP"}) / ($values[$i]{$line}{"INPUT"} / $usable_num{"INPUT"}));
	my $p_i = $values[$i]{$line}{"CLIP"} / $usable_num{"CLIP"};
	my $q_i = $values[$i]{$line}{"INPUT"} / $usable_num{"INPUT"};
	$entropy{$line}[$i] = sprintf( $p_i * log($p_i / $q_i) / log(2) );
    }
}
close(FL);

my %rbpnum;
for my $uid (keys %fold_enr) {
    my $finalid = $clip2finalid{$uid};
    my ($finaluid,$finalrbp,$finalcelltype) = split(/\-/,$finalid);
    unless ($finalrbp) {
	print STDERR "err coulnd't find rbp $uid $finalid\n";
    }
    my $current_rbpnum;
    if (exists $rbpnum{$finalrbp}) {
	$current_rbpnum = $rbpnum{$finalrbp};
    } else {
	$current_rbpnum = scalar(keys %rbpnum) + 1;
	$rbpnum{$finalrbp} = $current_rbpnum;
    }
    $fold_enr{$uid}[0] = "na";
    $entropy{$uid}[0] = "na";

#    print "".$uid2orig{$uid}."\t".$current_rbpnum."\t".join("\t",@{$fold_enr{$uid}})."\n";
    print OUT "".$finalid."\t".$current_rbpnum."\t".join("\t",@{$fold_enr{$uid}})."\n";
    print OUTENT "".$finalid."\t".$current_rbpnum."\t".join("\t",@{$entropy{$uid}})."\n";
}
close(OUT);
close(OUTENT);

sub read_readnumfi {
    my $file = shift;

    open(G,$file);
    for my $line (<G>) {
	chomp($line);
	if ($line =~ /^\#READINFO\tUsableReads\t(\d+)\t/) {
	    my $usable_reads = $1;
	    close(G);
	    return($usable_reads);
	}
    }
    close(G);
    print STDERR "shouldn't hit this - didn't find usable read number? $file\n";
}

sub parse_mpileup {
    my $fi = shift;
    my $uid = shift;
    my $type = shift;

    open(M,$fi) || die "no $fi\n";
    for my $line (<M>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	$values[$tmp[1]]{$uid}{$type}  = $tmp[3];
    }
    close(M);
}
	
