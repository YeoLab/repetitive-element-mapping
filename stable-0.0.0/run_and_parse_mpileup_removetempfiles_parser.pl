use warnings;
use strict;

my $desired_element = "RNU1||ENST00000383925.1";
my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/RNU1-1.fa";

#my $desired_element = "RNA18S||ENST00000606783.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/18S_only.fa";  
#my $desired_element = "RNA28S||ENST00000607521.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/28S_only.fa";    

#my $desired_element = "rRNA_extra||NR_046235.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/rRNA_extra.fa";
#my $desired_element = "RNA5-8S||ENST00000474885.1";
#my $desired_fasta = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/5.8S_only.fa";
my ($desfamily,$desenst) = split(/\|\|/,$desired_element);

my $working_dir = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/";
#my $working_dir = "/home/elvannostrand/scratch/Satellite_test4_20170428/";
#my $working_dir = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/";
#my $working_dir = "/home/elvannostrand/scratch/ALLCLIP_repmapping_20170201/";
my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20170304.txt";
#my $inputpairingfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/Satellite_inputpairing_20170419.txt";
#my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/ALLDATASETS_submittedonly.txt";
#my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_notyetsubmitted_20170325/NIPBL.txt";
#my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/Collaborations/Satellite/Satellite.txt2";
my $manifest_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_notyetsubmitted_20170325/FUS_list.txt";

my %desired_uids;
open(M,$manifest_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    my $uid = $tmp[0];
    next if ($uid eq "uID");
    $desired_uids{$uid} = 1;
}
close(M);

#my $inputpairingfi = "Satellite_inputpairing_20170419.txt";
#my $working_dir = "/home/elvannostrand/scratch/Satellite/";

#my $inputpairingfi = "Enching_list.txt";
#my $working_dir = "/home/elvannostrand/scratch/Enching_RepElement/";

print STDERR "doing $desired_element\n";
my $filist_fi = $manifest_fi;
my $output_file = $filist_fi.".".$desfamily.".csv";
open(OUT,">$output_file");
my $output_entfile = $filist_fi.".".$desfamily.".entropy.csv";
open(OUTENT,">$output_entfile");

my @file_list;
my %clip2input;
my %clip2finalid;
open(IN,$inputpairingfi);
for my $line (<IN>) {
    chomp($line);
    $line =~ s/\r//g;

    my @tmp = split(/\s+/,$line);

    next unless ($tmp[4] && exists $desired_uids{$tmp[4]});
    
    push @file_list,$tmp[1];
    push @file_list,$tmp[2];
    $clip2input{$tmp[1]} = $tmp[3];
    $clip2input{$tmp[2]} = $tmp[3];
    $clip2finalid{$tmp[1]} = $tmp[4]."-".$tmp[5]."-".$tmp[6]."_01";
    $clip2finalid{$tmp[2]} = $tmp[4]."-".$tmp[5]."-".$tmp[6]."_02";
#    last if ($tmp[4] eq "209");
}
close(IN);

my %uid2orig;
my @values;
my %fold_enr;
my %entropy;

for my $line (@file_list) {
    $uid2orig{$line} = $line;

    my $input = $clip2input{$line};

    print STDERR "doing $line\n";
    my %file;
    $file{"CLIP"} = $working_dir.$line.".combined_w_uniquemap.rmDup.sam";
    $file{"INPUT"} = $working_dir.$input.".combined_w_uniquemap.rmDup.sam";

    my %usable_num;
    for my $type ("CLIP","INPUT") {
	my $outbam = $file{$type}.".gz.tmp.".$desfamily.".bam";
	my $sorted = $outbam.".sorted.bam";
	my $tmp_sorting = $outbam.".tmp";
	my $mpileup_out = $sorted.".mpileup.short";

	my $readnum_fi = $file{$type}.".parsed";
	$usable_num{$type} = &read_readnumfi($readnum_fi);
	
	&parse_mpileup_short($mpileup_out,$line,$type);
	
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


my %rbpnum;
for my $uid (@file_list) {
    my $finalid = $clip2finalid{$uid};
    my ($finaluid,$finalrbp,$finalcelltype2) = split(/\-/,$finalid);
    my ($finalcelltype,$finalrep) = split(/\_/,$finalcelltype2);
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

    open(G,$file) || die "no $file\n";
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

sub parse_mpileup_short {
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
	
