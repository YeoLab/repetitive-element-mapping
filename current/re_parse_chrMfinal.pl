use warnings;
use strict;

my %convert_strand = ("+" => "-", "-" => "+");
my %enst2chrstr;
my $chrM_genelistfi = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/genelists.chrM.wchr.txt";
open(CHR,$chrM_genelistfi);
for my $line (<CHR>) {
    chomp($line);
    my ($ensg,$genename,$enst,$chr,$str) = split(/\t/,$line);
    $enst2chrstr{$enst}{chr} = $chr;
    $enst2chrstr{$enst}{str} = $str;
}
close(CHR);


my $fi = $ARGV[0];

if ($fi =~ /\.gz$/) {
    open(FI,"gunzip -c $fi |") || die "couldn't open $fi\n";
} else {
    open(FI,$fi) || die "couldn't open $fi\n";
}
print STDERR "processing $fi\n";

my $read_counts=0;
my %reprocess_counts;
while (<FI>) {
    my $r1 = $_;
    next if ($r1 =~ /^\@/);
    my $r2 = <FI>;
    chomp($r1);
    chomp($r2);
        
    my @tmp_r1 = split(/\t/,$r1);
    my @tmp_r2 = split(/\t/,$r2);
    my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
    my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);
        
    unless ($tmp_r1[0] eq $tmp_r2[0]) {
	print STDERR "paired end mismatch error: $tmp_r1[0] $tmp_r2[0]\n";
    }

    my $r1sam_flag = $tmp_r1[1];
    my $r2sam_flag = $tmp_r2[1];
    next if ($r1sam_flag == 77 || $r1sam_flag == 141);

        # 77 = R1, unmapped                      
        # 141 = R2, unmapped                     

        # 99 = R1, mapped, fwd strand --- frag on rev strand             
        # 101 = R1 unmapped, R2 mapped rev strand -- frag on rev strand  
        # 73 = R1, mapped, fwd strand --- frag on rev strand             
        # 147 = R2, mapped, rev strand -- frag on rev strand             
        # 153 = R2 mapped (R1 unmapped), rev strand -- frag on rev strand
        # 133 = R2 unmapped, R1 mapped fwd strand -- frag on rev strand  

        # 83 = R1, mapped, rev strand --- frag on fwd strand             
        # 69 = R1 unmapped, R2 mapped fwd strand -- frag on fwd strand   
        # 89 = R1 mapped rev strand, R2 unmapped -- frag on fwd strand   
        # 163 = R2, mapped, fwd strand -- frag on fwd strand             
        # 137 = R2 mapped (R1 unmapped), fwd strand -- frag on fwd strand
        # 165 = R2 unmapped, R1 rev strand -- frag on fwd strand         

   
    my $frag_strand;
### This section is for only properly paired reads   
    if ($r1sam_flag == 99 || $r1sam_flag == 355) {
	$frag_strand = "-";
    } elsif ($r1sam_flag == 83 || $r1sam_flag == 339) {
	$frag_strand = "+";
    } elsif ($r1sam_flag == 147 || $r1sam_flag == 403) {
	($r1,$r2) =($r2,$r1);
	@tmp_r1 = split(/\t/,$r1);
	@tmp_r2 = split(/\t/,$r2);
	$frag_strand = "-";
    } elsif ($r1sam_flag == 163 || $r1sam_flag == 419) {
	($r1,$r2) =($r2,$r1);
	@tmp_r1 = split(/\t/,$r1);
	@tmp_r2 = split(/\t/,$r2);
	$frag_strand = "+";
    }  else {
	next;
	print STDERR "R1 strand error $r1sam_flag\n";
    }

###
    my @read_name = split(/\:/,$tmp_r1[0]);
    my $randommer = $read_name[0];
    my $r1_cigar = $tmp_r1[5];
    my $r2_cigar = $tmp_r2[5];
    my $r1_chr = $tmp_r1[2];
    my $r1_start = $tmp_r1[3];
    my $r2_chr = $tmp_r2[2];
    my $r2_start = $tmp_r2[3];

    
    my $enst_type = $tmp_r1[$#tmp_r1 - 1];
    if ($enst_type eq "UniqueGenomic" && $r1_chr eq "chrM") {
#	$reprocess_counts{$enst_type."||chrM"}{all}++;
#	$reprocess_counts{$enst_type."||chrM"}{$frag_strand}++;
	$reprocess_counts{"chrM_reprocess"}{$frag_strand."strand"}++;

	unless ($frag_strand) {
	    print STDERR "error - no absolute strand?\n$r1\n$r2\n";
	}
    } elsif ($r1_chr eq "chrM") {
	print STDERR "error? shouldn't hit this $r1_chr $enst_type\n";
    } elsif ($enst_type eq "RepFamily" && ($r1_chr =~ /^chrM\|\|/ || $r1_chr =~ /^antisense_chrM\|\|/)) {
	my ($r1_chronly,$r1_repelement) = split(/\|\|/,$r1_chr);
	
	my $relative_strand = "+";

	if ($r1_repelement =~ /^antisense\_(.+)$/) {
	    $r1_repelement = $1;
	    $relative_strand = "-";
	}

	if ($r1_repelement =~ /^(.+)\_DOUBLEMAP$/) {
	    $r1_repelement = $1;
	}
	my $absolute_strand = $enst2chrstr{$r1_repelement}{str};
	if ($relative_strand eq "-") {
	    $absolute_strand = $convert_strand{$absolute_strand};
	}

	unless ($absolute_strand) {
            print STDERR "error - no absolute strand? r1rep $r1_repelement relative $relative_strand r1chr $r1_chr\n$r1\n$r2\n";

	}

#	$reprocess_counts{$enst_type."||chrM"}{$absolute_strand}++;
#	$reprocess_counts{$enst_type."||".$r1_chronly}{all}++;
	$reprocess_counts{"chrM_reprocess"}{$absolute_strand."strand"}++;
    }
    $read_counts++;
}
close(FI);

#for my $element ("chrM_reprocess_+strand","chrM_reprocess_+strand") {

for my $k ("chrM_reprocess") {
#for my $k (keys %reprocess_counts) {
#    for my $kk (keys %{$reprocess_counts{$k}}) {
    for my $kk ("+strand","-strand") {
	$reprocess_counts{$k}{$kk} = 0 unless (exists $reprocess_counts{$k}{$kk});
	print "".$k."_".$kk."\t".$reprocess_counts{$k}{$kk}."\t".$read_counts."\t".($reprocess_counts{$k}{$kk} / $read_counts)."\n";
    }
}
