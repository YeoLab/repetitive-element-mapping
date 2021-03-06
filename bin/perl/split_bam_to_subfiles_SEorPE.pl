#!/usr/bin/env perl

use warnings;
use strict;

my $sam_fi = $ARGV[0];
my @sam_fi_split = split(/\//,$sam_fi);
my $sam_fi_short = $sam_fi_split[$#sam_fi_split];

my $se_or_pe_flag = $ARGV[1];
if ($se_or_pe_flag eq "PE"){
    print STDERR "splitting in paired-end mode\n";
} elsif ($se_or_pe_flag eq "SE") {
    print STDERR "splitting in single-end mode\n";
} else {
    print STDERR "fatal error - SE or PE not defined\n";
    exit;
}

my %filehandles;
for my $b1 ("A","C","G","T","N") {
    for my $b2 ("A","C","G","T","N") {
    
	my $outfi = $b1.$b2.".".$sam_fi_short.".tmp";
	open(my $fh, '>', $outfi);
	$filehandles{$b1.$b2} = $fh;
    }
}



if ($sam_fi =~ /\.sam$/) {
    open(F,$sam_fi);
} elsif ($sam_fi =~ /\.bam$/) {
    open(F,"samtools view -h $sam_fi |") || die "no $sam_fi\n";
} else {
    print STDERR "weird - $sam_fi not either sam or bam file format - exit\n";
    exit;
}

while (<F>) {
    my $r1 = $_;

    next if ($r1 =~ /^\@/);

    chomp($r1);
    my @tmp_r1 = split(/\t/,$r1);
    my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
    my $r1sam_flag = $tmp_r1[1];

    my $r2 = <F> if ($se_or_pe_flag eq "PE");
    chomp($r2) if ($se_or_pe_flag eq "PE");
    my @tmp_r2 = split(/\t/,$r2) if ($se_or_pe_flag eq "PE");
    my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]) if ($se_or_pe_flag eq "PE");
    my $r2sam_flag = $tmp_r2[1] if ($se_or_pe_flag eq "PE");

    if ($se_or_pe_flag eq "PE") {
	unless ($r1name eq $r2name) {
	    print STDERR "paired end mismatch error: $sam_fi r1 $tmp_r1[0] r2 $tmp_r2[0]\n";
	}
	unless ($r1sam_flag) {
	    print STDERR "error $r1 $r2\n";
	}
    } 

    

    if ($se_or_pe_flag eq "PE") {
	next if ($r1sam_flag == 77 || $r1sam_flag == 141);
    } elsif ($se_or_pe_flag eq "SE") {
	next if ($r1sam_flag == 4);
    }
    my $frag_strand;
### This section is for only properly paired reads                                                                                                                
    if ($se_or_pe_flag eq "PE"){
	if ($r1sam_flag == 99 || $r1sam_flag == 355) {
	    $frag_strand = "-";
	} elsif ($r1sam_flag == 83 || $r1sam_flag == 339) {
	    $frag_strand = "+";
	} elsif ($r1sam_flag == 147 || $r1sam_flag == 403) {
	    $frag_strand = "-";
	    @tmp_r1 = split(/\t/,$r2);
	    @tmp_r2 = split(/\t/,$r1);
	} elsif ($r1sam_flag == 163 || $r1sam_flag == 419) {
	    $frag_strand = "+";
	    @tmp_r1 = split(/\t/,$r2);
	    @tmp_r2 = split(/\t/,$r1);
	}  else {
	    next;
	    print STDERR "R1 strand error $r1sam_flag\n";
	}
	my @read_name = split(/\:/,$tmp_r1[0]);
	my $randommer = $read_name[0];
	my $first2rand = substr($randommer,0,2);
	print { $filehandles{$first2rand} } "".$r1."\n".$r2."\n";



    } elsif ($se_or_pe_flag eq "SE") {
	if ($r1sam_flag == 16 || $r1sam_flag == 272) {
	    $frag_strand = "-";
	} elsif ($r1sam_flag eq "0" || $r1sam_flag == 256) {
	    $frag_strand = "+";
	}  else {
	    next;
	    print STDERR "R1 strand error $r1sam_flag\n";
	}

	my @read_name = split(/\_/,$tmp_r1[0]);
	my $randommer = pop(@read_name);
	my $first2rand = substr($randommer,0,2);
	print { $filehandles{$first2rand} } "".$r1."\n";

    }
###                                          

       # 77 = R1, unmapped                                                                      
    # 141 = R2, unmapped                                                                   
    # 99 = R1, mapped, fwd strand --- frag on rev strand                         ->    355 = not primary
    # 147 = R2, mapped, rev strand -- frag on rev strand                         ->    403 = not primary

    # 101 = R1 unmapped, R2 mapped rev strand -- frag on rev strand                                               
    # 73 = R1, mapped, fwd strand --- frag on rev strand                                                
    # 153 = R2 mapped (R1 unmapped), rev strand -- frag on rev strand                                   
    # 133 = R2 unmapped, R1 mapped fwd strand -- frag on rev strand                                                             

    # 83 = R1, mapped, rev strand --- frag on fwd strand                     ->    339 = not primary  
    # 163 = R2, mapped, fwd strand -- frag on fwd strand                     ->    419 = not primary  

    # 69 = R1 unmapped, R2 mapped fwd strand -- frag on fwd strand                                     
    # 89 = R1 mapped rev strand, R2 unmapped -- frag on fwd strand                          
    # 137 = R2 mapped (R1 unmapped), fwd strand -- frag on fwd strand                                                 
    # 165 = R2 unmapped, R1 rev strand -- frag on fwd strand                               



}
close(F);
