use warnings;
use strict;


my $sam_fi = $ARGV[0];
open(BOWTIE,$sam_fi);

my $rmdup_fi = $sam_fi.".rmDup.sam";
open(OUT,">$rmdup_fi");

my ($all_count,$duplicate_count,$unique_count) = (0,0,0);
my %fragment_hash;

while (<BOWTIE>) {
    my $r1 = $_;
    chomp($r1);
    if ($r1 =~ /^\@/) {
	print OUT "$r1\n";
	next;
    }

    my $r2 = <BOWTIE>;
    chomp($r2);

    my @tmp_r1 = split(/\t/,$r1);
    my @tmp_r2 = split(/\t/,$r2);
    my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
    my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);

    unless ($r1name eq $r2name) {
	print STDERR "paired end mismatch error: $sam_fi r1 $tmp_r1[0] r2 $tmp_r2[0]\n";
    }

    my $r1sam_flag = $tmp_r1[1];
    my $r2sam_flag = $tmp_r2[1];
    unless ($r1sam_flag) {
	print STDERR "error $r1 $r2\n";
    }
    next if ($r1sam_flag == 77 || $r1sam_flag == 141);

    my $frag_strand;
### This section is for only properly paired reads                                                                                                                                                                       
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




    my @read_name = split(/\:/,$tmp_r1[0]);
    my $randommer = $read_name[0];

    my $r1_cigar = $tmp_r1[5];
    my $r2_cigar = $tmp_r2[5];

      # 165 = R2 unmapped, R1 rev strand -- frag on fwd strand


    my $r1_chr = $tmp_r1[2];
    my $r1_start = $tmp_r1[3];

    my $r2_chr = $tmp_r2[2];
    my $r2_start = $tmp_r2[3];


    my $mismatch_flags = $tmp_r1[14];

    my %read_regions;
    @{$read_regions{"R1"}} = &parse_cigar_string($r1_start,$r1_cigar,$r1_chr,$frag_strand);
    @{$read_regions{"R2"}} = &parse_cigar_string($r2_start,$r2_cigar,$r2_chr,$frag_strand);

    my $hashing_value;
    if ($frag_strand eq "+") {
        my $r1_firstregion = $read_regions{"R1"}[scalar(@{$read_regions{"R1"}})-1];
        my ($r1_firstchr,$r1_firststr,$r1_firstpos) = split(/\:/,$r1_firstregion);
        my ($r1_firststart,$r1_firststop) = split(/\-/,$r1_firstpos);
        my $r2_firstregion = $read_regions{"R2"}[0];
        my ($r2_firstchr,$r2_firststr,$r2_firstpos) = split(/\:/,$r2_firstregion);
        my ($r2_firststart,$r2_firststop) = split(/\-/,$r2_firstpos);
        $hashing_value = $r2_firstchr.":".$r2_firststr.":".$r2_firststart."\t".$r1_firstchr.":".$r1_firststr.":".$r1_firststop;
    } elsif ($frag_strand eq "-") {
        my $r1_firstregion = $read_regions{"R1"}[0];
        my ($r1_firstchr,$r1_firststr,$r1_firstpos) = split(/\:/,$r1_firstregion);
        my ($r1_firststart,$r1_firststop) = split(/\-/,$r1_firstpos);
        my $r2_firstregion = $read_regions{"R2"}[scalar(@{$read_regions{"R2"}})-1];
        my ($r2_firstchr,$r2_firststr,$r2_firstpos) = split(/\:/,$r2_firstregion);
        my ($r2_firststart,$r2_firststop) = split(/\-/,$r2_firstpos);
        $hashing_value = $r1_firstchr.":".$r1_firststr.":".$r1_firststart."\t".$r2_firstchr.":".$r2_firststr.":".$r2_firststop;

    } else {
        print STDERR "strand error $frag_strand $r1 $r2\n";
    }

    my $full_frag_position = join("|",@{$read_regions{"R1"}})."\t".join("|",@{$read_regions{"R2"}});

    my $debug_flag = 0;


    my %full_fragment;
    $all_count++;
    if (exists $fragment_hash{$r1_chr."|".$frag_strand}{$hashing_value.":".$randommer}) {
        $duplicate_count++;
    } else {
        $fragment_hash{$r1_chr."|".$frag_strand}{$hashing_value.":".$randommer} = 1;
	print OUT "".$r1."\n".$r2."\n";

        $unique_count++;
    }
}

close(OUT);
close(BOWTIE);

#print STDERR "$sam_fi\nunique $unique_count\nduplicate $duplicate_count\nall $all_count\n";


sub parse_cigar_string {
    my $region_start_pos = shift;
    my $flags = shift;
    my $chr = shift;
    my $strand = shift;

    my $current_pos = $region_start_pos;
    my @regions;

    while ($flags =~ /(\d+)([A-Z])/g) {

        if ($2 eq "N") {
            #read has intron of N bases at location                  

            push @regions,$chr.":".$strand.":".$region_start_pos."-".($current_pos-1);

            $current_pos += $1;
            $region_start_pos = $current_pos;
        } elsif ($2 eq "M") {
            #read and genome match         
            $current_pos += $1;
        } elsif ($2 eq "S") {
            #beginning of read is soft-clipped; mapped pos is actually start pos of mapping not start of read            
        } elsif ($2 eq "I") {
            #read has insertion relative to genome; doesn't change genome position             
        } elsif ($2 eq "D") {
#           push @read_regions,$chr.":".$current_pos."-".($current_pos+=$1);                   
            $current_pos += $1;
            #read has deletion relative to genome; genome position has to increase             
        } else {
            print STDERR "flag $1 $2 $flags\n";

        }
    }
    push @regions,$chr.":".$strand.":".$region_start_pos."-".($current_pos-1);

    return(@regions);
}
