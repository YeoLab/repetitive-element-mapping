use warnings;
use strict;

my %hash;
my $fi1 = $ARGV[0];
open(F,$fi1);
while (<F>) {
    chomp($_);
    my $r1 = $_;    
    if ($r1 =~ /^\@/) {
	next;
    }

    my $r2 = <F>;
    chomp($r2);

    my @tmp_r1 = split(/\t/,$r1);
    my @tmp_r2 = split(/\t/,$r2);
    my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
    my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);

    unless ($r1name eq $r2name) {
	print STDERR "paired end mismatch error: r1 $tmp_r1[0] r2 $tmp_r2[0]\n";
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
    } elsif ($r1sam_flag == 163 || $r1sam_flag == 419) {
	$frag_strand = "+";
    }  else {
	next;
	print STDERR "R1 strand error $r1sam_flag\n";
    }
###
    
    $hash{$r1name} = $r1."\n".$r2;
}
close(F);

my $matched=0;

my $fi2 = $ARGV[1];
open(G,$fi2);
while (<G>) {
    chomp($_);
    my $r1 = $_;
    if ($r1 =~ /^\@/) {
        next;
    }

    my $r2 = <G>;
    chomp($r2);

    my @tmp_r1 = split(/\t/,$r1);
    my @tmp_r2 = split(/\t/,$r2);
    my ($r1name,$r1bc) = split(/\s+/,$tmp_r1[0]);
    my ($r2name,$r2bc) = split(/\s+/,$tmp_r2[0]);

    unless ($r1name eq $r2name) {
	print STDERR "paired end mismatch error: r1 $tmp_r1[0] r2 $tmp_r2[0]\n";
    }

    if (exists $hash{$r1name}) {
	if ($hash{$r1name} eq $r1."\n".$r2) {
	    $matched++;
	    delete $hash{$r1name};
	} else {
	    print "mismatch between fi1 and fi2:\n$hash{$r1name}\n".$r1."\n".$r2."\n";
	}
    } else {
	print "doesn't exist in fi1 ".$r1."\n".$r2."\n";
    }
}
close(G);

for my $k (keys %hash) {
    print "doesn't exist in fi2 ".$hash{$k}."\n";
}
print "matched = $matched\n";
