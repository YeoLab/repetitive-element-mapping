use warnings;
use strict;


my %uid2gabeid;
my $idconversion_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/scripts/inline_processing/RepElement_InputPairing_20170304.txt";
open(ID,$idconversion_fi);
for my $line (<ID>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    next unless ($tmp[4]);
    my $uid = $tmp[4];
    
    if ($tmp[0] eq "CORRECT" || $tmp[0] eq "manual") {
	$uid2gabeid{$uid}{"01"} = $tmp[1];
	$uid2gabeid{$uid}{"02"} = $tmp[2];
    }
    
}
close(ID);

my %submitted_uid;
my $submitted_manifest = "/home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/ALLDATASETS_submittedonly.txt";
open(SUB,$submitted_manifest);
for my $line (<SUB>) {
    chomp($line);
    my @tmp = split(/\s+/,$line);
    next if ($tmp[0] eq "uID");
    my $uid = $tmp[0];
    my $id = $tmp[1]."-".$tmp[2];
    
    $submitted_uid{$uid} = $id;


}
close(SUB);


my %data_hash;
my $file = "20170228_parsedoutput.ALLdatasets.csv.nopipes.csv";
my $output_file = $file.".top_categories_per_dataset.csv";
#my $output_grouped = $output_file.".grouped";
open(OUT,">$output_file");
#open(GROUP,">$output_grouped");

open(F,$file);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    
    my $uid;
    my $rep;
    if ($tmp[1] =~ /^(.+)\_(01)\_(.+)$/ || $tmp[1] =~ /^(.+)\_(02)\_(.+)$/) {
	$uid = $1."_".$3;
	$rep = $2;
    } else {
	next if ($tmp[0] eq "element");
	print STDERR "couldn't parse $tmp[1]\n";
    }
    
    $data_hash{$tmp[1]}{$tmp[0]} = $tmp[7];
#    $data_hash{$uid}{$rep}{$tmp[0]} = $tmp[7];
    
}
close(F);

for my $uid (keys %submitted_uid) {
#for my $uid (keys %data_hash) {
    my $rep1_id = $uid2gabeid{$uid}{"01"};
    my $rep2_id = $uid2gabeid{$uid}{"02"};

    my @rep1_sorted = sort {$data_hash{$rep1_id}{$b} <=> $data_hash{$rep1_id}{$a}} keys %{$data_hash{$rep1_id}};
    my @rep2_sorted = sort {$data_hash{$rep2_id}{$b} <=> $data_hash{$rep2_id}{$a}} keys %{$data_hash{$rep2_id}};
#    my @rep1_sorted = sort {$data_hash{$uid}{"01"}{$b} <=> $data_hash{$uid}{"01"}{$a}} keys %{$data_hash{$uid}{"01"}};
#    my @rep2_sorted = sort {$data_hash{$uid}{"02"}{$b} <=> $data_hash{$uid}{"02"}{$a}} keys %{$data_hash{$uid}{"02"}};
    
    
    my ($out1,$out2,$out_merge);
    if ($rep1_sorted[0] eq $rep2_sorted[0]) {
	$out1 = "$uid\t$submitted_uid{$uid}\t01\tMATCHED";
	$out2 = "$uid\t$submitted_uid{$uid}\t02\tMATCHED";
	$out_merge = "$uid\t$submitted_uid{$uid}\tmerged";
    } else {
	$out1 = "$uid\t$submitted_uid{$uid}\t01\tNOT";
	$out2 = "$uid\t$submitted_uid{$uid}\t02\tNOT";
	$out_merge = "$uid\t$submitted_uid{$uid}\tmerged";
    }

    for my $element (keys %{$data_hash{$rep1_id}},keys %{$data_hash{$rep2_id}}) {
	$data_hash{$rep1_id}{$element} = 0 unless (exists $data_hash{$rep1_id}{$element});
	$data_hash{$rep2_id}{$element} = 0 unless (exists $data_hash{$rep2_id}{$element});

	my $geometric_mean = log(sqrt( (2 ** $data_hash{$rep1_id}{$element}) * (2 ** $data_hash{$rep2_id}{$element}) ))/log(2);
	$data_hash{$uid}{"merge"}{$element} = $geometric_mean;
    }
    my @sorted_geommean = sort {$data_hash{$uid}{"merge"}{$b} <=> $data_hash{$uid}{"merge"}{$a}} keys %{$data_hash{$uid}{"merge"}};

    my %to_keep;
    $to_keep{$sorted_geommean[0]} = 1;
    my $top_g = $data_hash{$uid}{"merge"}{$sorted_geommean[0]};
    for (my $i=1;$i<@sorted_geommean;$i++) {
	if ($data_hash{$uid}{"merge"}{$sorted_geommean[$i]} >= .5 * $top_g) {
	    $to_keep{$sorted_geommean[$i]} = 1;
	}
    }
    $out_merge .= "\t".join("|",keys %to_keep);
	
    for my $i (0..4) {
	$out1 .= "\t".$rep1_sorted[$i]."\t".$data_hash{$rep1_id}{$rep1_sorted[$i]} if (exists $rep1_sorted[$i]);
	$out2 .= "\t".$rep2_sorted[$i]."\t".$data_hash{$rep2_id}{$rep2_sorted[$i]} if (exists $rep2_sorted[$i]);
	$out_merge .= "\t".$sorted_geommean[$i]."\t".$data_hash{$uid}{"merge"}{$sorted_geommean[$i]} if (exists $sorted_geommean[$i]);
    }
    print OUT "$out1\n$out2\n$out_merge\n";

}


close(OUT);


