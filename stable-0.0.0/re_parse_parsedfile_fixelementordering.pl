use warnings;
use strict;

my %enst2gene;
my %convert_enst2type;
my $filelist_file = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/MASTER_filelist.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat";
my $priority_last = &read_in_filelists($filelist_file);
my $filelist_file2 = "/home/elvannostrand/data/clip/CLIPseq_analysis/RNA_type_analysis/ALLRepBase_elements.id_table.FULL";
$priority_last = &read_in_filelists2($filelist_file2,$priority_last);


my %counts_primary;
my %counts;
#my $fi = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/539_02_NOL12.combined_w_uniquemap.rmDup.sam.parsed";
#my $fi = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/539_01_NOL12.combined_w_uniquemap.rmDup.sam.parsed";
my $fi = "/home/elvannostrand/scratch/ENCODE_20170429_newannotations_FINAL/539_INPUT_NOL12.combined_w_uniquemap.rmDup.sam.parsed";
my $outfi = $fi.".reparse_primary.txt";
open(OUT,">$outfi");
open(F,$fi);
for my $line (<F>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    if ($tmp[0] eq "ELEMENT") {
	
	my $family = $tmp[1];
	my $read_n = $tmp[2];
	my $read_rpr = $tmp[3];
	my $all_transcripts = $tmp[4];
	my $all_names = $tmp[5];

	if ($family =~ /^unique/) {
	    print OUT "$line\tUnique_gencode\n";
	    next;
	}
	
	my ($families_all,$transcripts_all) = split(/\|\|/,$all_transcripts);
	my @transcripts = split(/\|/,$transcripts_all);
	my @transcripts_sorted = sort {$a cmp $b } @transcripts;
	my $new_transcripts_all = join("|",@transcripts_sorted);
	
	my @family_split = split(/\|/,$family);
	my @family_sorted = sort { $a cmp $b } @family_split;
	my $new_family = join("|",@family_sorted);

	if (scalar(@family_sorted) > 1) {
	    print OUT "$line\tMultiple_families\n";
	    next;
	}

	my $new_hashval_all = $new_family."||".$new_transcripts_all;
	if (exists $counts{$new_hashval_all}{primary_family} && $counts{$new_hashval_all}{primary_family} ne $new_family) {
	    print STDERR "error family doesn't match? ".$counts{$new_hashval_all}{primary_family}." $new_family $line\n";
	}

	my %antisense_flags;
	my ($min_priority_n,$min_priority_enst);
	for my $enst (@transcripts_sorted) {
	    my $antisense_flag = "";
	    $enst =~ s/\(//g;
	    $enst =~ s/\)//g;

	    if ($enst =~ /^(.+)\_DOUBLEMAP/ || $enst =~ /^(.+)\_uniquegenomic/) {
		$enst = $1;
	    }
	    if ($enst =~ /^(.+)\_spliced/) { 
		$enst = $1;
	    }
	    unless (exists $convert_enst2type{$enst}) {
		print STDERR "missing $enst\n";
		next;
#		$priority_n = 0;
#		$type_label = ;
	    }
	    my ($type_label,$priority_n) = split(/\:/,$convert_enst2type{$enst});
	    
	    if ($min_priority_enst) {
		if ($priority_n < $min_priority_n) {
		    $min_priority_enst = $enst;
		    $min_priority_n = $priority_n;
		}
	    } else {
		$min_priority_enst = $enst;
		$min_priority_n = $priority_n;
	    }
	}
	

	$counts{$new_hashval_all}{primary_family} = $new_family;
	$counts{$new_hashval_all}{read_n} += $read_n;
	$counts{$new_hashval_all}{read_rpr} += $read_rpr;
	$counts{$new_hashval_all}{main_enst} = $min_priority_enst;
	
	$counts_primary{$min_priority_enst}{read_n} += $read_n;
	$counts_primary{$min_priority_enst}{read_rpr} += $read_rpr;
	$counts_primary{$min_priority_enst}{primary_family} = $new_family;
	
    }  else {
	print OUT "$line\n";
    }

}
close(F);



for my $new_hashval_all (keys %counts) {
    print OUT "ELEMENT\t".$counts{$new_hashval_all}{primary_family}."\t".$counts{$new_hashval_all}{read_n}."\t".$counts{$new_hashval_all}{read_rpr}."\t".$new_hashval_all."\t".$counts{$new_hashval_all}{main_enst}."|".$enst2gene{$counts{$new_hashval_all}{main_enst}}."\n";
}

for my $primary_enst (keys %counts_primary) {
    print OUT "PRIMARYELEMENT\t".$counts_primary{$primary_enst}{primary_family}."\t".$counts_primary{$primary_enst}{read_n}."\t".$counts_primary{$primary_enst}{read_rpr}."\t".$primary_enst."|".$enst2gene{$primary_enst}."\n";
}




sub read_in_filelists {
    my $fi = shift;
    my $priority_n = 0;
    open(F,$fi);
    for my $line (<F>) {
        chomp($line);
        my ($allenst,$allensg,$gid,$type_label,$typefile) = split(/\t/,$line);
#       my ($allensg,$gid,$allenst) = split(/\t/,$line);
#        $type_label =~ s/\_$//;
        unless ($allenst) {
            print STDERR "error missing enst $line $fi\n";
        }
        my @ensts = split(/\|/,$allenst);
        for my $enst (@ensts) {
            $enst2gene{$enst} = $gid;
            $enst2gene{"antisense_".$enst} = "antisense_".$gid;
            $convert_enst2type{$enst} = $type_label.":".$priority_n;
            $convert_enst2type{"antisense_".$enst} = "antisense_".$type_label.":".$priority_n;
            $priority_n++;
        }
    }
    close(F);
    return($priority_n);
}



sub read_in_filelists2 {
    my $fi = shift;
    my $priority_n = shift;
    open(F,$fi);
    for my $line (<F>) {
        chomp($line);
        my ($allenst,$allensg,$gid,$type_label,$typefile) = split(/\t/,$line);
#       my ($allensg,$gid,$allenst) = split(/\t/,$line);
        unless ($allenst) {
            print STDERR "error missing enst $line $fi\n";
        }
        my @ensts = split(/\|/,$allenst);
        $gid =~ s/\?$//;
        $gid =~ s/\_$//;
        $type_label =~ s/\?$//;
        for my $enst (@ensts) {
	    $enst =~ s/\(//;
	    $enst =~ s/\)//;

            $enst2gene{$enst} = $gid;
            $enst2gene{"antisense_".$enst} = "antisense_".$gid;
#            $convert_enst2type{$enst} = $type_label.":".$priority_n;
            $convert_enst2type{"antisense_".$enst} = "antisense_".$type_label.":".$priority_n;
            $convert_enst2type{$enst} = $type_label.":".$priority_n;

            $priority_n++;
        }
    }
    close(F);
}
