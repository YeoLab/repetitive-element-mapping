# repetitive-element-mapping
pipeline for mapping repetitive elements

# Requirements:
- Bowtie2=2.2.6
- perl=5.22.0 or perl=5.10.1 (perl 5.18+ introduces some randomness in
iterating over hashes, which may yield slightly different results). This is
due to the way this pipeline accesses reads via hashes, and determines ties
in quality by selecting the first read it encounters at the de-duplication
step. We have tried to mitigate this randomness by iterating over sorted
hash keys in 0.0.2. See NOTES.
- cwlref-runner=1.0

(Or run the ```source create_environment.sh``` script)

# Installation:
- For Yeo Lab: ```module load ecliprepmap```
- For all others:
    - ensure that the ```bin/```, ```bin/perl/```, and ```wf/```, and ```cwl/```
    are correctly in your path.

# Methods:
- (maprep.cwl - parse_bowtie2_output_realtime_includemultifamily) Runs bowtie2 using the following commands: ```bowtie2 -q --sensitive -a -p 3 --no-mixed --reorder -x $bowtie_db -1 $fastq_file1 -2 $fastq_file2 2> $bowtie_out```, where:
    - $fastq_file1 is read 1 of a trimmed CLIPSEQ expt
    - $fastq_file2 is read 2 of a trimmed CLIPSEQ expt
    - $bowtie_db is a bowtie database created using a manually curated set of repeat elements from RepBase or other.
    - For every proper read pair:
        - Mismatch score equals the sum of mismatch scores (as determined by the AS: field) for both reads
        - Reads are compared on the basis of their mismatch scores and read quality, keeping the best or prioritize if equal.
            - priority is based on input repeat database (e.g. primary transcripts are prioritized over pseudogenes). It is based on the ordering of the repeat database file (MASTER_filelist.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat)
        - Since one read may map to multiple elements, we must compare each element
            - If read maps to multiple elements of the same family with equal score:
                - Read mapping information (fastq line) is kept for the element with the highest priority.
                - Name of all elements is kept
            - Note: all reads that map to multiple families (not used for downstream analysis).
        - Create a SAM-like file
- (deduplicate.cwl - duplicate_removal_inline_paired_count_region_other_reads) Merge repeat analysis with unique genomic mapping and remove PCR duplicates
    - Use the randomer UMIs to remove duplicates based on quality
        - To save memory, both the SAM-like file and the uniquely mapped BAM file
        are split based on the first two nucleotides of the randomner umi. De-duplication
        is performed serially for AA, AC, AG, AT .. NN.
    - Between a read mapping to both unique genomic and repeat family, if
    unique mapping to genome is more than 2 mismatches per
    read = 2 * 2 * 6 alignment score better than to repeat element,
    throw out repeat element and use genome mapping. Otherwise keep the repeat
    element mapped read.

### Determining the most correct assignment among elements within one family:
- Between longer and shorter transcripts: keep the shorter one
- If a read maps to two places on the same transcript, keep the first
- Treat "rRNA extra hash" different. See: ```parse_bowtie2_output_realtime_includemultifamily.pl``` $rRNA_extra_hash
    - This is different because the rRNA precusor transcript contains 18S, 28S, 5.8S transcripts which are treated as separate families.
# Outputs:
.parsed file: tabbed file containing 4 comment lines and 4 or 6 columns:
- #READINFO (AllReads): total number of reads mapped
- #READINFO (UsableReads): total number of de-duplicated mapped reads
- #READINFO (GenomicReads): total number of uniquely mapped genomic reads
- #READINFO (RepFamilyReads): total number of repetitive family mapped reads
- total_or_element: this is either TOTAL or ELEMENT
    - Total: family (you can create this by summing all of the elements within this family)
        - total
        - family name
        - number of reads
        - reads per million
    - Element: element
        - element
        - family name (eg. family1). Can have multiple pipe-delimited families
        - number of reads
        - reads per million
        - family1||transcript1|transcript2|transcript3 (transcript 1/2/3 are all members of family1)
        - gene name of the transcripts

# Notes:
- Elements (second column) that contain pipes ```|``` are multifamily mapped
which means that the read is ambiguously mapped to each element. These can generally
be excluded from downstream analysis.
- You can sort by information content and log2 fold enrichment.
Information content is calculated as: log2(clip_rpr/input_rpr).
- High fold changes and information content (no guidelines yet to this cutoff)
typically indicate elements enriched in the experiment.