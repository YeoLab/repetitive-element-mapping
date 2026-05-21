#### System Prompt: Generate reference files for mice species mm10 and mm39

**Role:** Act as an expert bioinformatician.

### **Skills**
Use the skills provided in:
- /tscc/nfs/home/bay001/projects/codebase/repetitive-element-mapping/.claude/skills/

##### Context
*   **Tech Stack:** python, perl
*   **Architecture:** workflows written in Common Workflow Language (CWL)
*   **Background:** We are generating reference data for additional species

##### Scope
*   You may **read** from all files and softlinks in this folder, but may **only write** to examples/inputs/mm39 or examples/inputs/mm10.
    - **IFF** perl scripts contain hardcoded values (See Step 5) that prevent its acceptance of new annotations, you may modify them.
*   Review the entire prompt before starting work. You are allowed to push back on logic if you disagree, but be very clear and always ask permission.
*   Use the following modules to load the correct python, cwl and perl environments:
    - ```module load ecliprepmap/1.0.0```
*   Generate the following files for assemblies mm10 and mm39 (source examples/inputs/{assembly}):
    - gencode.{version}.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
    - bowtie2_index
    - UniqueGenomicElements.{assembly}.bed
    - MASTER_FILELIST.{date}.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.tsv

##### Constraints & Framework Rules
*   **Use established sources:** Original source files should all be provided inside examples/inputs/{assembly}/downloaded with one exception: Elements that start with "NR_" may exist in the refseq database. For mm10 and mm39-generated references, you must include the "NR_046233.2.fasta" sequence, which I have provided.
*   **Use existing generated refdata for hg38:** Pipeline reference files are provided inside examples/inputs/hg38/
*   **Always verify format is identical to reference data generated for hg38:** 

##### Logic & Implementation Steps 
**Step 1: generate parsed_ucsc_tableformat**
*   Write a python script that transforms the GTF (eg. gencode.v33.chr_patch_hapl_scaff.annotation.gtf) into a "parsed_ucsc_tableformat" file (eg. gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat). Skip comment lines "#".
*   Each parsed_ucsc_tableformat line should represent a unique transcript_id of feature (GTF column 3) "transcript"
*   parsed_ucsc_tableformat "exonCount" should equal the number of "exon" features for the transcript. This should equal the number of comma-delimited "exonStart" and "exonEnd" values which represent start and end (GTF columns 4 and 5) positions, respectively.
*   Test this python script by reproducing a new parsed_ucsc_tableformat from gencode.v33.chr_patch_hapl_scaff.annotation.gtf, which should be identical to the original source gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
  
**Step 2: generate bowtie_index**
*   Identify the sequences used to generate the bowtie_index from examples/inputs/hg38/bowtie2_index/MASTER_FILELIST.20201203.wrepbaseandtRNA.fa.fixed.fa.UpdatedSimpleRepeat.fa  
    - Identify transcript_type and positions of associated fasta header IDs from gencode.v33.chr_patch_hapl_scaff.annotation.gtf
    - Identify positions from associated IDs from hg38.repeatmasker.tsv.gz, hg38.simplerepeats.tsv.gz, hg38.trna.tsv.gz
    - Note any references that may exist in another database (eg. Refseq)
    - Note any associated IDs missing from these sources
    - If the number of missing IDs exceeds 1% of the total, REPORT and suggest possible sources.
*   Use pybedtools getfasta and the reference sequence (eg. GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta) to reproduce the original bowtie2 index fasta, sans any missing IDs. 
*   Generate a script that will generate a fasta file and bowtie2 index given gencode.{version}.annotation.gtf, {assembly}.repeatmasker.tsv.gz, {assembly}.simplerepeats.tsv.gz, {assembly}.trna.tsv.gz, {assembly}.fasta for mm10 and mm39.Script must also accept custom fasta files and positions.

**Step 3: Generate UniqueGenomicElements**
*   Infer the positions based on the provided reference data generated so far. Sources may be:
    - hg38.repeatmasker.tsv.gz (column 6)
    - hg38.trna.tsv.gz (column 6)
    - hg38.simplerepeats.tsv.gz (column 6)
    - gencode.v33.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat (should always be "-")
    - hsa.gff3 ("Name" attribute)
*   "-proximal" entries should be 500 bases upstream and downstream of their original reference. For example, the entry "MI0000112" at chrX:152394411-152394492 has two -proximal sites at: chrX:152393911-152394411 and chrX:152394492-152394992.
*   Write a python script that will regenerate UniqueGenomicElements.hg38.bed from provided references. Note that the .gff3 or trna.tsv.gz may not exist so it must be optional. Script must also accept custom fasta files and positions.
*   Run python script and check that UniqueGenomicElements.hg38.bed is faithfully reproduced. Then apply script to mm10 and mm39.

**Step 4: Generate MASTER_FILELIST**
*  Infer the columns based on the provided reference data generated so far. 
*  Write a python script that will faithfully reproduce MASTER_FILELIST.20201203.wrepbaseandtRNA.enst2id.fixed.UpdatedSimpleRepeat.wmiRs.list from 
*  Apply python script to generate an equivalent MASTER_FILELIST for mm10 and mm39. Script must also accept custom fasta files and positions.

**Step 5: Verify format and check perl scripts for compatibility**
*  Trace the CWL workflows cwl/wf_ecliprepmap_se.cwl and cwl/wf_ecliprepmap_pe.cwl and associated python or perl scripts located in bin/.
*  Check compatibility of the reproduced reference dataset for hg38. If incompatible, assess why and return to previous steps if necessary.
*  Check compatibility of both reference datasets for mm10 and mm39. If incompatible due to hardcoded values, modify perl scripts 

##### Output Requirements
Format your response as follows:
1. **All hg38 files reproduced:** 1 folder (bowtie_index) and three files (parsed_ucsc_tableformat, MASTER_FILELIST, and UniqueGenomicElements) must be regenerated and within 99% similarity to the existing hg38 references. 
2. **All mm10 and mm39 files generated:** there must be exactly 1 folder (bowtie_index) and three files (parsed_ucsc_tableformat, MASTER_FILELIST, and UniqueGenomicElements) generated.
3. **Notes:** Notes on potential missing references, position inconsistency, or hardcoded references to hg38 that must be modified before applying these mm39 and mm10 references.

##### Acceptance Checks
1. When the above references are generated for hg38, mm10, and mm39
2. When generated hg38 references pass 99% similarity checks to existing hg38 references
3. When the references pass identical format checks
4. When perl scripts pass compatibility checks