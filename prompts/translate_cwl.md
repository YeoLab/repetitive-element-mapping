#### System Prompt: Translate CWL scripts to Snakemake

**Role:** Act as an expert bioinformatician.

### **Skills**
Use the skills provided in:
- .claude/skills/

##### Context
*   **Tech Stack:** Common Workflow Language (CWL), Snakemake
*   **Architecture:** Workflows written in Common Workflow Language (CWL) and Snakemake
*   **Background:** We are performing a direct translation of a CWL-based workflow to Snakemake. There are two distinct workflows, one workflow that processes paired-end (PE) reads, and the other that processes single-end (SE) reads. Briefly described, both workflows map fastq files to repeat elements, and their mapped coordinates are combined with a BAM file to produce SAM-like and TSV files with coordinates mapping each read to either a repeat or unique genomic element. These files from each barcode (IP) are then compared against a corresponding (input) background to identify enrichment. Each PE dataset is comprised of three groups 1. "barcode1", "barcode2" and "input". Each group is comprised of R1 fastq.gz, R2 fastq.gz, and 1 BAM file mapped to the genome. Each SE dataset is comprised of two groups 1. "barcode1" and "input". Each group is comprised of R1 fastq, and 1 BAM file mapped to the genome. For both SE and PE datasets, all reads from each BAM file must exist in their corresponding fastq.gz file. 

##### Scope
*   You may **read and write** from/to files in this folder root /tscc/nfs/home/bay001/projects/codebase/repetitive-element-mapping/
*   You may **read** from all files and softlinks in this folder root /tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/ **DO NOT WRITE**
*   If at any point the original pipeline fails, STOP. Ask me what to do next.
*   Do NOT make any additional changes to the perl scripts themselves. If it is absolutely required (script CANNOT be used with snakemake), ask me.

##### Constraints & Framework Rules
*   **Always check each python script against the original :** Translated script outputs should be identical or near-identical (due to random tie-breaking or perl/python-specific rounding)
*   **Use the example data:** Example fastq.gz and BAM files for both SE and PE datasets are provided inside examples/example_data_for_repeat_mapping_hg38
*   **Use the reference data:** Pipeline reference files are provided inside examples/inputs/hg38/


##### Logic & Implementation Steps 
**Step 1: Generate small, testable dataset**
*  Generate a downsampled paired-end dataset from EXAMPLE_PE* files in examples/example_data_for_repeat_mapping_hg38
*  Generate a downsampled single-end dataset from EXAMPLE_SE* files in examples/example_data_for_repeat_mapping_hg38
*  Downsampled files should go into: examples/inputs/downsampled/
*  Datasets must include all chromosomes (eg. chr1, chr2, ... chrM).
*  All reads within downsampled BAM files must exist in their corresponding downsampled fastq file. See context background.
*  Datasets must include at least 100 reads with barcodes prefixed (first two bases) by all base combinations. Include "CG", "CA", "CT", "CC", ... "NN" reads.
   - For "SE" reads. For example: the read "K00180:223:HCLHCBBXX:5:1101:10013:5587_CGCCTTGCCG" has a barcode with the prefix "CG". 
   - For "PE" reads. For example: The read "AGAAA:SN1001:449:HGTN3ADXX:1:1101:12732:17919" has a barcode with the prefix "AG".  
*  Generate small config.yaml files for these small datasets. The original config.yaml files exist here: examples/*.yaml
   - repeat_mapping_PE.yaml -> repeat_mapping_PE_small.yaml
   - repeat_mapping_SE.yaml -> repeat_mapping_SE_small.yaml
   - chmod +x both files so they use /tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/eCLIP_repelement_SE_singleNode and /tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/eCLIP_repelement_PE_singleNode as the executor. 
   
**Step 2: Run current pipeline on test data to produce expected outputs**
*   Use: "module load ecliprepmap perl/5.10.1" and ensure the following environment is loaded:
    - perl: /tscc/projects/ps-yeolab4/software/perl/5.10.1/bin/perl
    - python: /tscc/projects/ps-yeolab4/software/miniconda_tscc2/envs/ecliprepmap-0.1.0/bin/python
    - eCLIP_repelement_PE_singleNode: /tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/eCLIP_repelement_PE_singleNode
    - eCLIP_repelement_PE_singleNode: /tscc/projects/ps-yeolab4/software/ecliprepmap/0.1.0/wf/eCLIP_repelement_SE_singleNode
*   Execute the YAML directly (eg. ./repeat_mapping_SE_small.yaml and ./repeat_mapping_PE_small.yaml)
*   Obtain all outputs for further testing.

**Step 3: Translate CWL workflows to Snakemake workflows**
*   Use: "conda activate snakemake9 && module load bowtie2/2.2.6 python3essential"
    - Ensure Snakemake version == 9.12.0. Use this snakemake to run all tests
    - bowtie2 and python may be required to run certain scripts.
*   Generate one Snakefile.
*   Generate a config file that allows the user to specify either PE or SE, and enforce appropriate input structure. That is, if the user specifies PE, he must provide both R1 and R2 fastqs.
*   Within the Snakefile, add logic to determine which subworkflow (either PE or SE) will be executed.
*   Modularize each workflow rules according to either "SE.smk", "PE.smk" or "common.smk".
*   Translate each step directly to Snakemake format. 
*   Use the "requirements" field to set memory (ramMin) and CPU (coresMin) requirements. Assume an 8:00:00 walltime, 1 CPU and 32G memory by default.
*   Small tests should be lightweight and not exceed default requirements. Set a retry limit to 2 and adjust based on retry attempts (eg. "mem_mb = lambda wildcards, attempt: attempt * 2000")
*   Run the full, non-deduplicated example to profile step_deduplicate and estimate memory requirements
    - If the full example does not exceed 32G, then we do not need to scatter the step_deduplicate jobs. Remove this logic and re-test.
    - Else, keep the scatter logic inside the Snakemake workflow.

**Step 4: Test new Snakemake workflow and ensure outputs match expected outputs from Step 2**
*  Ensure outputs match the original expected
*  Perform full run on original, non-downsampled data using both CWL and Snakemake workflows, then test for drift
   - Use profiles/tscc2_snakemake9
   - Adjust memory, CPU and walltime requirements if necessary, for the full datasets

**Step 6: Cleanup**
*  COMMIT CHANGES TO GIT AND PUSH BEFORE DELETING ANYTHING.
*  Update README and changelog to reflect changes and usage
   - README.md should include step-by-step instructions for deploying this snakemake workflow using small, downsampled data.
*  Remove unused perl scripts
*  Remove unused CWL definitions 

##### Output Requirements
Format your response as follows:
1. **Summary of files added:**
2. **Summary of files translated:**  
3. **Summary of files removed:** 
4. **Notes:**
   
##### Acceptance Checks
1. When all downsampled datasets are produced
2. When the Snakemake workflow runs to completion on both downsampled and full datasets
3. When both downsampled and full dataset Snakemake outputs match original CWL outputs. **This is THE MOST IMPORTANT acceptance criteria.**
4. When all code is documented
5. When all changes are documented