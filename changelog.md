# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)

## [0.2.0] - 2026-05-21 Snakemake workflow

### Added
- Snakemake workflow (`Snakefile`, `workflow/rules/SE.smk`, `workflow/rules/PE.smk`,
  `workflow/rules/common.smk`) as an alternative to the CWL launchers
- Conda environment definition `workflow/envs/dropin.yaml` (python=3.10, bowtie2≥2.5,
  samtools≥1.17, numpy, pandas)
- SLURM profile `profiles/tscc2_snakemake9/` for running on TSCC with Snakemake 9
- Snakemake config examples for SE/PE small (downsampled) and SE/PE full datasets
- Downsampled example inputs in `examples/inputs/downsampled/`

### Changed
- README updated with step-by-step Snakemake quickstart instructions

### Removed
- Unused Perl scripts: `bin/perl/duplicate_removal.pl`, `bin/perl/RepElement_pipeline_1dataset.pl`

## [0.1.0] - 2020-12-13 last file before GRCh38


### Removed
- Removed extraneous scripts and CWL documents leftover from previous versions
- Added a Dockerfile definition to each tool

## [Unreleased 0.0.4b] - 2019-05-01
### Changed
- Changed the filenames a bit for a few tools.

## [Unreleased 0.0.4a] - 2019-03-26
### Fixed
- fixed barcode issue in SE pipeline that didn't properly grab inline barcodes

## [0.0.4] - 2019-02-14
### Added
- Extra workflow steps:
  - calculate_fold_change_from_parsed_files.cwl: tool for taking IP/Input parsed files and generating fold change and entropy scores
  - reparse_samfile_updatedchrM_fixmultenstsort_PE.cwl: re-parses the paired-end SAMlike file
  
- Single-end processing:
  - map_repetitive_elements_se.cwl: maps repetitive elements from single-end reads
  - deduplicate_se.cwl: deduplicates single-end reads
  - splitbam_se.cwl: splits a BAM or SAM file into 25 parts (5bases[ATCGN]^2)
  - reparse_samfile_updatedchrM_fixmultenstsort_SE.cwl: re-parses the single-end SAMlike file
  - wf_ecliprepmap_se_1sample.cwl: single-sample (either IP or Input usually) repeat-mapping workflow
  - wf_ecliprepmap_se.cwl: repeat-mapping workflow using IP and Input single samples
- 
### Changed
- wf_ecliprepmap.cwl -> wf_ecliprepmap_pe.cwl
- wf_ecliprepmapsingle.cwl -> wf_ecliprepmap_pe_1barcode.cwl
- deduplicate.cwl -> deduplicate_pe.cwl
- maprep.cwl -> map_repetitive_elements_pe.cwl
- splitbam.cwl -> splitbam_pe.cwl

### Deprecated
- 

## 0.0.3 - 2018-02-12
### Added
- First sharable commit to github
- Generates repeat-mapped rmDup and preRmDup SAM-like files mapping reads to repeat element families.
- README detailing methods and output explanations

