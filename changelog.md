# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)

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

