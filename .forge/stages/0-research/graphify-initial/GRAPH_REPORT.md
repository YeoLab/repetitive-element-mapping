# Graph Report - .  (2026-05-14)

## Corpus Check
- Corpus is ~3,617 words - fits in a single context window. You may not need a graph.

## Summary
- 71 nodes · 75 edges · 13 communities (7 shown, 6 thin omitted)
- Extraction: 92% EXTRACTED · 8% INFERRED · 0% AMBIGUOUS · INFERRED: 6 edges (avg confidence: 0.86)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_CWL Pipeline Steps|CWL Pipeline Steps]]
- [[_COMMUNITY_Fold Change Metrics|Fold Change Metrics]]
- [[_COMMUNITY_Pipeline Orchestration|Pipeline Orchestration]]
- [[_COMMUNITY_Parsed Output Format|Parsed Output Format]]
- [[_COMMUNITY_Reference Data Generation|Reference Data Generation]]
- [[_COMMUNITY_Deduplication & Conflict Resolution|Deduplication & Conflict Resolution]]
- [[_COMMUNITY_Fold Change Script (AST)|Fold Change Script (AST)]]
- [[_COMMUNITY_Lab Meeting Slides (SAM Output)|Lab Meeting Slides (SAM Output)]]
- [[_COMMUNITY_SE Bowtie2 Parser|SE Bowtie2 Parser]]
- [[_COMMUNITY_Docker Image|Docker Image]]
- [[_COMMUNITY_Initial Release|Initial Release]]
- [[_COMMUNITY_Lab Slides (YAML Usage)|Lab Slides (YAML Usage)]]
- [[_COMMUNITY_Lab Slides (Output Structure)|Lab Slides (Output Structure)]]

## God Nodes (most connected - your core abstractions)
1. `CWL Workflow (production)` - 9 edges
2. `eCLIP Repetitive Element Mapping Pipeline` - 7 edges
3. `return_l2fc_entropy_from_parsed function` - 6 edges
4. `main function (calculate_fold_change)` - 5 edges
5. `eCLIP Repetitive Element Mapping Pipeline (ecliprepmap)` - 5 edges
6. `Snakemake Dropin Workflow` - 5 edges
7. `return_l2fc_entropy_from_parsed()` - 4 edges
8. `Mouse Assembly mm10 Reference Data` - 4 edges
9. `Mouse Assembly mm39 Reference Data` - 4 edges
10. `read_parsed()` - 3 edges

## Surprising Connections (you probably didn't know these)
- `YAML Template Usage Pattern (TEMPLATE.ecliprepmap)` --references--> `CWL Workflow (production)`  [INFERRED]
  lab_meeting_slides_rep_element_pipeline.pdf → CLAUDE.md
- `Pipeline Motivation v2: CLIP read mapping to repetitive elements` --references--> `eCLIP Repetitive Element Mapping Pipeline`  [INFERRED]
  lab_meeting_slides_rep_element_pipeline-20171018.pdf → README.md
- `eCLIP Repetitive Element Mapping Pipeline (ecliprepmap)` --references--> `main function (calculate_fold_change)`  [EXTRACTED]
  CLAUDE.md → python/calculate_fold_change_from_parsed_files.py
- `calculate_fold_change_from_parsed_files.cwl (added v0.0.4)` --references--> `main function (calculate_fold_change)`  [EXTRACTED]
  changelog.md → python/calculate_fold_change_from_parsed_files.py
- `return_l2fc_entropy_from_parsed function` --references--> `Multi-family Reads (pipe-separated ambiguous mappings)`  [EXTRACTED]
  python/calculate_fold_change_from_parsed_files.py → CLAUDE.md

## Hyperedges (group relationships)
- **CWL Per-barcode Sub-workflow Steps** — claude_md_map_repetitive_elements_cwl, claude_md_splitbam_cwl, claude_md_getpair_cwl, claude_md_deduplicate_cwl, claude_md_combine_cwl [EXTRACTED 1.00]
- **Reference Data Generation Workflow (mm10/mm39/hg38)** — generate_refdata_parsed_ucsc_tableformat, generate_refdata_bowtie2_index_gen, generate_refdata_unique_genomic_elements, generate_refdata_master_filelist, generate_refdata_perl_compat_check [EXTRACTED 1.00]
- **Fold Change Computation Pipeline** — calculate_fold_change_read_parsed, calculate_fold_change_return_l2fc_entropy, calculate_fold_change_fold_enrichment, calculate_fold_change_information_content, calculate_fold_change_pseudocount [EXTRACTED 1.00]
- **Snakemake Python Shim Layer over Perl Scripts** — claude_md_perl_compat, claude_md_split_bam_py, claude_md_merge_parsed_py, claude_md_parse_bowtie2_pe_pl, claude_md_parse_bowtie2_se_pl [EXTRACTED 0.95]
- **Supported Genome Assemblies** — generate_refdata_hg38_reference, generate_refdata_mm10_assembly, generate_refdata_mm39_assembly [EXTRACTED 1.00]

## Communities (13 total, 6 thin omitted)

### Community 0 - "CWL Pipeline Steps"
Cohesion: 0.18
Nodes (11): Single-end (SE) Pipeline Addition (v0.0.4), combine CWL step (merge_multiple_parsed_files), CWL Workflow (production), getpair CWL step, map_repetitive_elements CWL step, parse_bowtie2_output_realtime_includemultifamily_PE.pl, rRNA Special Handling (RNA45S precursor), splitbam CWL step (+3 more)

### Community 1 - "Fold Change Metrics"
Cohesion: 0.24
Nodes (10): Fold Enrichment Metric, Information Content Metric, main function (calculate_fold_change), nopipes TSV Output (unambiguous mappings), Pseudocount for Missing Input Reads, return_l2fc_entropy_from_parsed function, withpipes TSV Output (all mappings including multifamily), calculate_fold_change_from_parsed_files.cwl (added v0.0.4) (+2 more)

### Community 2 - "Pipeline Orchestration"
Cohesion: 0.22
Nodes (10): Bowtie2 Repeat Element Database, eCLIP Repetitive Element Mapping Pipeline (ecliprepmap), merge_multiple_parsed_files.simplified_20191022.py (pure Python port), _perl_compat.py (Perl shim resolver), se_foundation.smk (SE foundation tests), Snakemake Dropin Workflow, split_bam_to_subfiles_SEorPE.py (pure Python port), Bowtie2 Mapping Step (+2 more)

### Community 3 - "Parsed Output Format"
Cohesion: 0.22
Nodes (10): Parsed File Format (.parsed), read_parsed function, eCLIP Repetitive Element Mapping Pipeline, Parsed Output File Format, Van Nostrand et al. 2020 (Genome Biology, eCLIP 150 RBPs), Pipeline Motivation v2: CLIP read mapping to repetitive elements, Rep Element Pipeline Lab Meeting Overview v2 (2017-10-18), Pipeline Motivation: Integrate Eric's Analysis into eCLIP (+2 more)

### Community 4 - "Reference Data Generation"
Cohesion: 0.31
Nodes (9): Version 0.1.0 GRCh38 Release, Bowtie2 Index Generation (Step 2), Human Assembly hg38 Reference Data, MASTER_FILELIST Generation (Step 4), Mouse Assembly mm10 Reference Data, Mouse Assembly mm39 Reference Data, parsed_ucsc_tableformat Generation (Step 1), 500bp Proximal Region Flanking Strategy (+1 more)

### Community 5 - "Deduplication & Conflict Resolution"
Cohesion: 0.25
Nodes (8): Repeat vs Unique Genome Conflict Resolution, deduplicate CWL step, duplicate_removal_inline_paired...pl, Perl Version Hash Iteration Non-determinism, Rationale: 24-unit alignment score threshold for conflict resolution, Rationale: UMI splitting for memory efficiency, UMI-based 25-bin Splitting Strategy, UMI-based Deduplication Step

### Community 6 - "Fold Change Script (AST)"
Cohesion: 0.47
Nodes (5): main(), From 2 parsed rep element pipeline outputs (ip and input),     compute fold chan, Reads Eric's parsed file from the repetitive element pipeline.     Parameters, read_parsed(), return_l2fc_entropy_from_parsed()

## Knowledge Gaps
- **29 isolated node(s):** `Reads Eric's parsed file from the repetitive element pipeline.     Parameters`, `From 2 parsed rep element pipeline outputs (ip and input),     compute fold chan`, `Fold Enrichment Metric`, `Pseudocount for Missing Input Reads`, `splitbam CWL step` (+24 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **6 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `eCLIP Repetitive Element Mapping Pipeline (ecliprepmap)` connect `Pipeline Orchestration` to `CWL Pipeline Steps`, `Fold Change Metrics`, `Parsed Output Format`?**
  _High betweenness centrality (0.296) - this node is a cross-community bridge._
- **Why does `eCLIP Repetitive Element Mapping Pipeline` connect `Parsed Output Format` to `Pipeline Orchestration`, `Deduplication & Conflict Resolution`?**
  _High betweenness centrality (0.213) - this node is a cross-community bridge._
- **Why does `CWL Workflow (production)` connect `CWL Pipeline Steps` to `Pipeline Orchestration`, `Deduplication & Conflict Resolution`?**
  _High betweenness centrality (0.156) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `CWL Workflow (production)` (e.g. with `Single-end (SE) Pipeline Addition (v0.0.4)` and `YAML Template Usage Pattern (TEMPLATE.ecliprepmap)`) actually correct?**
  _`CWL Workflow (production)` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Reads Eric's parsed file from the repetitive element pipeline.     Parameters`, `From 2 parsed rep element pipeline outputs (ip and input),     compute fold chan`, `Fold Enrichment Metric` to the rest of the system?**
  _29 weakly-connected nodes found - possible documentation gaps or missing edges._