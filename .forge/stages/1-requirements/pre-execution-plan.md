## Pre-Execution Plan: 1-requirements

1. **Three most likely failure modes**:
   - Missing mm10/mm39 source files: generate_refdata.md assumes files exist under examples/inputs/{assembly}/downloaded/. If absent, requirements will be unverifiable. Watch for empty directory listings.
   - Format spec derivation errors: MASTER_FILELIST and UniqueGenomicElements formats must be inferred from hg38 reference; incorrect column mapping creates broken AC. Watch for column count mismatches when examining hg38 reference files.
   - Scope creep on Perl modifications: prompt says "modify only if hardcoded values prevent acceptance of new annotations." Must nail down exactly which values are hardcoded and what constitutes an acceptable threshold.

2. **First verification steps**:
   - List examples/inputs/ to confirm all three assemblies exist (hg38, mm10, mm39)
   - Inventory downloaded/ subdirectories for each assembly
   - Inspect hg38 reference output files to extract column schemas

3. **Context dependencies**:
   - prompts/generate_refdata.md (primary task spec)
   - examples/inputs/hg38/ (format ground truth)
   - examples/inputs/mm10/downloaded/ and examples/inputs/mm39/downloaded/ (source file availability)
   - bin/perl/*.pl (for scope-of-modification assessment)
   - cwl/ workflow files (for compatibility verification context)
   - .forge/stages/0-research/graphify-initial/GRAPH_REPORT.md (knowledge graph context)
