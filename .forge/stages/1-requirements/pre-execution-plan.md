## Pre-Execution Plan: 1-requirements

1. **Three most likely failure modes**:
   - Missing CWL step details: translate_cwl.md describes the pipeline at a high level; the architect needs exact CWL tool inputs/outputs. Mitigation: requirements stage must read all CWL files in cwl/ directly.
   - Ambiguous SE/PE config structure: translate_cwl.md says "enforce appropriate input structure" but doesn't define the config schema. Mitigation: surface this as an explicit AC for the architect.
   - Scatter/deduplication memory threshold unclear: translate_cwl.md says "run full non-deduplicated example to profile step_deduplicate" — this is a runtime decision, not a pre-implementation AC. Mitigation: flag as a deferred decision requiring profiling.

2. **First verification steps**: Confirm that the 11-14 context files cover: (a) exact CWL tool specs, (b) SE vs PE workflow differences, (c) test data locations, (d) expected output file formats and acceptance checks.

3. **Context dependencies**:
   - prompts/translate_cwl.md (the PRD)
   - cwl/ directory (all CWL files — the source to translate)
   - examples/repeat_mapping_PE.yaml and examples/repeat_mapping_SE.yaml (input format)
   - examples/example_data_for_repeat_mapping_hg38/ (test data)
   - test-provenance/ (reference outputs for validation)
   - profiles/tscc2_snakemake9/ (SLURM profile constraints)
