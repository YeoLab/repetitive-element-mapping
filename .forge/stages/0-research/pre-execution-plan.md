## Pre-Execution Plan: 0-research

1. **Three most likely failure modes**:
   - Missing downloaded source files: The prompt assumes mm10/mm39 source files exist under examples/inputs/{assembly}/downloaded. If they are absent, all downstream steps stall. Watch for file listing returning empty.
   - Format drift between assemblies: GTF annotation formats for mm10 (GRCm38) and mm39 (GRCm39) may have schema differences from hg38 (GRCh38) that break the parsing logic. Watch for column count mismatches in the GTF.
   - Module load failures: `module load ecliprepmap/1.0.0` may not be available in the current shell context or may have different PATH exports than expected. Watch for command-not-found errors after module load.

2. **First verification steps**:
   - List examples/inputs/ to confirm mm10, mm39, hg38 directories exist
   - Check examples/inputs/{assembly}/downloaded/ for source files
   - Confirm examples/inputs/hg38/ contains the four reference outputs (bowtie2_index, parsed_ucsc_tableformat, UniqueGenomicElements, MASTER_FILELIST)

3. **Context dependencies**:
   - prompts/generate_refdata.md (concept document — already read)
   - examples/inputs/hg38/ reference outputs (format ground truth)
   - bin/ perl and python scripts (compatibility check in Step 5)
   - cwl/ workflow definitions (trace in Step 5)
   - CLAUDE.md (project conventions)
