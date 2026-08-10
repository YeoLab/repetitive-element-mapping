## Pre-Execution Plan: 6-security

1. **Three most likely failure modes**:
   - subprocess shell injection: map_repetitive_elements_*.py calls bowtie2 with user-supplied paths. Risk if shell=True with f-strings.
   - Path traversal: scripts open files from sys.argv or config without sanitization. Low severity for CLI tools but worth documenting.
   - SAM parsing robustness: malformed SAM lines (wrong field count, non-integer CIGAR) could raise unhandled exceptions. Not a security issue per se but a reliability concern.

2. **First verification steps**: Check subprocess calls for shell=True; check all open() calls for user-supplied paths; check sys.argv handling.

3. **Context dependencies**:
   - workflow/scripts/map_repetitive_elements_pe.py (subprocess bowtie2 call)
   - workflow/scripts/map_repetitive_elements_se.py (subprocess bowtie2 call)
   - workflow/scripts/split_bam_to_subfiles.py (file I/O)
   - workflow/scripts/merge_parsed_files.py (file I/O)
   - workflow/scripts/deduplicate.py (most complex, GTF/BED parsing)
