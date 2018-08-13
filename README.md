# `gopher-pipelines` Recfator
This repository contains the Python refactoring of `gohper-pipelines`.

## TODO
These are some discrete tasks that need to be accomplished within major
divisions of the pipeline:

- **General Usability**
    - ~~Finish bulk RNAseq argument parsing~~ (Done)
- **Pipeline Object**
    - Check directories exist and can be written into
    - ~~Merge user-specified options with default options~~ (Done)
    - Write `-l` option to `qsub` based on user input for `threads` and `mem`.
    - Build `qsub` command (see `pipeline.sh`)
    - Add `--overwrite` option
    - Write program versions and date/time info into samplesheet data
    - `BulkRNASeqPipeline`: check HISAT2 index exists and is complete
    - `BulkRNASeqPipeline`: list available organisms in risdb
- **Samplesheet Object**
    - ~~Finalize data in samplesheet~~ (Done)
    - Parse UMGC samplesheet (And verify the consistency of the format)
    - Build samplesheet from directory
    - Build `gopher-pipelines` samplesheet input for `qsub`
    - Write samplesheet to disk
- **PBS Scripts**
    - ~~Add checkpointing for `resume`-like behavior~~ (Done)
    - Add option for `overwrite`, with a default of false
    - `EdgeR` scripts for raw counts, TPM, FPKM, etc.
- **Documentation**
    - Write manual draft
