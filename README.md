# `gopher-pipelines` Recfator
This repository contains the Python refactoring of `gohper-pipelines`.

## TODO
These are some discrete tasks that need to be accomplished within major
divisions of the pipeline:

- **General Usability**
    - Finish bulk RNAseq argument parsing
- **Pipeline Object**
    - Check directories exist and can be written into
    - Merge user-specified options with default options
    - `BulkRNASeqPipeline`: check HISAT2 index exists and is complete
- **Samplesheet Object**
    - Parse UMGC samplesheet (And verify the consistency of the format)
    - Build samplesheet from directory
    - Build `gopher-pipelines` samplesheet input for `qsub`
    - Write samplesheet to disk
- **PBS Scripts**
- **Documentation**
    - Write manual draft
