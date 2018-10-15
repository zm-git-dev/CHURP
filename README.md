# `gopher-pipelines` Recfator
This repository contains the Python refactoring of `gohper-pipelines`.

## Testing
A test dataset is located at `/panfs/roc/scratch/konox006/gopher-pipelines_test`.
It contains tutorial RNAseq reads in `fastq/` and the hg19/GRCh37.p13 release
of the human genome in `genome/`.

The fastq files are simply the `Tutorial_File_R1.fastq` and
`Tutorial_File_R2.fastq` file copied over with new names. All R1 files are
identical and R2 files are identical. This shouldn't matter for testing
purposes.

The genome is the complete GRCh37 pre-built index from 
`ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch37_snp_tran.tar.gz`. It has
been built with SNP and splice site databases.

## TODO
These are some discrete tasks that need to be accomplished within major
divisions of the pipeline:

- **General Usability**
    - ~~Finish bulk RNAseq argument parsing~~ (Done)
    - ~~User-friendly error messages~~ (Done)
- **Pipeline Object**
    - ~~Check directories exist and can be written into~~ (Done)
    - ~~Merge user-specified options with default options~~ (Done)
    - ~~Write `-l` option to `qsub` based on user input for `threads` and `mem`.~~ (Done)
    - ~~Build `qsub` command (see `pipeline.sh`)~~ (Done)
    - ~~Add `--overwrite` option~~ (Done)
    - ~~Write program versions and date/time info into samplesheet data~~
    - ~~`BulkRNASeqPipeline`: check HISAT2 index exists and is complete~~ (Done)
    - `BulkRNASeqPipeline`: list available organisms in risdb
- **Samplesheet Object**
    - ~~Finalize data in samplesheet~~ (Done)
    - ~~Parse UMGC samplesheet (And verify the consistency of the format)~~ (Not doing)
    - ~~Build samplesheet from directory~~ (Done)
    - ~~Build `gopher-pipelines` samplesheet input for `qsub`~~ (Done)
    - ~~Write samplesheet to disk~~ (Done)
- **PBS Scripts**
    - ~~Add checkpointing for `resume`-like behavior~~ (Done)
    - ~~Choose software versions to use for `hisat2`, `samtools`, `R`, `trimmomatic`, and `fastqc`~~ (Done)
    - Add option for `overwrite`, with a default of false
    - `EdgeR` scripts for raw counts, TPM, FPKM, etc.
- **Documentation**
    - Write manual draft
    - Update diagrams
    - Change versioning to a three-part number (major.minor.patch; semantic versioning)
    - Make workflow diagram for each pipeline
    - Start a CHANGELOG
