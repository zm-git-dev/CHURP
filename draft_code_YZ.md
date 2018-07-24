# New_GopherPipeline

1. Written in Python
2. Use JobArray to submit jobs

## Pseudo Code
```
def parse_argument:
    Required:
        FastqFolder
        HiSat2index
        Annotationfile
    Optional:
        Trimmomatic
        OutputDir
        SingleCell (a separate function?)

def setup_workdir:
    Input:
        FastqFolder    
    Output:
        Create /scratch.global/username-pipeline directory
        Create a samplesheet
        Create symbolic links

def submit_jobs_per_sample:
    Input:
        SampleSheet
        JobTemplate
    Output:
        subprocess.call(qsub_command, shell=True)    

def summarize_count:
    Input:
        BAM files
    Output:
        Count Table / FPKM Table / Other R ready dataframe

def generate_log:
    Output:    
        Check progress
        Report Error
        SetUp resume operation

def main():
    parse_argument
    setup_workdir
    submit_jobs_per_sample
    summarize_count
    generate_log
```

## Job Template
- Standard PBS setting: #PBS -l nodes=1:ppn=8,mem=16gb,walltime=12:00:00
- Command line arguments: SampelID (JobArrayID), HiSat2index, AnnotationFile, WorkDIR
- Optional Command line arguments
```
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=12:00:00
cd ${WorkDIR}/${SampleID}
module load fastqc hisat2 picard
fastqc -t 2 R1.fastq R2.fastq
if TRIM: 
    trimmomatic
    fastqc 
hisat2 
picard
```

## Required modules in Python 3.x:
- OS/subprocess
- ArgParse
- Glob
- Error Checking (?)
