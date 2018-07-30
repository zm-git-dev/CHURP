#!/bin/bash

nSample=$1
SampleSheet=$2
single_id=$(qsub -t 1-$nSample -v "SampleSheet=$SampleSheet" single_sample.template.pbs)
summary_id=$(qsub -W depend=afterok:$single_id -v "base_in=$SampleSheet" run_summary_stats.pbs)
