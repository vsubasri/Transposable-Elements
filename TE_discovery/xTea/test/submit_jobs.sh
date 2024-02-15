#!/bin/bash

qsub -l vmem=48g,mem=48g,walltime=48:00:00 /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/test/NA12878/L1/run_xTEA_pipeline.sh
qsub -l vmem=48g,mem=48g,walltime=48:00:00 /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/test/NA12878/Alu/run_xTEA_pipeline.sh
qsub -l vmem=48g,mem=48g,walltime=48:00:00 /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/test/NA12878/SVA/run_xTEA_pipeline.sh
qsub -l vmem=48g,mem=48g,walltime=48:00:00 /hpf/largeprojects/davidm/vsubasri/transposable_elements/xTea/test/NA12878/HERV/run_xTEA_pipeline.sh
