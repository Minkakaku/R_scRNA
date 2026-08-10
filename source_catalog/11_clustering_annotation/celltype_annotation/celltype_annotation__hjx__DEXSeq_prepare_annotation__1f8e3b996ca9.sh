#!/usr/bin/bash

gtf="../Mus_musculus.GRCm39.109.gtf"
script="./dexseq_prepare_annotation2.py"

conda run -n htseq python ${script} -f ./Mus_musculus.GRCm39.109.flat.gtf ${gtf} ./Mus_musculus.GRCm39.109.flat.gff
