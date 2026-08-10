#!/usr/bin/bash

gtf="../Homo_sapiens.GRCh38.109.gtf"
script="./dexseq_prepare_annotation2.py"

conda run -n htseq python ${script} -f ./Homo_sapiens.GRCh38.109.flat.gtf ${gtf} ./Homo_sapiens.GRCh38.109.flat.gff
