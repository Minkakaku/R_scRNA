docker run --name diopy_sh--rm \
-v /home/sh:/home/sh \
-v /data/ :/data/ \
-u sh \
sh/monocle:20241226 Rscript /data/NAS02/sh/Projects/FWSC20241287/01.Mission/script/h5ad2rds/dior_run.R -r $1 -o $2 -f $3 -s 'cell_type'