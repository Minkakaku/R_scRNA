docker run --name diopy_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u wzb \
wangzhenbo/monocle:20240909 Rscript /data/NAS03/backup/wzb/wzb/pipline/dio/dior_run.R -r $1 -o $2 -f $3 -s 'cell_type'