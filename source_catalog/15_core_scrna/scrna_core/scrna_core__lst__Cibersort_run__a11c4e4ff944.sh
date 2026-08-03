docker run --name cib_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
wangzhenbo/diopy:20231121 Rscript /home/wzb/pipline/dio/dior_run.R -r $1 -o $2 -f $3 -s 'cell_type'