docker run --name Roe_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u wzb \
wangzhenbo/seurat_pip:v1.1 Rscript /lvdata/wzb/pipline/Roe/roe.R