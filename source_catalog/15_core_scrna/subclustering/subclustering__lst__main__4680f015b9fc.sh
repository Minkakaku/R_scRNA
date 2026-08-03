docker run --name proportion_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
wangzhenbo/monocle:20230719 Rscript /lvdata/wzb/pipline/proportion/cell_proportion.R \
'/lvdata/wzb/scRNA/FW2024-225_04/Recluster/NK/celltype/group/meta.xls' \
cell_type \
group \
'/lvdata/wzb/scRNA/FW2024-225_04/Recluster/NK/celltype/group' \
'/lvdata/wzb/scRNA/FW2024-225_04/order.xls' \
'/lvdata/wzb/scRNA/FW2024-225_04/Recluster/NK/celltype/group/cell_typecolor_dict.xls' 
