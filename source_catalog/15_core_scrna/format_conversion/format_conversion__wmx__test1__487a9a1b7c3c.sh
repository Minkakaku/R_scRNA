# 测试数据1
# rename 这个参数是取决于loom文件的barcode中包含的sample名是不是与--input_h5ad参数读入的文件中的sample对应的，如果不对应就需要给出对应关系
python Velocity.py \
  --input_h5ad /mnt/NexData_file/home/wzb/backup/scRNA/AE2512157482/2026.04.01/Malignant_cells/celltype/result.h5ad \
  --loom_dir /mnt/NexData_file/nas/outputdata/cellranger_out/AE2512157482/loom/ \
  --rename_file /mnt/NexData_file/nas/outputdata/cellranger_out/AE2512157482/rename.csv \
  --output_dir /mnt/NexData_file/home/wmx/Script/Velocity/test1/test_1/ \
  --sample_col sample \
  --groupby cell_type \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5 

# 不计算且绘制基因
# 加上  --no_rank_velocity_genes
python Velocity.py \
  --input_h5ad /mnt/NexData_file/home/wzb/backup/scRNA/AE2512157482/2026.04.01/Malignant_cells/celltype/result.h5ad \
  --loom_dir /mnt/NexData_file/nas/outputdata/cellranger_out/AE2512157482/loom/ \
  --rename_file /mnt/NexData_file/nas/outputdata/cellranger_out/AE2512157482/rename.csv \
  --output_dir /mnt/NexData_file/home/wmx/Script/Velocity/test1/test_1_nogene/ \
  --sample_col sample \
  --groupby cell_type \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5 \
  --no_rank_velocity_genes
 
 # 样品不对应报错，见log_1.txt
 python Velocity.py \
  --input_h5ad /mnt/NexData_file/home/wzb/backup/scRNA/AE2512157482/2026.04.01/Malignant_cells/celltype/result.h5ad \
  --loom_dir /mnt/NexData_file/nas/outputdata/cellranger_out/AE2512157482/loom/ \
  --output_dir /mnt/NexData_file/home/wmx/Script/Velocity/test1/test_1_error/ \
  --sample_col sample \
  --groupby cell_type \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5 \
  --no_rank_velocity_genes
   # --rename_file /mnt/NexData_file/nas/outputdata/cellranger_out/AE2512157482/rename.csv \