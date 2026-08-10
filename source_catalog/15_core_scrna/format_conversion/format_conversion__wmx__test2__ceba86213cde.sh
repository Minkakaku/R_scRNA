# 测试数据2 
# rename 这个参数是取决于loom文件的barcode中包含的sample名是不是与--input_h5ad参数读入的文件中的sample对应的，如果不对应就需要给出对应关系
python Velocity.py \
  --input_h5ad "/mnt/NexData_file/home/zwl/Project/FWSC20251318/subcelltype/Myeloid/celltype/anno.h5ad" \
  --loom_dir "/mnt/NexData_file/nas/outputdata/cellranger_out/FWSC20251318/" \
  --output_dir "/mnt/NexData_file/home/wmx/Script/Velocity/test2/test_2/" \
  --sample_col sample \
  --rename_file "/mnt/NexData_file/home/wmx/Script/Velocity/test2_rename.csv" \
  --groupby cell_type \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5 

# 不计算且绘制基因
# 加上  --no_rank_velocity_genes
python Velocity.py \
  --input_h5ad "/mnt/NexData_file/home/zwl/Project/FWSC20251318/subcelltype/Myeloid/celltype/anno.h5ad" \
  --loom_dir "/mnt/NexData_file/nas/outputdata/cellranger_out/FWSC20251318/" \
  --output_dir "/mnt/NexData_file/home/wmx/Script/Velocity/test2/test_2_nogene/" \
  --sample_col sample \
  --rename_file "/mnt/NexData_file/home/wmx/Script/Velocity/test2_rename.csv" \
  --groupby cell_type \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5 \
  --no_rank_velocity_genes

# 这里就是因为不对所以报错，详情见log_2.txt
python Velocity.py \
  --input_h5ad "/mnt/NexData_file/home/zwl/Project/FWSC20251318/subcelltype/Myeloid/celltype/anno.h5ad" \
  --loom_dir "/mnt/NexData_file/nas/outputdata/cellranger_out/FWSC20251318/" \
  --output_dir "/mnt/NexData_file/home/wmx/Script/Velocity/test2/test_2_error/" \
  --sample_col sample \
  --groupby cell_type \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5 
  # --rename_file "/mnt/NexData_file/home/wmx/Script/Velocity/test2_rename.csv" \


