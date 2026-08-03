# python /mnt/NexData_file/home/wmx/Script/Velocity/Velocity.py \
#   --input_h5ad "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/CD8T_checked.h5ad" \
#   --rename_file "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/rename.csv" \
#   --loom_dir "/mnt/NexData_file/nas/outputdata/cellranger_out/SC2507033447/2026.06.11" \
#   --output_dir "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Velocity/" \
#   --sample_col sample \
#   --groupby cell_type \
#   --add_sample_prefix \
#   --basis umap \
#   --prefix scvelo \
#   --top_gene_n 5

python /mnt/NexData_file/home/wmx/Script/Velocity/Velocity.py \
  --input_h5ad "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/CD8T_checked.h5ad" \
  --rename_file "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/rename.csv" \
  --loom_dir "/mnt/NexData_file/nas/outputdata/cellranger_out/SC2507033447/2026.06.11" \
  --output_dir "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Velocity_Cluster/" \
  --sample_col sample \
  --groupby Cluster \
  --add_sample_prefix \
  --basis umap \
  --prefix scvelo \
  --top_gene_n 5
