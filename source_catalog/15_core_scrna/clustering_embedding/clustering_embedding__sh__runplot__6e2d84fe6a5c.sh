python "/home/sh/sc_pipline/Marker/Markerplot.py" --input_dir "/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/Myeloid cells/result.h5ad" --marker_file "/home/sh/sc_pipline/Marker/human-Myeloid_population.txt" --species human --output_dir "/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/Myeloid cells/Markerplot/" --groupby 'leiden'

python "/mnt/NexData_file/home/lst/scrna-code/scanpy_fun/plot_celltype_proportion.py" \
--h5ad "/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/Myeloid cells/Myeloid.h5ad" \
--outdir "/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/" \
--group_key group