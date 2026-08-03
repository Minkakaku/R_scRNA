#Rscript scTenifoldKnk.R -i All_filter_state_samp_new_with_CNs10-15_N_mcaf.rds -o ./N_MMP14_results --max_cells 2000 --nComp 3 --fdr_threshold 0.1 --min_cells 5 -g MMP14 && \
#Rscript scTenifoldKnk.R -i All_filter_state_samp_new_with_CNs10-15_A_mcaf.rds -o ./A_MMP14_results --max_cells 2000 --nComp 3 --fdr_threshold 0.1 --min_cells 5 -g MMP14 && \
#Rscript scTenifoldKnk.R -i All_filter_state_samp_new_with_CNs10-15_C_mcaf.rds -o ./C_MMP14_results --max_cells 2000 --nComp 3 --fdr_threshold 0.1 --min_cells 5 -g MMP14 && \
#Rscript scTenifoldKnk.R -i All_filter_state_samp_new_with_CNs10-15_N_mcaf.rds -o ./N_THY1_results --max_cells 2000 --nComp 3 --fdr_threshold 0.1 --min_cells 5 -g THY1 && \
#Rscript scTenifoldKnk.R -i All_filter_state_samp_new_with_CNs10-15_A_mcaf.rds -o ./A_THY1_results --max_cells 2000 --nComp 3 --fdr_threshold 0.1 --min_cells 5 -g THY1 && \
Rscript scTenifoldKnk.R -i All_filter_state_samp_new_with_CNs10-15_C_mcaf.rds -o ./C_THY1_results --max_cells 2000 --nComp 3 --fdr_threshold 0.1 --min_cells 5 -g THY1 
