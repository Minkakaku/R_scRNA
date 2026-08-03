library(data.table)

mouse_anno <- fread(file.choose())
bacteria_anno <- fread(file.choose())
total_annotation <- rbindlist(list(mouse_anno, bacteria_anno), fill = TRUE)

mouse_produ <- fread(file.choose())
bacteria_produ <- fread(file.choose())
total_production <- rbindlist(list(mouse_produ, bacteria_produ), fill = TRUE)

output_filename_anno <- "annotation.tsv"
output_filename_produ <- "product.tsv"


data.table::fwrite(total_annotation, 
                   file = output_filename_anno, 
                   sep = "\t", 
                   quote = FALSE, 
                   row.names = FALSE)

data.table::fwrite(total_production, 
                   file = output_filename_produ, 
                   sep = "\t", 
                   quote = FALSE, 
                   row.names = FALSE)