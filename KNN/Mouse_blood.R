source("mippie_nc.R")

mippie_subnet <- mippie_nc(query.file = "gene_symbol_Mouse_embryonic_blood.txt",
path.to.mippie = "mippie_ppi_v1_0.tsv",
path.to.proteins = "mippie_proteins_v1_0.tsv", order = 1, output.file = "mipple_subset_blood.tsv");
