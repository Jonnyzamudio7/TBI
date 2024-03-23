library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")

merged_metagenomes <- import_biom("~/GIT/TBI/new_otu_table.biom")
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


percentages <- transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )


library(phyloseq)

# Suponiendo que `percentages` es tu objeto Phyloseq ya creado

# Aglomerar por Filo
percentages_phylum <- tax_glom(percentages, taxrank = 'Phylum')
# Convertir a dataframe
phylum_df_percentages <- psmelt(percentages_phylum)
# Exportar a CSV
write.csv(phylum_df_percentages, 'phylum_df_percentages.csv', row.names = FALSE)

# Aglomerar por Orden
percentages_orden <- tax_glom(percentages, taxrank = 'Order')
# Convertir a dataframe
order_df_percentages <- psmelt(percentages_orden)
# Exportar a CSV
write.csv(order_df_percentages, 'order_df_percentages.csv', row.names = FALSE)

# Aglomerar por Género
percentages_genus <- tax_glom(percentages, taxrank = 'Genus')
# Convertir a dataframe
genus_df_percentages <- psmelt(percentages_genus)
# Exportar a CSV
write.csv(genus_df_percentages, 'genus_df_percentages.csv', row.names = FALSE)
